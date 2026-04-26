program Transient_CD_Hybrid
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny, nz
    integer*4 i, j, k
    integer*4 time_step, max_step
    integer*4 iter, maxiter, iter_limit
    integer*4 output_freq
    parameter(nx = 100, ny = 100, nz = 100)
    parameter(maxiter = 5000, max_step = 1000)
    parameter(output_freq = 25)
    
    ! Program convergence parameter
    real*8 tol, maxres, err, rhs, lhs

    ! Flow field
    real*8, allocatable :: phi(:,:,:), phiold(:,:,:), phi_new(:,:,:)
    
    ! Flow parameters
    real*8 rho, u, v, w
    real*8 area_w, area_e, area_n, area_s, area_t, area_b
    real*8 ga_w, ga_e, ga_n, ga_s, ga_t, ga_b
    real*8 dx, dy, dz, Lx, Ly, Lz, volume, delta_t
    real*8 Fe, Fw, Fn, Fs, Ft, Fb, delta_F
    real*8 De, Dw, Dn, Ds, Db, Dt
    real*8 aw, ae, an, as, at, ab, ap, ap_o, ap_prime
    real*8 w_node, e_node, n_node, s_node, t_node, b_node, p_past
    parameter(tol = 1.0D-7)
    parameter(rho = 1.0D0, u = 0.0D0, v = 0.1D0, w = 0.0D0)
    parameter(ga_w = 0.01D0, ga_e = 0.01D0, ga_n = 0.01D0, ga_s = 0.01D0, ga_b = 0.01D0, ga_t = 0.01D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0)
    parameter(delta_t = 0.05D0)
    
    ! Output string for VTK files
    character(len=50) :: vtk_filename
    
    ! Temporal discretization scheme
    real*8 theta
    ! 1.0 = Fully Implicit | 0.5 = Crank-Nicolson | 0.0 = Explcit
    parameter(theta = 1.0D0)
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    dz = Lz / dble(nz)
    
    area_w = dy * dz
    area_e = dy * dz
    area_n = dx * dz
    area_s = dx * dz
    area_t = dx * dy
    area_b = dx * dy
    
    volume = dx * dy * dz
    
    Fw = rho * u * area_w
    Fe = rho * u * area_e
    Fn = rho * v * area_n
    Fs = rho * v * area_s
    Ft = rho * w * area_t
    Fb = rho * w * area_b
    delta_F = Fe - Fw + Fn - Fs + Ft - Fb
    
    De = (ga_e / dx) * area_e
    Dw = (ga_w / dx) * area_w
    Dn = (ga_n / dy) * area_n
    Ds = (ga_s / dy) * area_s
    Dt = (ga_t / dz) * area_t
    Db = (ga_b / dz) * area_b
    
    ap_o = rho * volume / delta_t
    aw = max(Fw, Dw + (Fw/2), 0.0D0)
    ae = max(-Fe, De - (Fe/2), 0.0D0)
    as = max(Fs, Ds + (Fs/2), 0.0D0)
    an = max(-Fn, Dn - (Fn/2), 0.0D0)
    ab = max(Fb, Db + (Fb/2), 0.0D0)
    at = max(-Ft, Dt - (Ft/2), 0.0D0)
    
    ap = ap_o + (ae + aw + as + an + at + ab + delta_F) * theta
    ap_prime = ap_o - (ae + aw + an + as + at + ab + delta_F) * (1.0D0 - theta)
    
! Allocate memory 
    allocate(phi(nx, ny, nz), phiold(nx, ny, nz), phi_new(nx, ny, nz))
    
! Initialize flow field
    phi = 0.0D0
    phiold = 0.0D0
    phi_new = 0.0D0
    ! Applying boundary condition
    !$omp parallel do collapse(2)
    do i = 1, nx
        do k = 1, nz
            phi(i, 1, k) = 100.0D0
            phiold(i, 1, k) = 100.0D0
        enddo
    enddo
    !$omp end parallel do
    
! Check for inner iteration number
    if (theta == 0.0D0) then
        iter_limit = 1
    else
        iter_limit = maxiter
    endif
    
! Start time step
    do time_step = 1, max_step
        
        ! Save old field
        !$omp parallel do collapse(3)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    phiold(i, j, k) = phi(i, j, k)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        ! Inner iteration 
        do iter = 1, iter_limit
            maxres = 0.0D0
            
            !$omp parallel do collapse(3) private(i, j, k, e_node, w_node, n_node, s_node, t_node, b_node, &
            !$omp p_past, rhs, lhs, err) reduction(max:maxres)
            do i = 2, nx-1
                do j = 2, ny-1
                    do k = 2, nz-1
                        w_node = aw * (theta * phi(i-1, j, k) + (1.0D0 - theta) * phiold(i-1, j, k))
                        e_node = ae * (theta * phi(i+1 ,j, k) + (1.0D0 - theta) * phiold(i+1, j, k))
                        s_node = as * (theta * phi(i, j-1, k) + (1.0D0 - theta) * phiold(i, j-1, k))
                        n_node = an * (theta * phi(i, j+1, k) + (1.0D0 - theta) * phiold(i, j+1, k))
                        b_node = ab * (theta * phi(i, j, k-1) + (1.0D0 - theta) * phiold(i, j, k-1))
                        t_node = at * (theta * phi(i, j, k+1) + (1.0D0 - theta) * phiold(i, j, k+1))
                        p_past = ap_prime * phiold(i, j, k)
                    
                        rhs = w_node + e_node + s_node + n_node + t_node + b_node + p_past
                        lhs = ap * phi(i, j, k)
                    
                        err = abs(rhs - lhs)
                        if (err >= maxres) maxres = err
                    
                        phi_new(i, j, k) = rhs / ap
                    enddo
                enddo
            enddo
            !$omp end parallel do
            
            !$omp parallel do collapse(3) private(i, j, k)
            do i = 2, nx-1
                do j = 2, ny-1
                    do k = 2, nz-1
                        phi(i, j, k) = phi_new(i, j, k)
                    enddo
                enddo
            enddo
            !$omp end parallel do
            
            if (theta > 0.0D0 .and. maxres <= tol) exit
        enddo
        
        if (theta == 0.0D0) then
            write(*,*) 'Time step:', time_step
        else
            write(*,*) 'Time step:', time_step, '|Iters:', iter, '|Res:', maxres
        endif
        
        ! Outputting solution 
        if (mod(time_step, output_freq) == 0) then
            write(vtk_filename, '(A, I4.4, A)') 'flow_', time_step, '.vtk'
            call write_vtk(vtk_filename, phi, nx, ny, nz, dx, dy, dz)
        endif
        
    enddo
    

! Deallocate memory
    deallocate(phi, phiold, phi_new)
    
    stop
end program Transient_CD_Hybrid


! This subroutine aims at outputting matrices as Legacy VTK format
subroutine write_vtk(filename, phi, nx, ny, nz, dx, dy, dz)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny, nz
    real*8, intent(in) :: dx, dy, dz
    real*8, intent(in) :: phi(nx, ny, nz)
    integer*4 :: i, j, k
    
    ! Open file
    open(unit=20, file=filename, status='replace')
    
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') 'CFD Transient Output'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET STRUCTURED_POINTS'
    
    ! Grid properties
    write(20, '(A, 3I8)') 'DIMENSIONS ', nx, ny, nz
    write(20, '(A, 3F10.6)') 'ORIGIN ', 0.0, 0.0, 0.0
    write(20, '(A, 3F10.6)') 'SPACING ', dx, dy, dz
    
    ! Data assignment 
    write(20, '(A, I12)') 'POINT_DATA ', nx * ny * nz
    write(20, '(A)') 'SCALARS phi float 1'
    write(20, '(A)') 'LOOKUP_TABLE default'
    
    ! Write phi values 
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                write(20, '(E16.6)') phi(i, j, k)
            enddo
        enddo
    enddo
    
    close(20)
end subroutine write_vtk