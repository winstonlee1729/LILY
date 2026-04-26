program Transient_Diffusion
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration numbers
    integer*4 nx, ny, nz
    integer*4 i, j, k
    integer*4 time_step, max_step
    integer*4 iter, maxiter
    integer*4 output_freq
    parameter(nx = 100, ny = 100, nz = 100)
    parameter(max_step = 3000, maxiter = 5000)
    parameter(output_freq = 50) 
    
    ! Parameters for convergence
    real*8 maxres, err, rhs, lhs
    real*8, parameter :: tol = 1.0D-7
    
    ! Flow fields 
    real*8, allocatable :: phi(:,:,:), phiold(:,:,:), temp(:,:,:)
    
    ! Flow parameters
    real*8 ga_w, ga_e, ga_n, ga_s, ga_t, ga_b
    real*8 area_w, area_e, area_n, area_s, area_t, area_b
    real*8 dx, dy, dz, Lx, Ly, Lz, dt
    real*8 rho, Cp, volume, rad_sq, cx, cy, cz, R
    parameter(ga_e = 5.0D0, ga_w = 5.0D0, ga_n = 5.0D0, ga_s = 5.0D0, ga_t = 5.0D0, ga_b = 5.0D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0)
    parameter(dt = 0.05D0)
    parameter(rho = 1000.0D0, Cp = 1000.0D0)
    parameter(cx = 0.5D0, cy = 0.5D0, cz = 0.5D0, R = 0.05D0)
    
    ! Output string for VTK files
    character(len=50) :: vtk_filename
    
    ! Equation coefficients
    real*8 ap, ae, aw, an, as, at, ab, ap_o
    real*8 w_node, e_node, n_node, s_node, t_node, b_node, p_past
    
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
    
    ap_o = rho * Cp * (volume / dt)
    aw = (ga_w * area_w) / dx
    ae = (ga_e * area_e) / dx
    an = (ga_n * area_n) / dy
    as = (ga_s * area_s) / dy
    at = (ga_t * area_t) / dz
    ab = (ga_b * area_b) / dz
    ap = aw + ae + an + as + at + ab + ap_o
    
! Allocate memory
    allocate(phi(nx, ny, nz), phiold(nx, ny, nz), temp(nx, ny, nz))
    
! Initialization
    phi = 0.0D0
    phiold = 0.0D0
    temp = 0.0D0
    ! Initial condition
    !$omp parallel do collapse(3) private(i, j, k, rad_sq)
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                rad_sq = (dble(i) * dx - cx)**2 + (dble(j) * dy - cy)**2 + (dble(k) * dz - cz)**2
                if (rad_sq <= R**2) then
                    phi(i, j, k) = 100.0D0
                    phiold(i, j, k) = 100.0D0
                    temp(i, j, k) = 100.0D0
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do
    
! Start time loop
    do time_step = 1, max_step
        
        ! Save old field
        !$omp parallel do collapse(3) private(i, j, k)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    phiold(i, j, k) = phi(i, j, k)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        ! Inner iteration
        do iter = 1, maxiter
            maxres = 0.0D0
            
            !$omp parallel do collapse(3) private(i, j, k, e_node, w_node, n_node, s_node, t_node, b_node, &
            !$omp p_past, rhs, lhs, err) reduction(max:maxres)
            do i = 2, nx-1
                do j = 2, ny-1
                    do k = 2, nz-1
                        e_node = ae * phi(i+1, j, k)
                        w_node = aw * phi(i-1, j, k)
                        n_node = an * phi(i, j+1, k)
                        s_node = as * phi(i, j-1, k)
                        t_node = at * phi(i, j, k+1)
                        b_node = ab * phi(i, j, k-1)
                        p_past = ap_o * phiold(i, j, k)
                        
                        rhs = e_node + w_node + n_node + s_node + t_node + b_node + p_past
                        lhs = ap * phi(i, j, k)
                        
                        err = abs(rhs - lhs)
                        if(err >= maxres) maxres = err
                        
                        temp(i, j, k) = rhs / ap
                    enddo
                enddo
            enddo
            !$omp end parallel do
            
            ! Update main field
            !$omp parallel do collapse(3) private(i, j, k)
            do i = 2, nx-1
                do j = 2, ny-1
                    do k = 2, nz-1
                        phi(i, j, k) = temp(i, j, k)
                    enddo
                enddo
            enddo
            !$omp end parallel do
            
            if (maxres <= tol) exit
        enddo
        
        write(*,*) 'Time step:', time_step, '|Iters:', iter, '|Res:', maxres
        
        ! Outputting solution to VTK
        if (mod(time_step, output_freq) == 0) then
            write(vtk_filename, '(A, I4.4, A)') 'diff_', time_step, '.vtk'
            call write_vtk(trim(vtk_filename), phi, nx, ny, nz, dx, dy, dz) 
        endif
        
    enddo
    
! Deallocate memory
    deallocate(phi, phiold, temp)
    
    stop
end program Transient_Diffusion


! This subroutine aims at generating Legacy VTK files
subroutine write_vtk(filename, phi, nx, ny, nz, dx, dy, dz)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny, nz
    real*8, intent(in) :: dx, dy, dz
    real*8, intent(in) :: phi(nx, ny, nz)
    integer*4 :: i, j, k
    
    ! Open file
    open(unit=20, file=filename, status='replace')
    
    ! Header
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') '3D Transient Diffusion Output'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET STRUCTURED_POINTS'
    
    ! Dimensions
    write(20, '(A, I0, 1X, I0, 1X, I0)') 'DIMENSIONS ', nx, ny, nz
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ', 0.0, 0.0, 0.0
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'SPACING ', dx, dy, dz
    
    ! Data assignment 
    write(20, '(A, I0)') 'POINT_DATA ', nx * ny * nz
    write(20, '(A)') 'SCALARS phi float 1'
    write(20, '(A)') 'LOOKUP_TABLE default'
    
    ! Write phi values 
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                write(20, '(F12.6)') phi(i, j, k)
            enddo
        enddo
    enddo
    
    close(20)
end subroutine write_vtk