program CD_Upwind
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny, nz
    integer*4 iter, maxiter
    integer*4 i, j, k
    parameter(nx = 100, ny = 100, nz = 100)
    parameter(maxiter = 100000)
    
    ! Calculate program time
    real*8 start_time, end_time
    
    ! Flow field
    real*8, allocatable :: phi(:,:,:), phiold(:,:,:)
    
    ! Parameters for convergence
    real*8 maxres, err
    real*8, parameter :: tol = 1.0D-5
    
    ! Flow parameters
    real*8 u, v, w, rho
    real*8 Lx, Ly, Lz, dx, dy, dz
    real*8 area_w, area_e, area_n, area_s, area_t, area_b
    real*8 ga_w, ga_e, ga_n, ga_s, ga_t, ga_b
    parameter(Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0)
    parameter(u = 0.0D0, v = 0.1D0, w = 0.0D0)
    parameter(rho = 1.0D0)
    parameter(ga_w = 1.0D0, ga_e = 1.0D0, ga_n = 1.0D0, ga_s = 1.0D0, ga_t = 1.0D0, ga_b = 1.0D0)
    
    ! Equation coefficients
    real*8 De, Dw, Dn, Ds, Dt, Db
    real*8 Fe, Fw, Fn, Fs, Ft, Fb, delta_F
    real*8 ap, aw, ae, an, as, at, ab
    real*8 x_term, y_term, z_term
    real*8 xc, yc, zc
    
! Allocate memory
    allocate(phi(nx, ny, nz), phiold(nx, ny, nz))
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    dz = Lz / dble(nz)
    
    area_e = dy * dz
    area_w = dy * dz
    area_n = dx * dz
    area_s = dx * dz
    area_b = dx * dy
    area_t = dx * dy
    
    Fe = rho * u * area_e
    Fw = rho * u * area_w
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
    
    aw = Dw + max(Fw, 0.0D0)
    ae = De + max(-Fe, 0.0D0)
    as = Ds + max(Fs, 0.0D0)
    an = Dn + max(-Fn, 0.0D0)
    ab = Db + max(Fb, 0.0D0)
    at = Dt + max(-Ft, 0.0D0)
    
    ap = aw + ae + as + an + ab + at + delta_F
    
! Initialization
    phi = 0.0D0
    phiold = 0.0D0
    ! Applying boundary condition
    !$omp parallel do collapse(2)
    do i = 1, nx
        do k = 1, nz
            phi(i, 1, k) = 100.0D0       
            phiold(i, 1, k) = 100.0D0
        enddo
    enddo
    !$omp end parallel do
    
! Start the clock
    start_time = omp_get_wtime()
    
! Jacobi Iteration
    do iter = 1, maxiter
        maxres = 0.0D0
        
        ! Save old flow field
        !$omp parallel do collapse(3)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    phiold(i, j, k) = phi(i, j, k)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        !$omp parallel do collapse(3) reduction(max:maxres) private(err, x_term, y_term, z_term)
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    x_term = aw * phiold(i-1, j, k) + ae * phiold(i+1, j, k)
                    y_term = an * phiold(i, j+1, k) + as * phiold(i, j-1, k)
                    z_term = at * phiold(i, j, k+1) + ab * phiold(i, j, k-1)
                    
                    phi(i, j, k) = (x_term + y_term + z_term) / ap
                    
                    err = abs(phi(i, j, k) - phiold(i, j, k))
                    maxres = max(maxres, err)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        !if (mod(iter, 100) == 0) then
            write(*,*) 'Iteration: ', iter, '| Max residual:', maxres
        !endif
        
        if (maxres <= tol) exit
        
    enddo
    
! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total iterations: ', iter
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'
    
! Open the output file
    open(unit=10, file='flow_field.txt', status='replace')
    
! Outputting solution
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                xc = (dble(i) - 0.5D0) * dx
                yc = (dble(j) - 0.5D0) * dy
                zc = (dble(k) - 0.5D0) * dz
                write(10, '(4E16.6)') xc, yc, zc, phi(i, j, k)
            enddo
        enddo
    enddo
    
! Close the output file
    close(10)
    
! Deallocate memory
    deallocate(phi, phiold)
    
    stop
end program CD_Upwind
        
        
        
    
    
    
    
    
    
    
    
    