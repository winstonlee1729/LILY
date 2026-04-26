program heat_diffusion
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny, nz
    integer*4 i, j, k
    integer*4 iter, maxiter
    parameter(nx = 100, ny = 100, nz = 100)
    parameter(maxiter = 10000000)
    
    ! Flow parameter
    real*8, parameter :: Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0
    real*8 dx, dy, dz
    
    ! Flow field
    real*8, allocatable :: T(:,:,:), Told(:,:,:)
    
    ! Program running time
    real*8 start_time, end_time
    
    ! Parameters for program convergence
    real*8 maxres, err
    real*8, parameter :: tol = 1.0D-5
    
    ! Equation coefficients
    real*8 x_term, y_term, z_term
    
    ! For outputting solution
    real*8 xc, yc, zc
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    dz = Lz / dble(nz)
    
! Allocate memory
    allocate(T(nx, ny, nz), Told(nx, ny, nz))
    
! Initialization
    T = 0.0D0
    Told = 0.0D0
    !$omp parallel do collapse(2)
    do i = 1, nx
        do k = 1, nz
            T(i, 1, k) = 100.0D0       
            Told(i, 1, k) = 100.0D0
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
                    Told(i, j, k) = T(i, j, k)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        !$omp parallel do collapse(3) reduction(max:maxres) private(err, x_term, y_term, z_term)
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    x_term = Told(i-1, j, k) + Told(i+1, j, k)
                    y_term = Told(i, j+1, k) + Told(i, j-1, k)
                    z_term = Told(i, j, k+1) + Told(i, j, k-1)
                    
                    T(i, j, k) = (1.0D0 / 6.0D0) * (x_term + y_term + z_term)
                    
                    err = abs(T(i, j, k) - Told(i, j, k))
                    maxres = max(maxres, err)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        
        write(*,*) 'Iteration: ', iter, '| Max residual:', maxres
        
        if (maxres <= tol) exit
        
    enddo
    
! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total iterations: ', iter
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'
    
! Open the output file
    open(unit=10, file='temperature_field.txt', status='replace')
    
! Outputting solution
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                xc = (dble(i) - 0.5D0) * dx
                yc = (dble(j) - 0.5D0) * dy
                zc = (dble(k) - 0.5D0) * dz
                write(10, '(4E16.6)') xc, yc, zc, T(i, j, k)
            enddo
        enddo
    enddo
    
! Close the output file
    close(10)
    
! Deallocate memory
    deallocate(T, Told)
    
    stop
end program heat_diffusion
        
        
        
    
    
    
    
    
    
    
    
    
    