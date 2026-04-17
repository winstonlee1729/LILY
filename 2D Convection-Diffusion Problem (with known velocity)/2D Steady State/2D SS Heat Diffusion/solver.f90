program heat_diffusion
    use omp_lib
    implicit none
      
! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny
    integer*4 maxiter, iter
    integer*4 i, j
    parameter(nx = 100, ny = 100)
    parameter(maxiter = 10000)
    
    ! Flow parameter
    real*8, parameter :: Lx = 1.0D0, Ly = 1.0D0
    real*8 dx, dy

    ! Flow field
    real*8, allocatable :: T(:, :), Told(:, :)
    
    ! Program convergence parameters
    real*8 tol, maxres, err
    parameter(tol = 1.0D-5)
    
! Defining variable
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
      
! Allocate memory
    allocate(T(nx, ny), Told(nx, ny))
    
! Initialization
    do i = 1, nx
        do j = 1, ny
            T(i, j) = 0.0D0
        enddo
    enddo
    ! Applying boundary condition
    do i = 1, nx
        T(i, 1) = 100.0D0
    enddo
      
! Gauss-Seidel Iteration
    do iter = 1, maxiter
        maxres = 0.0D0
          
        ! Save the old temperature field to calculate the residual
        Told = T
        
        !$omp parallel do collapse(2) reduction(max:maxres) private(err)
        do i = 2, nx-1
            do j = 2, ny-1
                T(i, j) = 0.25D0 * (T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1))
                  
                err = abs (T(i, j) - Told(i, j))
                if (err >= maxres) maxres = err
                  
            enddo
        enddo
        !$omp end parallel do
        
        if (mod(iter, 100) == 0) then
            write(*,*) 'Iteration: ', iter, '|Max Residual = ', maxres
        endif
          
        if (maxres <= tol) then
            write(*,*) 'Converged at iteration ', iter
            exit
        endif
          
    enddo

! Open the output file
    open(unit=10, file='flow_field.txt', status='unknown')
      
! Outputting solution
    do i = 1, nx
        do j = 1, ny
            write(10, '(3E15.6)') (dble(i) - 0.5D0) * dx, (dble(j) - 0.5D0) * dy, T(i, j)
        enddo
    enddo
      
! Close the file
    close(10)
      
    stop
end program heat_diffusion