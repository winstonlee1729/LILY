program CD_Hybrid
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny 
    integer*4 iter, maxiter 
    integer*4 i, j
    parameter(nx = 100, ny = 100, maxiter = 1000000)
    
    ! Calculate program time
    real*8 start_time, end_time
    
    ! Flow field
    real*8, allocatable :: phi(:, :), phiold(:, :)
    
    ! Parameters for convergence
    real*8 tol, maxres, err
    parameter(tol = 1.0D-5)
    
    ! Flow parameters
    real*8 u, v, rho 
    real*8 Lx, Ly, dx, dy, thickness
    real*8 area_w, area_e, area_n, area_s
    real*8 ga_w, ga_e, ga_n, ga_s
    parameter(Lx = 1.0D0, Ly = 1.0D0, thickness = 1.0D0)
    parameter(u = 0.0D0, v = 0.1D0)
    parameter(rho = 1.0D0)
    parameter(ga_w = 1.0D0, ga_e = 1.0D0, ga_n = 1.0D0, ga_s = 1.0D0)
    
    ! Equation coefficients
    real*8 De, Dw, Dn, Ds 
    real*8 Fe, Fw, Fn, Fs, delta_F
    real*8 ap, aw, ae, an, as
    
! Allocate memory
    allocate(phi(nx, ny), phiold(nx, ny))
      
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    area_e = dy * thickness
    area_w = dy * thickness
    area_n = dx * thickness
    area_s = dx * thickness
    
    Fe = rho * u * area_e
    Fw = rho * u * area_w
    Fn = rho * v * area_n
    Fs = rho * v * area_s
      
    De = (ga_e / dx) * area_e
    Dw = (ga_w / dx) * area_w
    Dn = (ga_n / dy) * area_n
    Ds = (ga_s / dy) * area_s
    delta_F = Fe - Fw + Fn - Fs
      
    aw = Dw + max(Fw, 0.0D0)
    ae = De + max(-Fe, 0.0D0)
    as = Ds + max(Fs, 0.0D0)
    an = Dn + max(-Fn, 0.0D0)
    ap = aw + ae + as + an + delta_F
      
      
! Initialization
    !$omp parallel do collapse(2)
    do i = 1, nx
        do j = 1, ny
            phi(i, j) = 0.0D0
        enddo
    enddo
    !$omp end parallel do
    ! Applying boundary condition
    do i = 1, nx
        phi(i, 1) = 100.0D0
    enddo
      
! Start the clock
    start_time = omp_get_wtime()
      
! Gauss-Seidel Iteration
    do iter = 1, maxiter
        maxres = 0.0D0
          
        ! Save old flow field to calculate the residual
        !$omp parallel do collapse(2)
        do i = 1, nx
            do j = 1, ny
                phiold(i, j) = phi(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        !$omp parallel do collapse(2) reduction(max:maxres) private(err)
        do i = 2, nx-1
            do j = 2, ny-1
                phi(i,j) = (aw * phi(i-1, j) + ae * phi(i+1, j) + an * phi(i, j+1) + as * phi(i, j-1)) / ap
                  
                err = abs(phi(i,j) - phiold(i,j))
                if (err >= maxres) maxres = err
            enddo
        enddo
        !$omp end parallel do
        
        if (mod(iter, 100) == 0) then
            write(*,*) 'Iteration: ', iter, '|Max Residual: ', maxres
        endif
          
          
        if (maxres <= tol) exit
          
    enddo
      
! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total Iterations: ', iter
    write(*,*) 'Total Time (s): ', end_time - start_time
    write(*,*) '------------------------------------'
      
! Open the output file
    open(unit=10, file='flow_field.txt', status='unknown')
      
! Outputting solution
do i = 1, nx
    do j = 1, ny
        write(10, '(3E15.6)') (dble(i) - 0.5D0) * dx, (dble(j) - 0.5D0) * dy, phi(i, j)
    enddo
enddo
      
! Close the output file
    close(10)
      
      
    stop
end program CD_Hybrid