program CD_QUICK
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
    real*8, allocatable :: phi(:, :), phiold(:, :), tmp(:, :)
    real*8 phi_new
    
    ! Parameters for convergence
    real*8 tol, maxres, err, rf
    parameter(tol = 1.0D-5)
    parameter(rf = 0.6D0)
    
    ! Flow parameters
    real*8 u, v, rho 
    real*8 Lx, Ly, dx, dy, thickness
    real*8 area_w, area_e, area_n, area_s
    real*8 ga_w, ga_e, ga_n, ga_s
    real*8 De, Dw, Dn, Ds 
    real*8 Fe, Fw, Fn, Fs, delta_F
    parameter(Lx = 1.0D0, Ly = 1.0D0, thickness = 1.0D0)
    parameter(u = 0.0D0, v = 0.1D0)
    parameter(rho = 1.0D0)
    parameter(ga_w = 1.0D0, ga_e = 1.0D0, ga_n = 1.0D0, ga_s = 1.0D0)
    
    ! Equation coefficients
    real*8 alf_w, alf_e, alf_s, alf_n
    real*8 aw_up, ae_up, as_up, an_up, ap_up
    real*8 aw_qk, ae_qk, as_qk, an_qk, ap_qk
    real*8 aww_qk, aee_qk, ass_qk, ann_qk
    
! Allocate memory
    allocate(phi(nx, ny), phiold(nx, ny), tmp(nx,ny))
      
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
      
    ! For alpha value
    ! For alpha w
    if (Fw > 0) then
        alf_w = 1.0D0
    else
        alf_w = 0.0D0
    endif
    ! For alpha e
    if (Fe > 0) then
        alf_e = 1.0D0
    else
        alf_e = 0.0D0
    endif
    ! For alpha s
    if (Fs > 0) then
        alf_s = 1.0D0
    else
        alf_s = 0.0D0
    endif
    ! For alpha n
    if (Fn > 0) then
        alf_n = 1.0D0
    else
        alf_n = 0.0D0
    endif
    
    ! Upwind scheme coefficient
    aw_up = max(Fw, Dw + (Fw / 2.0D0), 0.0D0)
    ae_up = max(-Fe, De - (Fe / 2.0D0), 0.0D0)
    as_up = max(Fs, Ds + (Fs / 2.0D0), 0.0D0)
    an_up = max(-Fn, Dn - (Fn / 2.0D0), 0.0D0)
    ap_up = aw_up + ae_up + as_up + an_up + delta_F
      
    ! QUICK scheme coefficient 
    aw_qk = Dw + (6.0D0/8.0D0)*alf_w*Fw + (1.0D0/8.0D0)*alf_e*Fe + (3.0D0/8.0D0)*(1.0D0-alf_w)*Fw
    aww_qk = (-1.0D0/8.0D0)*alf_w*Fw
    
    ae_qk = De - (3.0D0/8.0D0)*alf_e*Fe - (6.0D0/8.0D0)*(1.0D0-alf_e)*Fe - (1.0D0/8.0D0)*(1.0D0-alf_w)*Fw
    aee_qk = (1.0D0/8.0D0)*(1.0D0-alf_e)*Fe

    as_qk = Ds + (6.0D0/8.0D0)*alf_s*Fs + (1.0D0/8.0D0)*alf_n*Fn + (3.0D0/8.0D0)*(1.0D0-alf_s)*Fs
    ass_qk = (-1.0D0/8.0D0)*alf_s*Fs
    
    an_qk = Dn - (3.0D0/8.0D0)*alf_n*Fn - (6.0D0/8.0D0)*(1.0D0-alf_n)*Fn - (1.0D0/8.0D0)*(1.0D0-alf_s)*Fs
    ann_qk = (1.0D0/8.0D0)*(1.0D0-alf_n)*Fn

    ap_qk = aw_qk + ae_qk + as_qk + an_qk + aww_qk + aee_qk + ass_qk + ann_qk + delta_F
      
! Initialization
    !$omp parallel do collapse(2)
    do i = 1, nx
        do j = 1, ny
            phi(i, j) = 0.0D0
            phiold(i, j) = 0.0D0
        enddo
    enddo
    !$omp end parallel do
    ! Applying boundary condition
    do i = 1, nx
        phi(i, 1) = 100.0D0
        phiold(i, 1) = 100.0D0
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
        

        ! Pass 1: Red cells
        !$omp parallel do private(i, j, err, phi_new) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                if (mod(i+j, 2) == 0) then
                      
                    if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1) then
                        phi_new = (aw_up*phi(i-1, j) + ae_up*phi(i+1, j) + as_up*phi(i, j-1) + an_up*phi(i, j+1)) / ap_up
                    else
                        phi_new = (aw_qk*phi(i-1, j) + ae_qk*phi(i+1, j) + as_qk*phi(i, j-1) + an_qk*phi(i, j+1) &
                        + aww_qk*phi(i-2, j) + aee_qk*phi(i+2, j) + ass_qk*phi(i, j-2) + ann_qk*phi(i, j+2)) / ap_qk
                    endif
                      
                    phi(i, j) = (1.0D0 - rf) * phiold(i, j) + rf * phi_new
                      
                    err = abs(phi(i, j) - phiold(i, j))
                    if (err > maxres) maxres = err    
                endif
            enddo
        enddo
        !$omp end parallel do


        ! Pass 2: Black cells
        !$omp parallel do private(i, j, err, phi_new) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                if (mod(i+j, 2) /= 0) then
                      
                    if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1) then
                        phi_new = (aw_up*phi(i-1, j) + ae_up*phi(i+1, j) + as_up*phi(i, j-1) + an_up*phi(i, j+1)) / ap_up
                    else
                        phi_new = (aw_qk*phi(i-1, j) + ae_qk*phi(i+1, j) + as_qk*phi(i, j-1) + an_qk*phi(i, j+1) &
                        + aww_qk*phi(i-2, j) + aee_qk*phi(i+2, j) + ass_qk*phi(i, j-2) + ann_qk*phi(i, j+2)) / ap_qk
                    endif
                      
                    phi(i, j) = (1.0D0 - rf) * phiold(i, j) + rf * phi_new
                      
                    err = abs(phi(i, j) - phiold(i, j))
                    if (err > maxres) maxres = err    
                endif
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
end program CD_QUICK