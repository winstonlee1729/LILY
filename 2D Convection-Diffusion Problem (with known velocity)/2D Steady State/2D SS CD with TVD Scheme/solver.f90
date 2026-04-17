program CD_TVD
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
    parameter(Lx = 1.0D0, Ly = 1.0D0, thickness = 1.0D0)
    parameter(u = 0.0D0, v = 0.1D0)
    parameter(rho = 1.0D0)
    parameter(ga_w = 1.0D0, ga_e = 1.0D0, ga_n = 1.0D0, ga_s = 1.0D0)
    
    ! Equation coefficients
    real*8 De, Dw, Dn, Ds 
    real*8 Fe, Fw, Fn, Fs, delta_F
    real*8 alf_w, alf_e, alf_s, alf_n
    real*8 aw_up, ae_up, as_up, an_up, ap_up
    real*8 aw_t, ae_t, as_t, an_t, ap_t, Su_DC
    real*8 rw_p, re_p, rs_p, rn_p
    real*8 rw_m, re_m, rs_m, rn_m
    real*8 psi_w, psi_e, psi_s, psi_n, tmp
    real*8 term1, term2, term3, term4, C
    parameter(C = 1.0D-15)

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

    ! Flow direction 
    if (Fw > 0) then; alf_w = 1.0D0; else; alf_w = 0.0D0; endif
    if (Fe > 0) then; alf_e = 1.0D0; else; alf_e = 0.0D0; endif
    if (Fs > 0) then; alf_s = 1.0D0; else; alf_s = 0.0D0; endif
    if (Fn > 0) then; alf_n = 1.0D0; else; alf_n = 0.0D0; endif
      
    ! Upwind scheme coefficients
    aw_up = max(Fw, Dw + (Fw / 2.0D0), 0.0D0)
    ae_up = max(-Fe, De - (Fe / 2.0D0), 0.0D0)
    as_up = max(Fs, Ds + (Fs / 2.0D0), 0.0D0)
    an_up = max(-Fn, Dn - (Fn / 2.0D0), 0.0D0)
    ap_up = aw_up + ae_up + as_up + an_up + delta_F
      
    ! TVD coefficients 
    aw_t = Dw + max(Fw, 0.0D0)
    ae_t = De + max(-Fe, 0.0D0)
    as_t = Ds + max(Fs, 0.0D0)
    an_t = Dn + max(-Fn, 0.0D0)
    ap_t = aw_t + ae_t + as_t + an_t + delta_F
      
      
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
        !$omp parallel do private(i, j, err, phi_new, re_p, rw_p, rs_p, rn_p &
        !$omp ,re_m, rw_m, rs_m, rn_m, tmp, psi_w, psi_e, psi_s, psi_n, su_DC) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                if (mod(i+j, 2) == 0) then
                      
                    if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1) then
                        phi_new = (aw_up*phi(i, j-1) + ae_up*phi(i, j+1) + as_up*phi(i+1, j) + an_up*phi(i-1, j)) / ap_up
                    else
                        re_p = (phi(i, j) - phi(i-1, j)) / (phi(i+1, j) - phi(i, j) + C)
                        rw_p = (phi(i-1, j) - phi(i-2, j)) / (phi(i, j) - phi(i-1, j) + C)
                        rs_p = (phi(i, j-1) - phi(i, j-2)) / (phi(i, j) - phi(i, j-1) + C)
                        rn_p = (phi(i, j) - phi(i, j-1)) / (phi(i, j+1) - phi(i, j) + C)
                          
                        re_m = (phi(i+2, j) - phi(i+1, j)) / (phi(i+1, j) - phi(i, j) + C)
                        rw_m = (phi(i+1, j) - phi(i, j)) / (phi(i, j) - phi(i-1, j) + C)
                        rs_m = (phi(i, j+1) - phi(i, j)) / (phi(i, j) - phi(i, j-1) + C)
                        rn_m = (phi(i, j+2) - phi(i, j+1)) / (phi(i, j-1) - phi(i, j) + C)
                         
                        ! Use UMIST model for limiter function
                        if (Fw > 0.0D0) then
                            tmp = min(2.0D0*rw_p, (1.0D0+3.0D0*rw_p)/4.0D0, (3.0D0+rw_p)/4.0D0, 2.0D0)
                            psi_w = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rw_m, (1.0D0+3.0D0*rw_m)/4.0D0, (3.0D0+rw_m)/4.0D0, 2.0D0)
                            psi_w = max(0.0D0, tmp)
                        endif
                          
                        if (Fe > 0.0D0) then
                            tmp = min(2.0D0*re_p, (1.0D0+3.0D0*re_p)/4.0D0, (3.0D0+re_p)/4.0D0, 2.0D0)
                            psi_e = max(0.0D0, tmp) 
                        else
                            tmp = min(2.0D0*re_m, (1.0D0+3.0D0*re_m)/4.0D0, (3.0D0+re_m)/4.0D0, 2.0D0)
                            psi_e = max(0.0D0, tmp) 
                        endif
                          
                        if (Fs > 0.0D0) then
                            tmp = min(2.0D0*rs_p, (1.0D0+3.0D0*rs_p)/4.0D0, (3.0D0+rs_p)/4.0D0, 2.0D0)
                            psi_s = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rs_m, (1.0D0+3.0D0*rs_m)/4.0D0, (3.0D0+rs_m)/4.0D0, 2.0D0)
                            psi_s = max(0.0D0, tmp)
                        endif

                        if (Fn > 0.0D0) then
                            tmp = min(2.0D0*rn_p, (1.0D0+3.0D0*rn_p)/4.0D0, (3.0D0+rn_p)/4.0D0, 2.0D0)
                            psi_n = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rn_m, (1.0D0+3.0D0*rn_m)/4.0D0, (3.0D0+rn_m)/4.0D0, 2.0D0)
                            psi_n = max(0.0D0, tmp)
                        endif
                          
                        tmp = (1.0D0 - alf_e) * psi_e - alf_e * psi_e
                        term1 = 0.5D0 * Fe * (phi(i+1, j) - phi(i, j)) * tmp
                          
                        tmp = alf_w * psi_w - (1.0D0 - alf_w) * psi_w
                        term2 = 0.5D0 * Fw * (phi(i, j) - phi(i-1, j)) * tmp
                      
                        tmp = (1.0D0 - alf_n) * psi_n - alf_n * psi_n
                        term3 = 0.5D0 * Fn * (phi(i, j+1) - phi(i, j)) * tmp
                          
                        tmp = alf_s * psi_s - (1.0D0 - alf_s) * psi_s
                        term4 = 0.5D0 * Fs * (phi(i, j) - phi(i, j-1)) * tmp
                          
                        Su_DC = term1 + term2 + term3 + term4
                          
                        phi_new = (aw_t * phi(i-1, j) + ae_t * + phi(i+1, j) + as_t * phi(i, j-1) &
                                + an_t * phi(i, j+1) + Su_DC) / ap_t
                    endif
                      
                    phi(i,j) = (1.0D0 - rf) * phiold(i, j) + rf * phi_new
                      
                    err = abs(phi(i, j) - phiold(i, j))
                    if (err > maxres) maxres = err
                      
                endif
            enddo
        enddo
        !$omp end parallel do
        
        ! Pass 2: Black cells 
        !$omp parallel do private(i, j, err, phi_new, re_p, rw_p, rs_p, rn_p &
        !$omp ,re_m, rw_m, rs_m, rn_m, tmp, psi_w, psi_e, psi_s, psi_n, su_DC) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                if (mod(i+j, 2) /= 0) then
                      
                    if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1) then
                        phi_new = (aw_up*phi(i, j-1) + ae_up*phi(i, j+1) + as_up*phi(i+1, j) + an_up*phi(i-1, j)) / ap_up
                    else
                        re_p = (phi(i, j) - phi(i-1, j)) / (phi(i+1, j) - phi(i, j) + C)
                        rw_p = (phi(i-1, j) - phi(i-2, j)) / (phi(i, j) - phi(i-1, j) + C)
                        rs_p = (phi(i, j-1) - phi(i, j-2)) / (phi(i, j) - phi(i, j-1) + C)
                        rn_p = (phi(i, j) - phi(i, j-1)) / (phi(i, j+1) - phi(i, j) + C)
                          
                        re_m = (phi(i+2, j) - phi(i+1, j)) / (phi(i+1, j) - phi(i, j) + C)
                        rw_m = (phi(i+1, j) - phi(i, j)) / (phi(i, j) - phi(i-1, j) + C)
                        rs_m = (phi(i, j+1) - phi(i, j)) / (phi(i, j) - phi(i, j-1) + C)
                        rn_m = (phi(i, j+2) - phi(i, j+1)) / (phi(i, j-1) - phi(i, j) + C)
                         
                        ! Use UMIST model for limiter function
                        if (Fw > 0.0D0) then
                            tmp = min(2.0D0*rw_p, (1.0D0+3.0D0*rw_p)/4.0D0, (3.0D0+rw_p)/4.0D0, 2.0D0)
                            psi_w = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rw_m, (1.0D0+3.0D0*rw_m)/4.0D0, (3.0D0+rw_m)/4.0D0, 2.0D0)
                            psi_w = max(0.0D0, tmp)
                        endif
                          
                        if (Fe > 0.0D0) then
                            tmp = min(2.0D0*re_p, (1.0D0+3.0D0*re_p)/4.0D0, (3.0D0+re_p)/4.0D0, 2.0D0)
                            psi_e = max(0.0D0, tmp) 
                        else
                            tmp = min(2.0D0*re_m, (1.0D0+3.0D0*re_m)/4.0D0, (3.0D0+re_m)/4.0D0, 2.0D0)
                            psi_e = max(0.0D0, tmp) 
                        endif
                          
                        if (Fs > 0.0D0) then
                            tmp = min(2.0D0*rs_p, (1.0D0+3.0D0*rs_p)/4.0D0, (3.0D0+rs_p)/4.0D0, 2.0D0)
                            psi_s = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rs_m, (1.0D0+3.0D0*rs_m)/4.0D0, (3.0D0+rs_m)/4.0D0, 2.0D0)
                            psi_s = max(0.0D0, tmp)
                        endif

                        if (Fn > 0.0D0) then
                            tmp = min(2.0D0*rn_p, (1.0D0+3.0D0*rn_p)/4.0D0, (3.0D0+rn_p)/4.0D0, 2.0D0)
                            psi_n = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rn_m, (1.0D0+3.0D0*rn_m)/4.0D0, (3.0D0+rn_m)/4.0D0, 2.0D0)
                            psi_n = max(0.0D0, tmp)
                        endif
                          
                        tmp = (1.0D0 - alf_e) * psi_e - alf_e * psi_e
                        term1 = 0.5D0 * Fe * (phi(i+1, j) - phi(i, j)) * tmp
                          
                        tmp = alf_w * psi_w - (1.0D0 - alf_w) * psi_w
                        term2 = 0.5D0 * Fw * (phi(i, j) - phi(i-1, j)) * tmp
                      
                        tmp = (1.0D0 - alf_n) * psi_n - alf_n * psi_n
                        term3 = 0.5D0 * Fn * (phi(i, j+1) - phi(i, j)) * tmp
                          
                        tmp = alf_s * psi_s - (1.0D0 - alf_s) * psi_s
                        term4 = 0.5D0 * Fs * (phi(i, j) - phi(i, j-1)) * tmp
                          
                        Su_DC = term1 + term2 + term3 + term4
                          
                        phi_new = (aw_t * phi(i-1, j) + ae_t * + phi(i+1, j) + as_t * phi(i, j-1) &
                                + an_t * phi(i, j+1) + Su_DC) / ap_t
                    endif
                      
                    phi(i,j) = (1.0D0 - rf) * phiold(i, j) + rf * phi_new
                      
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
end program CD_TVD