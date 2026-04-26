program CD_TVD
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
    real*8 phi_new
    
    ! Parameters for convergence
    real*8 maxres, err
    real*8, parameter :: tol = 1.0D-5, rf = 0.6D0
    
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
    real*8 alf_w, alf_e, alf_n, alf_s, alf_t, alf_b
    real*8 aw_hb, ae_hb, as_hb, an_hb, at_hb, ab_hb, ap_hb
    real*8 aw_t, ae_t, as_t, an_t, at_t, ab_t, ap_t, Su_DC
    real*8 rw_p, re_p, rs_p, rn_p, rt_p, rb_p
    real*8 rw_m, re_m, rs_m, rn_m, rt_m, rb_m
    real*8 psi_w, psi_e, psi_n, psi_s, psi_t, psi_b, tmp
    real*8 term1, term2, term3, term4, term5, term6
    real*8 x_term, y_term, z_term, C
    real*8 xc, yc, zc
    parameter(C = 1.0D-15)
    
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
    
    ! Flow direction
    if (Fw > 0) then; alf_w = 1.0D0; else; alf_w = 0.0D0; endif
    if (Fe > 0) then; alf_e = 1.0D0; else; alf_e = 0.0D0; endif
    if (Fs > 0) then; alf_s = 1.0D0; else; alf_s = 0.0D0; endif
    if (Fn > 0) then; alf_n = 1.0D0; else; alf_n = 0.0D0; endif
    if (Ft > 0) then; alf_t = 1.0D0; else; alf_t = 0.0D0; endif
    if (Fb > 0) then; alf_b = 1.0D0; else; alf_b = 0.0D0; endif
    
    ! Hybrid scheme coefficients
    aw_hb = max(Fw, Dw + (Fw / 2.0D0), 0.0D0)
    ae_hb = max(-Fe, De - (Fe / 2.0D0), 0.0D0)
    as_hb = max(Fs, Ds + (Fs / 2.0D0), 0.0D0)
    an_hb = max(-Fn, Dn - (Fn / 2.0D0), 0.0D0)
    ab_hb = max(Fb, Db + (Fb / 2.0D0), 0.0D0)
    at_hb = max(-Ft, Dt - (Ft / 2.0D0), 0.0D0)
    
    ap_hb = aw_hb + ae_hb + as_hb + an_hb + ab_hb + at_hb + delta_F
    
    ! TVD scheme coefficients 
    aw_t = Dw + max(Fw, 0.0D0)
    ae_t = De + max(-Fe, 0.0D0)
    as_t = Ds + max(Fs, 0.0D0)
    an_t = Dn + max(-Fn, 0.0D0)
    ab_t = Db + max(Fb, 0.0D0)
    at_t = Dt + max(-Ft, 0.0D0)
    
    ap_t = aw_t + ae_t + as_t + an_t + at_t + ab_t + delta_F
    
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
        
        ! Pass 1: Red cells
        !$omp parallel do private(i, j, k, err, phi_new, re_p, rw_p, rn_p, rs_p, rt_p, rb_p &
        !$omp ,re_m, rw_m, rn_m, rs_m, rt_m, rb_m, tmp, psi_w, psi_e, psi_s, psi_n, psi_t, psi_b &
        !$omp ,x_term, y_term, z_term, term1, term2, term3, term4, term5, term6, Su_DC) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if(mod(i+j+k, 2) == 0) then
                        if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1 .or. k==2 .or. k==nz-1) then
                            x_term = aw_hb * phi(i-1, j, k) + ae_hb * phi(i+1, j, k)
                            y_term = as_hb * phi(i, j-1, k) + an_hb * phi(i, j+1, k)
                            z_term = ab_hb * phi(i, j, k-1) + at_hb * phi(i, j, k+1)
                            phi_new = (x_term + y_term + z_term) / ap_hb
                        else

                            re_p = (phi(i, j, k) - phi(i-1, j, k)) / (phi(i+1, j, k) - phi(i, j, k) + C)
                            rw_p = (phi(i-1, j, k) - phi(i-2, j, k)) / (phi(i, j, k) - phi(i-1, j, k) + C)
                            rn_p = (phi(i, j, k) - phi(i, j-1, k)) / (phi(i, j+1, k) - phi(i, j, k) + C)
                            rs_p = (phi(i, j-1, k) - phi(i, j-2, k)) / (phi(i, j, k) - phi(i, j-1, k) + C)
                            rt_p = (phi(i, j, k) - phi(i, j, k-1)) / (phi(i, j, k+1) - phi(i, j, k) + C)
                            rb_p = (phi(i, j, k-1) - phi(i, j, k-2)) / (phi(i, j, k) - phi(i, j, k-1) + C)
                            
                            re_m = (phi(i+2, j, k) - phi(i+1, j, k)) / (phi(i+1, j, k) - phi(i, j, k) + C)
                            rw_m = (phi(i+1, j, k) - phi(i, j, k)) / (phi(i, j, k) - phi(i-1, j, k) + C)
                            rn_m = (phi(i, j+2, k) - phi(i, j+1, k)) / (phi(i, j+1, k) - phi(i, j, k) + C) 
                            rs_m = (phi(i, j+1, k) - phi(i, j, k)) / (phi(i, j, k) - phi(i, j-1, k) + C)
                            rt_m = (phi(i, j, k+2) - phi(i, j, k+1)) / (phi(i, j, k+1) - phi(i, j, k) + C)
                            rb_m = (phi(i, j, k+1) - phi(i, j, k)) / (phi(i, j, k) - phi(i, j, k-1) + C)
                            
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
                            
                            if (Fb > 0.0D0) then
                                tmp = min(2.0D0*rb_p, (1.0D0+3.0D0*rb_p)/4.0D0, (3.0D0+rb_p)/4.0D0, 2.0D0)
                                psi_b = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rb_m, (1.0D0+3.0D0*rb_m)/4.0D0, (3.0D0+rb_m)/4.0D0, 2.0D0)
                                psi_b = max(0.0D0, tmp)
                            endif

                            if (Ft > 0.0D0) then
                                tmp = min(2.0D0*rt_p, (1.0D0+3.0D0*rt_p)/4.0D0, (3.0D0+rt_p)/4.0D0, 2.0D0)
                                psi_t = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rt_m, (1.0D0+3.0D0*rt_m)/4.0D0, (3.0D0+rt_m)/4.0D0, 2.0D0)
                                psi_t = max(0.0D0, tmp)
                            endif
                            
                            tmp = (1.0D0 - alf_e) * psi_e - alf_e * psi_e
                            term1 = 0.5D0 * Fe * (phi(i+1, j, k) - phi(i, j, k)) * tmp
                          
                            tmp = alf_w * psi_w - (1.0D0 - alf_w) * psi_w
                            term2 = 0.5D0 * Fw * (phi(i, j, k) - phi(i-1, j, k)) * tmp
                      
                            tmp = (1.0D0 - alf_n) * psi_n - alf_n * psi_n
                            term3 = 0.5D0 * Fn * (phi(i, j+1, k) - phi(i, j, k)) * tmp
                          
                            tmp = alf_s * psi_s - (1.0D0 - alf_s) * psi_s
                            term4 = 0.5D0 * Fs * (phi(i, j, k) - phi(i, j-1, k)) * tmp
                            
                            tmp = (1.0D0 - alf_t) * psi_t - alf_t * psi_t
                            term5 = 0.5D0 * Ft * (phi(i, j, k+1) - phi(i, j, k)) * tmp
                          
                            tmp = alf_b * psi_b - (1.0D0 - alf_b) * psi_b
                            term6 = 0.5D0 * Fb * (phi(i, j, k-1) - phi(i, j, k)) * tmp
                            
                          
                            Su_DC = term1 + term2 + term3 + term4 + term5 + term6
                            
                            x_term = aw_t * phi(i-1, j, k) + ae_t * phi(i+1, j, k)
                            y_term = an_t * phi(i, j+1, k) + as_t * phi(i, j-1, k)
                            z_term = ab_t * phi(i, j, k-1) + at_t * phi(i, j, k+1)
                            
                            phi_new = (x_term + y_term + z_term + Su_DC) / ap_t
                        endif
                        
                        phi(i, j, k) = (1.0D0 - rf) * phiold(i, j, k) + rf * phi_new
                        
                        err = abs(phi(i, j, k) - phiold(i, j, k))
                        if (err > maxres) maxres = err
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do

        ! Pass 2: Black cells
        !$omp parallel do private(i, j, k, err, phi_new, re_p, rw_p, rn_p, rs_p, rt_p, rb_p &
        !$omp ,re_m, rw_m, rn_m, rs_m, rt_m, rb_m, tmp, psi_w, psi_e, psi_s, psi_n, psi_t, psi_b &
        !$omp ,x_term, y_term, z_term, term1, term2, term3, term4, term5, term6, Su_DC) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if(mod(i+j+k, 2) /= 0) then
                        if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1 .or. k==2 .or. k==nz-1) then
                            x_term = aw_hb * phi(i-1, j, k) + ae_hb * phi(i+1, j, k)
                            y_term = as_hb * phi(i, j-1, k) + an_hb * phi(i, j+1, k)
                            z_term = ab_hb * phi(i, j, k-1) + at_hb * phi(i, j, k+1)
                            phi_new = (x_term + y_term + z_term) / ap_hb
                        else

                            re_p = (phi(i, j, k) - phi(i-1, j, k)) / (phi(i+1, j, k) - phi(i, j, k) + C)
                            rw_p = (phi(i-1, j, k) - phi(i-2, j, k)) / (phi(i, j, k) - phi(i-1, j, k) + C)
                            rn_p = (phi(i, j, k) - phi(i, j-1, k)) / (phi(i, j+1, k) - phi(i, j, k) + C)
                            rs_p = (phi(i, j-1, k) - phi(i, j-2, k)) / (phi(i, j, k) - phi(i, j-1, k) + C)
                            rt_p = (phi(i, j, k) - phi(i, j, k-1)) / (phi(i, j, k+1) - phi(i, j, k) + C)
                            rb_p = (phi(i, j, k-1) - phi(i, j, k-2)) / (phi(i, j, k) - phi(i, j, k-1) + C)
                            
                            re_m = (phi(i+2, j, k) - phi(i+1, j, k)) / (phi(i+1, j, k) - phi(i, j, k) + C)
                            rw_m = (phi(i+1, j, k) - phi(i, j, k)) / (phi(i, j, k) - phi(i-1, j, k) + C)
                            rn_m = (phi(i, j+2, k) - phi(i, j+1, k)) / (phi(i, j+1, k) - phi(i, j, k) + C) 
                            rs_m = (phi(i, j+1, k) - phi(i, j, k)) / (phi(i, j, k) - phi(i, j-1, k) + C)
                            rt_m = (phi(i, j, k+2) - phi(i, j, k+1)) / (phi(i, j, k+1) - phi(i, j, k) + C)
                            rb_m = (phi(i, j, k+1) - phi(i, j, k)) / (phi(i, j, k) - phi(i, j, k-1) + C)
                            
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
                            
                            if (Fb > 0.0D0) then
                                tmp = min(2.0D0*rb_p, (1.0D0+3.0D0*rb_p)/4.0D0, (3.0D0+rb_p)/4.0D0, 2.0D0)
                                psi_b = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rb_m, (1.0D0+3.0D0*rb_m)/4.0D0, (3.0D0+rb_m)/4.0D0, 2.0D0)
                                psi_b = max(0.0D0, tmp)
                            endif

                            if (Ft > 0.0D0) then
                                tmp = min(2.0D0*rt_p, (1.0D0+3.0D0*rt_p)/4.0D0, (3.0D0+rt_p)/4.0D0, 2.0D0)
                                psi_t = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rt_m, (1.0D0+3.0D0*rt_m)/4.0D0, (3.0D0+rt_m)/4.0D0, 2.0D0)
                                psi_t = max(0.0D0, tmp)
                            endif
                            
                            tmp = (1.0D0 - alf_e) * psi_e - alf_e * psi_e
                            term1 = 0.5D0 * Fe * (phi(i+1, j, k) - phi(i, j, k)) * tmp
                          
                            tmp = alf_w * psi_w - (1.0D0 - alf_w) * psi_w
                            term2 = 0.5D0 * Fw * (phi(i, j, k) - phi(i-1, j, k)) * tmp
                      
                            tmp = (1.0D0 - alf_n) * psi_n - alf_n * psi_n
                            term3 = 0.5D0 * Fn * (phi(i, j+1, k) - phi(i, j, k)) * tmp
                          
                            tmp = alf_s * psi_s - (1.0D0 - alf_s) * psi_s
                            term4 = 0.5D0 * Fs * (phi(i, j, k) - phi(i, j-1, k)) * tmp
                            
                            tmp = (1.0D0 - alf_t) * psi_t - alf_t * psi_t
                            term5 = 0.5D0 * Ft * (phi(i, j, k+1) - phi(i, j, k)) * tmp
                          
                            tmp = alf_b * psi_b - (1.0D0 - alf_b) * psi_b
                            term6 = 0.5D0 * Fb * (phi(i, j, k-1) - phi(i, j, k)) * tmp
                            
                          
                            Su_DC = term1 + term2 + term3 + term4 + term5 + term6
                            
                            x_term = aw_t * phi(i-1, j, k) + ae_t * phi(i+1, j, k)
                            y_term = an_t * phi(i, j+1, k) + as_t * phi(i, j-1, k)
                            z_term = ab_t * phi(i, j, k-1) + at_t * phi(i, j, k+1)
                            
                            phi_new = (x_term + y_term + z_term + Su_DC) / ap_t
                        endif
                        
                        phi(i, j, k) = (1.0D0 - rf) * phiold(i, j, k) + rf * phi_new
                        
                        err = abs(phi(i, j, k) - phiold(i, j, k))
                        if (err > maxres) maxres = err
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do

        write(*,*) 'Iteration: ', iter, '|Max Residual: ', maxres

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
end program CD_TVD