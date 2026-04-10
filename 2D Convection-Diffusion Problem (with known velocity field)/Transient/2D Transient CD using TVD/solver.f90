program Transient_CD_TVD
    use omp_lib
    implicit none
    
    ! Declarations of variables
    integer*4 nx, ny, iter, maxiter, iter_limit, i, j
    integer*4 time_step, max_step, output_freq
    parameter(nx = 300, ny = 300, maxiter = 5000, max_step = 2000)
    parameter(output_freq = 25)
    
    real*8, allocatable :: phi(:, :), phiold(:, :), temp(:, :)
    real*8 tol, maxres, err, rhs, lhs
    real*8 u, v, rho, dx, dy, volume, dt, theta, Lx, Ly
    real*8 Dw, De, Ds, Dn, Fw, Fe, Fs, Fn, delta_F
    real*8 ga_w, ga_e, ga_n, ga_s
    real*8 area_w, area_e, area_s, area_n, thickness
    real*8 alf_w, alf_e, alf_s, alf_n
    
    ! Upwind variables 
    real*8 aw_up, ae_up, as_up, an_up
    real*8 ap_up_s, ap_up, ap_prime_up
    
    ! TVD variables 
    real*8 aw_t, ae_t, as_t, an_t
    real*8 ap_t_s, ap_t, ap_prime_t, ap_o
    real*8 Su_DC
    
    real*8 w_node, e_node, n_node, s_node, p_past
    
    ! TVD ratios and limiters
    real*8 rw_p, re_p, rs_p, rn_p
    real*8 rw_m, re_m, rs_m, rn_m
    real*8 psi_w, psi_e, psi_s, psi_n, tmp
    real*8 term1, term2, term3, term4, C
    
    parameter(tol = 1.0D-7)
    parameter(u = 0.0D0, v = 0.01D0, rho = 1.0D0)
    parameter(ga_w = 0.01D0, ga_e = 0.01D0, ga_n = 0.01D0, ga_s = 0.01D0)
    parameter(dt = 0.05D0, C = 1.0D-15)
    parameter(Lx = 1.0D0, Ly = 1.0D0, thickness = 1.0D0)
    
    ! 1.0 = Fully Implicit | 0.5 = Crank-Nicolson | 0.0 = Explcit
    parameter(theta = 0.5D0)
    
    ! Allocate memory
    allocate(phi(nx,ny), phiold(nx,ny), temp(nx,ny))
    
    ! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    volume = dx * dy * thickness
    area_w = dy * thickness
    area_e = area_w
    area_n = dx * thickness
    area_s = area_n
    
    Fw = rho * u * area_w
    Fe = rho * u * area_e
    Fs = rho * v * area_s
    Fn = rho * v * area_n
    
    Dw = (ga_w / dx) * area_w
    De = (ga_e / dx) * area_e
    Ds = (ga_s / dy) * area_s
    Dn = (ga_n / dy) * area_n
    
    delta_F = (Fe - Fw) + (Fn - Fs)
    ap_o = rho * volume / dt

    ! Flow direction 
    if (Fw > 0.0D0) then; alf_w = 1.0D0; else; alf_w = 0.0D0; endif
    if (Fe > 0.0D0) then; alf_e = 1.0D0; else; alf_e = 0.0D0; endif
    if (Fs > 0.0D0) then; alf_s = 1.0D0; else; alf_s = 0.0D0; endif
    if (Fn > 0.0D0) then; alf_n = 1.0D0; else; alf_n = 0.0D0; endif
    
    ! Upwind coefficients (Safety ring)
    aw_up = max(Fw, Dw + (Fw/2.0D0), 0.0D0)
    ae_up = max(-Fe, De - (Fe/2.0D0), 0.0D0)
    as_up = max(Fs, Ds + (Fs/2.0D0), 0.0D0)
    an_up = max(-Fn, Dn - (Fn/2.0D0), 0.0D0)
    
    ap_up_s = aw_up + ae_up + as_up + an_up + delta_F
    ap_up = ap_o + ap_up_s * theta
    ap_prime_up = ap_o - ap_up_s * (1.0D0 - theta)

    ! TVD coefficients
    aw_t = Dw + max(Fw, 0.0D0)
    ae_t = De + max(-Fe, 0.0D0)
    as_t = Ds + max(Fs, 0.0D0)
    an_t = Dn + max(-Fn, 0.0D0)
    
    ap_t_s = aw_t + ae_t + as_t + an_t + delta_F
    ap_t = ap_o + ap_t_s * theta
    ap_prime_t = ap_o - ap_t_s * (1.0D0 - theta)
    
    ! Initialization of field
    phi = 0.0D0
    phiold = 0.0D0
    temp = 0.0D0
    
    ! Applying boundary condition
    !$omp parallel do
    do i = 1, nx
        phi(i, 1) = 100.0D0
        phiold(i, 1) = 100.0D0
        temp(i, 1) = 100.0D0
    enddo
    !$omp end parallel do
    
    ! Check for inner iteration number 
    if (theta == 0.0D0) then
        iter_limit = 1
    else
        iter_limit = maxiter
    endif
    
    ! Open the output file before the time loop starts
    open(unit=10, file='phi_out.txt', status='replace')

    ! Start time step
    do time_step = 1, max_step
        
        ! Save old field
        !$omp parallel do collapse(2)
        do i = 1, nx
            do j = 1, ny
                phiold(i, j) = phi(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        ! Inner iteration
        do iter = 1, iter_limit
            maxres = 0.0D0
            
            !$omp parallel do private(j, w_node, e_node, s_node, n_node, p_past, rhs, lhs, err, &
            !$omp re_p, rw_p, rs_p, rn_p, re_m, rw_m, rs_m, rn_m, tmp, psi_w, psi_e, psi_s, psi_n, &
            !$omp term1, term2, term3, term4, Su_DC) reduction(max:maxres)
            do i = 2, nx-1
                do j = 2, ny-1
                    
                    ! Boundary safety ring
                    if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1) then
                        w_node = aw_up * (theta * phi(i-1, j) + (1.0D0 - theta) * phiold(i-1, j))
                        e_node = ae_up * (theta * phi(i+1, j) + (1.0D0 - theta) * phiold(i+1, j))
                        s_node = as_up * (theta * phi(i, j-1) + (1.0D0 - theta) * phiold(i, j-1))
                        n_node = an_up * (theta * phi(i, j+1) + (1.0D0 - theta) * phiold(i, j+1))
                        p_past = ap_prime_up * phiold(i, j)
                        
                        rhs = w_node + e_node + s_node + n_node + p_past
                        lhs = ap_up * phi(i, j)
                        temp(i, j) = rhs / ap_up
                        
                    ! TVD 
                    else
                        ! Positive flow ratios 
                        re_p = (phi(i,j) - phi(i-1,j)) / (phi(i+1,j) - phi(i,j) + C)
                        rw_p = (phi(i-1,j) - phi(i-2,j)) / (phi(i,j) - phi(i-1,j) + C)
                        rn_p = (phi(i,j) - phi(i,j-1)) / (phi(i,j+1) - phi(i,j) + C)
                        rs_p = (phi(i,j-1) - phi(i,j-2)) / (phi(i,j) - phi(i,j-1) + C)
                        
                        ! Negative flow ratios
                        re_m = (phi(i+1,j) - phi(i+2,j)) / (phi(i,j) - phi(i+1,j) + C)
                        rw_m = (phi(i,j) - phi(i+1,j)) / (phi(i-1,j) - phi(i,j) + C)
                        rn_m = (phi(i,j+1) - phi(i,j+2)) / (phi(i,j) - phi(i,j+1) + C)
                        rs_m = (phi(i,j) - phi(i,j+1)) / (phi(i,j-1) - phi(i,j) + C)
                        
                        ! Use UMIST for limter function
                        ! West limiter
                        if (Fw > 0.0D0) then
                            tmp = min(2.0D0*rw_p, (1.0D0+3.0D0*rw_p)/4.0D0, (3.0D0+rw_p)/4.0D0, 2.0D0)
                            psi_w = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rw_m, (1.0D0+3.0D0*rw_m)/4.0D0, (3.0D0+rw_m)/4.0D0, 2.0D0)
                            psi_w = max(0.0D0, tmp)
                        endif
                        
                        ! East limiter
                        if (Fe > 0.0D0) then
                            tmp = min(2.0D0*re_p, (1.0D0+3.0D0*re_p)/4.0D0, (3.0D0+re_p)/4.0D0, 2.0D0)
                            psi_e = max(0.0D0, tmp) 
                        else
                            tmp = min(2.0D0*re_m, (1.0D0+3.0D0*re_m)/4.0D0, (3.0D0+re_m)/4.0D0, 2.0D0)
                            psi_e = max(0.0D0, tmp) 
                        endif
                        
                        ! South limiter
                        if (Fs > 0.0D0) then
                            tmp = min(2.0D0*rs_p, (1.0D0+3.0D0*rs_p)/4.0D0, (3.0D0+rs_p)/4.0D0, 2.0D0)
                            psi_s = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rs_m, (1.0D0+3.0D0*rs_m)/4.0D0, (3.0D0+rs_m)/4.0D0, 2.0D0)
                            psi_s = max(0.0D0, tmp)
                        endif
                        
                        ! North limiter
                        if (Fn > 0.0D0) then
                            tmp = min(2.0D0*rn_p, (1.0D0+3.0D0*rn_p)/4.0D0, (3.0D0+rn_p)/4.0D0, 2.0D0)
                            psi_n = max(0.0D0, tmp)
                        else
                            tmp = min(2.0D0*rn_m, (1.0D0+3.0D0*rn_m)/4.0D0, (3.0D0+rn_m)/4.0D0, 2.0D0)
                            psi_n = max(0.0D0, tmp)
                        endif
                        
                        ! Deferred correction source terms
                        tmp = (1.0D0 - alf_e) * psi_e - alf_e * psi_e
                        term1 = 0.5D0 * Fe * (phi(i+1,j) - phi(i,j)) * tmp
                        
                        tmp = alf_w * psi_w - (1.0D0 - alf_w) * psi_w
                        term2 = 0.5D0 * Fw * (phi(i,j) - phi(i-1,j)) * tmp
                        
                        tmp = (1.0D0 - alf_n) * psi_n - alf_n * psi_n
                        term3 = 0.5D0 * Fn * (phi(i,j+1) - phi(i,j)) * tmp
                        
                        tmp = alf_s * psi_s - (1.0D0 - alf_s) * psi_s
                        term4 = 0.5D0 * Fs * (phi(i,j) - phi(i,j-1)) * tmp
                        
                        Su_DC = term1 + term2 + term3 + term4
                        
                        
                        w_node = aw_t * (theta * phi(i-1, j) + (1.0D0 - theta) * phiold(i-1, j))
                        e_node = ae_t * (theta * phi(i+1, j) + (1.0D0 - theta) * phiold(i+1, j))
                        s_node = as_t * (theta * phi(i, j-1) + (1.0D0 - theta) * phiold(i, j-1))
                        n_node = an_t * (theta * phi(i, j+1) + (1.0D0 - theta) * phiold(i, j+1))
                        p_past = ap_prime_t * phiold(i, j)
                        
                        rhs = w_node + e_node + s_node + n_node + p_past + Su_DC
                        lhs = ap_t * phi(i, j)
                        temp(i, j) = rhs / ap_t
                    endif
                    
                    ! Evaluate local residual
                    err = abs(rhs - lhs)
                    if (err >= maxres) maxres = err
                    
                enddo
            enddo
            !$omp end parallel do
            
            ! Update main phi field for the next iteration
            !$omp parallel do collapse(2)
            do i = 2, nx-1
                do j = 2, ny-1
                    phi(i, j) = temp(i, j)
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
        
        ! Outputting solution at intervals (or every time step if output_freq = 1)
        if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
            do j = 1, ny
                write(10, '(1000E15.6)') (phi(i, j), i = 1, nx)
            enddo
        endif
        
    enddo
    
    close(10)
    
    deallocate(phi, phiold, temp)
    
    stop
end program Transient_CD_TVD