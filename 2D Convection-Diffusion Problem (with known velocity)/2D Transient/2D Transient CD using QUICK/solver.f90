program Transient_CD_QUICK
    use omp_lib
    implicit none
    
! The velocity field is known beforehand
! Declaration of variables
    ! Grid size, time step, and iteration number
    integer*4 nx, ny
    integer*4 i, j
    integer*4 time_step, max_step
    integer*4 iter, maxiter, iter_limit
    integer*4 output_freq
    parameter(nx = 300, ny = 300)
    parameter(maxiter = 5000, max_step = 1000)
    parameter(output_freq = 25)
    
    ! Program convergence parameters
    real*8 tol, maxres, err, rhs, lhs

    ! Flow field
    real*8, allocatable :: phi(:,:), phiold(:,:), temp(:,:)
    
    ! Flow parameters
    real*8 rho, u, v 
    real*8 area_w, area_e, area_n, area_s
    real*8 ga_w, ga_e, ga_n, ga_s
    real*8 dx, dy, Lx, Ly, volume, thickness, dt
    real*8 Fe, Fw, Fn, Fs, delta_F
    real*8 De, Dw, Dn, Ds
    parameter(tol = 1.0D-7)
    parameter(rho = 1.0D0, u = 0.0D0, v = 0.001D0)
    parameter(ga_w = 0.01D0, ga_e = 0.01D0, ga_n = 0.01D0, ga_s = 0.01D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, dt = 0.05D0, thickness = 1.0D0)
    
    ! Temporal discretization scheme
    real*8 theta
    ! 1.0 = Fully Implicit | 0.5 = Crank-Nicolson | 0.0 = Explcit
    parameter(theta = 0.5D0)
    
    ! Alpha variables for flow direction
    real*8 alf_w, alf_e, alf_s, alf_n
    
    ! Upwind coefficients
    real*8 aw_up, ae_up, as_up, an_up
    real*8 ap_up_s, ap_up, ap_o, ap_prime_up
    
    ! QUICK coefficients
    real*8 aw_qk, ae_qk, an_qk, as_qk
    real*8 aww_qk, aee_qk, ann_qk, ass_qk
    real*8 ap_qk_s, ap_qk, ap_prime_qk
    
    ! Nodal values for solver
    real*8 w_node, e_node, n_node, s_node, p_past
    real*8 ww_node, ee_node, nn_node, ss_node
    
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    
    area_w = dy * thickness
    area_e = area_w
    area_n = dx * thickness
    area_s = area_n
    
    volume = dx * dy * thickness
    
    Fw = rho * u * area_w
    Fe = rho * u * area_e
    Fn = rho * v * area_n
    Fs = rho * v * area_s
    
    Dw = (ga_w * area_w) / dx
    De = (ga_e * area_e) / dx
    Dn = (ga_n * area_n) / dy
    Ds = (ga_s * area_s) / dy
    
    delta_F = Fe - Fw + Fn - Fs
    ap_o = rho * volume / dt
    
    ! Flow direction 
    if (Fw > 0.0D0) then; alf_w = 1.0D0; else; alf_w = 0.0D0; endif
    if (Fe > 0.0D0) then; alf_e = 1.0D0; else; alf_e = 0.0D0; endif
    if (Fs > 0.0D0) then; alf_s = 1.0D0; else; alf_s = 0.0D0; endif
    if (Fn > 0.0D0) then; alf_n = 1.0D0; else; alf_n = 0.0D0; endif
    
    ! Upwind coefficients
    aw_up = max(Fw, Dw + (Fw / 2.0D0), 0.0D0)
    ae_up = max(-Fe, De - (Fe / 2.0D0), 0.0D0)
    as_up = max(Fs, Ds + (Fs / 2.0D0), 0.0D0)
    an_up = max(-Fn, Dn - (Fn / 2.0D0), 0.0D0)
    
    ap_up_s = aw_up + ae_up + as_up + an_up + delta_F
    ap_up = ap_o + ap_up_s * theta
    ap_prime_up = ap_o - ap_up_s * (1.0D0 - theta)
    
    ! QUICK coefficients 
    aw_qk = Dw + (6.0D0/8.0D0)*alf_w*Fw + (1.0D0/8.0D0)*alf_e*Fe + (3.0D0/8.0D0)*(1.0D0-alf_w)*Fw
    aww_qk = (-1.0D0/8.0D0)*alf_w*Fw

    ae_qk = De - (3.0D0/8.0D0)*alf_e*Fe - (6.0D0/8.0D0)*(1.0D0-alf_e)*Fe - (1.0D0/8.0D0)*(1.0D0-alf_w)*Fw
    aee_qk = (1.0D0/8.0D0)*(1.0D0-alf_e)*Fe

    as_qk = Ds + (6.0D0/8.0D0)*alf_s*Fs + (1.0D0/8.0D0)*alf_n*Fn + (3.0D0/8.0D0)*(1.0D0-alf_s)*Fs
    ass_qk = (-1.0D0/8.0D0)*alf_s*Fs
 
    an_qk = Dn - (3.0D0/8.0D0)*alf_n*Fn - (6.0D0/8.0D0)*(1.0D0-alf_n)*Fn - (1.0D0/8.0D0)*(1.0D0-alf_s)*Fs
    ann_qk = (1.0D0/8.0D0)*(1.0D0-alf_n)*Fn
    
    ap_qk_s = aw_qk + ae_qk + as_qk + an_qk + aww_qk + aee_qk + ass_qk + ann_qk + delta_F
    ap_qk = ap_o + ap_qk_s * theta
    ap_prime_qk = ap_o - ap_qk_s * (1.0D0 - theta)
    
! Allocate memory
    allocate(phi(nx, ny), phiold(nx, ny), temp(nx, ny))
    
! Initialization
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
    
! Open the output file
    open(unit=10, file='phi_out.txt', status='replace')
    
! Check for inner iteration number 
    if (theta == 0.0D0) then
        iter_limit = 1
    else
        iter_limit = maxiter
    endif
    
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
            
            !$omp parallel do private(j, e_node, w_node, n_node, s_node, p_past, &
            !$omp ww_node, ee_node, nn_node, ss_node, rhs, lhs, err) reduction(max:maxres)
            do i = 2, nx-1
                do j = 2, ny-1
                    
                    ! Safety ring
                    if (i == 2 .or. i == nx-1 .or. j == 2 .or. j == ny-1) then
                        w_node = aw_up * (theta * phi(i-1, j) + (1.0D0 - theta) * phiold(i-1, j))
                        e_node = ae_up * (theta * phi(i+1, j) + (1.0D0 - theta) * phiold(i+1, j))
                        s_node = as_up * (theta * phi(i, j-1) + (1.0D0 - theta) * phiold(i, j-1))
                        n_node = an_up * (theta * phi(i, j+1) + (1.0D0 - theta) * phiold(i, j+1))
                        p_past = ap_prime_up * phiold(i, j)
                        
                        rhs = w_node + e_node + s_node + n_node + p_past
                        lhs = ap_up * phi(i, j)
                        temp(i, j) = rhs / ap_up
                    
                    else
                    ! QUICK 
                        w_node = aw_qk * (theta * phi(i-1, j) + (1.0D0 - theta) * phiold(i-1, j))
                        e_node = ae_qk * (theta * phi(i+1, j) + (1.0D0 - theta) * phiold(i+1, j))
                        s_node = as_qk * (theta * phi(i, j-1) + (1.0D0 - theta) * phiold(i, j-1))
                        n_node = an_qk * (theta * phi(i, j+1) + (1.0D0 - theta) * phiold(i, j+1))
                        
                        ww_node = aww_qk * (theta * phi(i-2, j) + (1.0D0 - theta) * phiold(i-2, j))
                        ee_node = aee_qk * (theta * phi(i+2, j) + (1.0D0 - theta) * phiold(i+2, j))
                        ss_node = ass_qk * (theta * phi(i, j-2) + (1.0D0 - theta) * phiold(i, j-2))
                        nn_node = ann_qk * (theta * phi(i, j+2) + (1.0D0 - theta) * phiold(i, j+2))
                        
                        p_past = ap_prime_qk * phiold(i, j)
                        
                        rhs = w_node + e_node + s_node + n_node + ww_node + ee_node + ss_node + nn_node + p_past
                        lhs = ap_qk * phi(i, j)
                        temp(i, j) = rhs / ap_qk
                    endif
                    
                    ! Evaluate local residual
                    err = abs(rhs - lhs)
                    if (err >= maxres) maxres = err
                    
                enddo
            enddo
            !$omp end parallel do
            
            ! Update the main phi field with the newly calculated values
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
        
        ! Outputting solution
        if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
            do j = 1, ny
                write(10, '(1000E15.6)') (phi(i, j), i = 1, nx)
            enddo
        endif
        
    enddo
    
! Close the file
    close(10)
    
! Deallocate memory
    deallocate(phi, phiold, temp)
    
    stop
end program Transient_CD_QUICK