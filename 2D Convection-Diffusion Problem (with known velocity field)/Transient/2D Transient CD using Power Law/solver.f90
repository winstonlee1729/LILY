program Transient_CD_Power_Law
    use omp_lib
    implicit none
    
    ! The velocity field is known beforehand
    ! Declaration of variables
    integer*4 nx, ny, i, j
    integer*4 time_step, max_step, iter, maxiter, iter_limit
    integer*4 output_freq
    parameter(nx = 300, ny = 300)
    parameter(maxiter = 5000, max_step = 2000)
    parameter(output_freq = 25)
    
    real*8 tol, maxres, err, rhs, lhs, theta
    real*8, allocatable :: phi(:,:), phiold(:,:)
    real*8 rho, u, v, area_w, area_e, area_n, area_s
    real*8 ga_w, ga_e, ga_n, ga_s
    real*8 dx, dy, Lx, Ly, volume, thickness, dt
    real*8 Fe, Fw, Fn, Fs, De, Dw, Dn, Ds, delta_F
    real*8 Pe_e, Pe_w, Pe_n, Pe_s
    real*8 aw, ae, an, as, ap, ap_o, ap_prime
    real*8 w_node, e_node, n_node, s_node, p_past
    parameter(tol = 1.0D-7)
    parameter(rho = 1.0D0, u = 0.0D0, v = 0.001D0)
    parameter(ga_w = 0.01D0, ga_e = 0.01D0, ga_n = 0.01D0, ga_s = 0.01D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, dt = 0.05D0, thickness = 1.0D0)
    
    ! 1.0 = Fully Implicit | 0.5 = Crank-Nicolson | 0.0 = Explcit
    parameter(theta = 0.5D0)
    
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
    
    Pe_e = Fe / De
    Pe_w = Fw / Dw
    Pe_n = Fn / Dn
    Pe_s = Fs / Ds
    delta_F = Fe - Fw + Fn - Fs
    
    ap_o = rho * volume / dt
    aw = Dw * max(0.0D0, (1.0D0 - 0.1D0*abs(Pe_w))**5) + max(Fw, 0.0D0)
    ae = De * max(0.0D0, (1.0D0 - 0.1D0*abs(Pe_e))**5) + max(-Fe, 0.0D0)
    as = Ds * max(0.0D0, (1.0D0 - 0.1D0*abs(Pe_s))**5) + max(Fs, 0.0D0)
    an = Dn * max(0.0D0, (1.0D0 - 0.1D0*abs(Pe_n))**5) + max(-Fn, 0.0D0)
    ap = ap_o + (ae + aw + as + an +  delta_F) * theta
    ap_prime = ap_o - (ae + aw + an + as + delta_F) * (1 - theta)
    
    ! Allocate memory
    allocate(phi(nx, ny), phiold(nx, ny))
    
    ! Initialize flow field
    phi = 0.0D0
    phiold = 0.0D0
    
    ! Applying boundary condition
    !$omp parallel do
    do i = 1, nx
        phi(i, 1) = 100.0D0
        phiold(i, 1) = 100.0D0
    enddo
    !$omp end parallel do
    
    ! Open the file
    open(unit=10, file='phi_out.txt', status='replace')
    
    ! Check for inner iteration number
    if (theta == 0.0D0) then
        iter_limit = 1
    else
        iter_limit = maxiter
    endif
    
    
    ! Start time step 
    do time_step = 1, max_step
        
        ! Save old filed
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
            
            !$omp parallel do private(j, e_node, w_node, n_node, s_node, p_past, rhs, lhs, err) reduction(max:maxres)
            do i = 2, nx-1
                do j = 2, ny-1
                    w_node = aw * (theta * phi(i-1, j) + (1 - theta) * phiold(i-1, j))
                    e_node = ae * (theta * phi(i+1 ,j) + (1 - theta) * phiold(i+1, j))
                    s_node = as * (theta * phi(i, j-1) + (1 - theta) * phiold(i, j-1))
                    n_node = an * (theta * phi(i, j+1) + (1 - theta) * phiold(i, j+1))
                    p_past = ap_prime * phiold(i, j)
                    
                    rhs = w_node + e_node + s_node + n_node + p_past
                    lhs = ap * phi(i, j)
                    
                    err = abs(rhs - lhs)
                    if (err >= maxres) maxres = err
                    
                    phi(i, j) = rhs / ap
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
    deallocate(phi, phiold)
    
    stop
    end program Transient_CD_Power_Law