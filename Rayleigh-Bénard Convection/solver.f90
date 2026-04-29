program Rayleigh_Benard_Convection
    use omp_lib
    implicit none

! Declaration of variables
    ! Grid size and iteration numbers
    integer*4 nx, ny
    integer*4 i, j, iter, p_iter, time_step
    integer*4 maxiter, max_step, output_freq
    parameter(nx = 200, ny = 100)
    parameter(maxiter = 500, max_step = 5000, output_freq = 10)

    ! For program convergence
    real*8, parameter :: tol = 1.0D-5, sor_tol = 1.0D-7
    
    ! Flow parameter
    real*8 rho, Ra, Pr
    real*8 T_hot, T_cold, T_ref
    real*8 mu, thermal_alpha
    real*8 Lx, Ly, dx, dy, volume, delta_t
    parameter(rho = 1.0D0, Ra = 1.0D5, Pr = 0.71D0)
    parameter(T_hot = 100.0D0, T_cold = 0.0D0, T_ref = 0.5D0)
    parameter(Lx = 2.0D0, Ly = 1.0D0)
    parameter(delta_t = 0.005D0)

    ! Under-relaxation factors 
    real*8, parameter :: alph_u = 0.7D0, alph_v = 0.7D0, alph_p = 0.3D0, alph_T = 0.8D0
    
    ! SOR relaxation factor
    real*8, parameter :: omega = 1.7D0
    
    ! Theta for temporal discretization scheme
    real*8, parameter :: theta = 1.0D0

    ! Equation coefficients
    real*8 Du, Dv, DTu, DTv
    real*8 Fe, Fw, Fn, Fs, delta_F
    real*8 aw, ae, as, an, ap, ap_o
    real*8 aw_s, ae_s, an_s, as_s, ap_s, ap_prime
    real*8 w_node, e_node, n_node, s_node, p_past
    real*8 bmax, p_err, p_new, b_mean, p_mean
    real*8 T_face, buoy
    real*8 DC_e, DC_w, DC_n, DC_s, S_dc
    
    ! Initial condition generation
    real*8 x_coord, y_coord, pi
    
    ! For calculating program time
    real*8 start_time, end_time
    
    ! For outputting solution
    character(len=60) :: vtk_filename

    ! Flow field matrices
    real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), duu(:,:), un(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vold(:,:), dvv(:,:), vn(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:), bprime(:,:)
    real*8, allocatable :: T(:,:), Tstar(:,:), Told(:,:), Tn(:,:)
    
! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1), duu(0:nx+1, 0:ny+1), un(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1), dvv(0:nx+1, 0:ny+1), vn(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1))
    allocate(T(0:nx+1, 0:ny+1), Tstar(0:nx+1, 0:ny+1), Told(0:nx+1, 0:ny+1), Tn(0:nx+1, 0:ny+1))

! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    volume = dx * dy
    
    ap_o = rho * volume / delta_t
    mu = sqrt(Pr / Ra)
    thermal_alpha = 1.0D0 / sqrt(Ra * Pr)
    
    Du  = mu * dy / dx         
    Dv  = mu * dx / dy
    DTu = rho * thermal_alpha * dy / dx
    DTv = rho * thermal_alpha * dx / dy

! Initialization
    u = 0.0D0
    v = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0
    pstar = 0.0D0
    uold = 0.0D0
    vold = 0.0D0
    un = 0.0D0
    vn = 0.0D0
    duu = 0.0D0
    dvv = 0.0D0
    pprime = 0.0D0
    bprime = 0.0D0
    
    pi = acos(-1.0D0)
    
    ! For Rayleigh-Benard initial temperature
    do j = 0, ny+1
        y_coord = (dble(j) - 0.5D0) * dy
        do i = 0, nx+1
            x_coord = (dble(i) - 0.5D0) * dx
            T(i, j) = T_hot - (T_hot - T_cold) * (y_coord / Ly)
            T(i, j) = T(i, j) + 0.05D0 * (T_hot - T_cold) * cos(3.0D0 * pi * x_coord / Lx)
            
            if (T(i,j) < T_cold) T(i,j) = T_cold
            if (T(i,j) > T_hot)  T(i,j) = T_hot
        enddo
    enddo
    
    Tstar = T
    Told = T
    Tn = T
    
    ! Initial boundary conditions
    !$omp parallel do 
    do j = 0, ny+1
        T(0, j) = T(1, j)                  ! Left adiabatic
        T(nx+1, j) = T(nx, j)              ! Right adiabatic
    enddo
    !$omp end parallel do
    !$omp parallel do
    do i = 0, nx+1
        T(i, 0) = 2.0D0 * T_hot - T(i, 1)  ! Bottom wall hot
        T(i, ny+1) = 2.0D0 * T_cold - T(i, ny) ! Top wall cold
    enddo
    !$omp end parallel do
    
! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step
        
        ! Save old fields
        !$omp parallel do collapse(2)
        do j = 0, ny+1
            do i = 0, nx+1
                un(i, j) = u(i, j)
                vn(i, j) = v(i, j)
                Tn(i, j) = T(i, j)
            enddo
        enddo
        !$omp end parallel do

        ! SIMPLE algo. starts here.
        do iter = 1, maxiter
            
            ! Store old velocity and temperature fields
            !$omp parallel do collapse(2)
            do j = 0, ny+1
                do i = 0, nx+1
                    uold(i, j) = u(i, j)
                    vold(i, j) = v(i, j)
                    Told(i, j) = T(i, j)
                enddo
            enddo
            !$omp end parallel do

            ! Temperature predictor
            !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, aw_s, ae_s, as_s, an_s, ap_s, ap_prime, ap, &
            !$omp w_node, e_node, s_node, n_node, p_past, DC_e, DC_w, DC_n, DC_s, S_dc)
            do j = 1, ny
                do i = 1, nx
                    Fe = rho * dy * u(i+1, j)
                    Fw = rho * dy * u(i, j)
                    Fn = rho * dx * v(i, j+1)
                    Fs = rho * dx * v(i, j)
                    
                    DC_e = 0.0D0
                    DC_w = 0.0D0
                    DC_n = 0.0D0
                    DC_s = 0.0D0
                    
                    if (i > 1 .and. i < nx) then
                        if (Fe > 0.0D0) then
                            DC_e = Fe * (0.375D0 * (Told(i+1, j) - Told(i, j)) - 0.125D0 * (Told(i, j) - Told(i-1, j)))
                        else
                            DC_e = Fe * (0.375D0 * (Told(i, j) - Told(i+1, j)) - 0.125D0 * (Told(i+1, j) - Told(min(i+2, nx+1), j)))
                        endif
                        
                        if (Fw > 0.0D0) then
                            DC_w = Fw * (0.375D0 * (Told(i, j) - Told(i-1, j)) - 0.125D0 * (Told(i-1, j) - Told(max(i-2, 0), j)))
                        else
                            DC_w = Fw * (0.375D0 * (Told(i-1, j) - Told(i, j)) - 0.125D0 * (Told(i, j) - Told(i+1, j)))
                        endif
                    endif
                    
                    if (j > 1 .and. j < ny) then
                        if (Fn > 0.0D0) then
                            DC_n = Fn * (0.375D0 * (Told(i, j+1) - Told(i, j)) - 0.125D0 * (Told(i, j) - Told(i, j-1)))
                        else
                            DC_n = Fn * (0.375D0 * (Told(i, j) - Told(i, j+1)) - 0.125D0 * (Told(i, j+1) - Told(i, min(j+2, ny+1))))
                        endif
                        
                        if (Fs > 0.0D0) then
                            DC_s = Fs * (0.375D0 * (Told(i, j) - Told(i, j-1)) - 0.125D0 * (Told(i, j-1) - Told(i, max(j-2, 0))))
                        else
                            DC_s = Fs * (0.375D0 * (Told(i, j-1) - Told(i, j)) - 0.125D0 * (Told(i, j) - Told(i, j+1)))
                        endif
                    endif
                    
                    S_dc = DC_w - DC_e + DC_s - DC_n
                    
                    aw_s = DTu + max(Fw, 0.0D0)
                    ae_s = DTu + max(-Fe, 0.0D0)
                    as_s = DTv + max(Fs, 0.0D0)
                    an_s = DTv + max(-Fn, 0.0D0)
                    
                    ap_s = aw_s + ae_s + as_s + an_s 
                    ap = (ap_o + theta * ap_s) / alph_T
                    ap_prime = ap_o - (1.0D0 - theta) * ap_s
                    
                    w_node = aw_s * (theta * Told(i-1, j) + (1.0D0 - theta) * Tn(i-1, j))
                    e_node = ae_s * (theta * Told(i+1, j) + (1.0D0 - theta) * Tn(i+1, j))
                    n_node = an_s * (theta * Told(i, j+1) + (1.0D0 - theta) * Tn(i, j+1))
                    s_node = as_s * (theta * Told(i, j-1) + (1.0D0 - theta) * Tn(i, j-1))
                    p_past = ap_prime * Tn(i, j)
                    
                    Tstar(i, j) = (w_node + e_node + n_node + s_node + p_past + S_dc + &
                                   (1.0D0 - alph_T) * ap * Told(i, j)) / ap
                enddo
            enddo
            !$omp end parallel do
            
            ! Applying boundary conditions for T*
            !$omp parallel do collapse(2)
            do j = 1, ny
                do i = 1, nx
                    T(i, j) = Tstar(i, j)
                enddo
            enddo
            !$omp end parallel do
            !$omp parallel do
            do j = 0, ny+1
                T(0, j) = T(1, j)           ! Left adiabatic
                T(nx+1, j) = T(nx, j)       ! Right adiabatic
            enddo
            !$omp end parallel do
            !$omp parallel do
            do i = 0, nx+1
                T(i, 0) = 2.0D0 * T_hot - T(i, 1)      ! Bottom hot
                T(i, ny+1) = 2.0D0 * T_cold - T(i, ny) ! Top cold
            enddo
            !$omp end parallel do

            ! u-momentum predictor 
            !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, aw_s, ae_s, as_s, an_s, ap_s, ap_prime, ap, &
            !$omp w_node, e_node, s_node, n_node, p_past, DC_e, DC_w, DC_n, DC_s, S_dc)
            do j = 1, ny
                do i = 2, nx
                    
                    Fe = 0.5D0 * rho * dy * (u(i, j) + u(i+1, j))
                    Fw = 0.5D0 * rho * dy * (u(i-1, j) + u(i, j))
                    Fn = 0.5D0 * rho * dx * (v(i-1, j+1) + v(i, j+1))
                    Fs = 0.5D0 * rho * dx * (v(i-1, j) + v(i, j))
                    
                    DC_e = 0.0D0
                    DC_w = 0.0D0
                    DC_n = 0.0D0
                    DC_s = 0.0D0
                    
                    if (i > 2 .and. i < nx) then
                        if (Fe > 0.0D0) then
                            DC_e = Fe * (0.375D0 * (uold(i+1, j) - uold(i, j)) - 0.125D0 * (uold(i, j) - uold(i-1, j)))
                        else
                            DC_e = Fe * (0.375D0 * (uold(i, j) - uold(i+1, j)) - 0.125D0 * (uold(i+1, j) - uold(min(i+2, nx+1), j)))
                        endif
                        
                        if (Fw > 0.0D0) then
                            DC_w = Fw * (0.375D0 * (uold(i, j) - uold(i-1, j)) - 0.125D0 * (uold(i-1, j) - uold(max(i-2, 0), j)))
                        else
                            DC_w = Fw * (0.375D0 * (uold(i-1, j) - uold(i, j)) - 0.125D0 * (uold(i, j) - uold(i+1, j)))
                        endif
                    endif
                    
                    if (j > 1 .and. j < ny) then
                        if (Fn > 0.0D0) then
                            DC_n = Fn * (0.375D0 * (uold(i, j+1) - uold(i, j)) - 0.125D0 * (uold(i, j) - uold(i, j-1)))
                        else
                            DC_n = Fn * (0.375D0 * (uold(i, j) - uold(i, j+1)) - 0.125D0 * (uold(i, j+1) - uold(i, min(j+2, ny+1))))
                        endif
                        
                        if (Fs > 0.0D0) then
                            DC_s = Fs * (0.375D0 * (uold(i, j) - uold(i, j-1)) - 0.125D0 * (uold(i, j-1) - uold(i, max(j-2, 0))))
                        else
                            DC_s = Fs * (0.375D0 * (uold(i, j-1) - uold(i, j)) - 0.125D0 * (uold(i, j) - uold(i, j+1)))
                        endif
                    endif
                    
                    S_dc = DC_w - DC_e + DC_s - DC_n
                    
                    aw_s = Du + max(Fw, 0.0D0)
                    ae_s = Du + max(-Fe, 0.0D0)
                    as_s = Dv + max(Fs, 0.0D0)
                    an_s = Dv + max(-Fn, 0.0D0)
                    
                    ap_s = aw_s + ae_s + as_s + an_s 
                    ap = (ap_o + theta * ap_s) / alph_u
                    ap_prime = ap_o - (1.0D0 - theta) * ap_s
                    duu(i, j) = dy / ap
                    
                    w_node = aw_s * (theta * uold(i-1, j) + (1.0D0 - theta) * un(i-1, j))
                    e_node = ae_s * (theta * uold(i+1, j) + (1.0D0 - theta) * un(i+1, j))
                    n_node = an_s * (theta * uold(i, j+1) + (1.0D0 - theta) * un(i, j+1))
                    s_node = as_s * (theta * uold(i, j-1) + (1.0D0 - theta) * un(i, j-1))
                    p_past = ap_prime * un(i, j)
                    
                    ustar(i, j) = (w_node + e_node + n_node + s_node + p_past + S_dc + &
                                   (1.0D0 - alph_u) * ap * uold(i, j) + &
                                   (pstar(i-1, j) - pstar(i, j)) * dy) / ap
                enddo
            enddo
            !$omp end parallel do

            ! v-momentum predictor
            !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, aw_s, ae_s, as_s, an_s, ap_s, ap_prime, ap, &
            !$omp w_node, e_node, s_node, n_node, p_past, T_face, buoy, DC_e, DC_w, DC_n, DC_s, S_dc)
            do j = 2, ny
                do i = 1, nx
                    Fe = 0.5D0 * rho * dy * (u(i+1, j-1) + u(i+1, j))
                    Fw = 0.5D0 * rho * dy * (u(i, j-1) + u(i, j))
                    Fn = 0.5D0 * rho * dx * (v(i, j) + v(i, j+1))
                    Fs = 0.5D0 * rho * dx * (v(i, j-1) + v(i, j))
                    
                    DC_e = 0.0D0
                    DC_w = 0.0D0
                    DC_n = 0.0D0
                    DC_s = 0.0D0
                    
                    if (i > 1 .and. i < nx) then
                        if (Fe > 0.0D0) then
                            DC_e = Fe * (0.375D0 * (vold(i+1, j) - vold(i, j)) - 0.125D0 * (vold(i, j) - vold(i-1, j)))
                        else
                            DC_e = Fe * (0.375D0 * (vold(i, j) - vold(i+1, j)) - 0.125D0 * (vold(i+1, j) - vold(min(i+2, nx+1), j)))
                        endif
                        
                        if (Fw > 0.0D0) then
                            DC_w = Fw * (0.375D0 * (vold(i, j) - vold(i-1, j)) - 0.125D0 * (vold(i-1, j) - vold(max(i-2, 0), j)))
                        else
                            DC_w = Fw * (0.375D0 * (vold(i-1, j) - vold(i, j)) - 0.125D0 * (vold(i, j) - vold(i+1, j)))
                        endif
                    endif
                    
                    if (j > 2 .and. j < ny) then
                        if (Fn > 0.0D0) then
                            DC_n = Fn * (0.375D0 * (vold(i, j+1) - vold(i, j)) - 0.125D0 * (vold(i, j) - vold(i, j-1)))
                        else
                            DC_n = Fn * (0.375D0 * (vold(i, j) - vold(i, j+1)) - 0.125D0 * (vold(i, j+1) - vold(i, min(j+2, ny+1))))
                        endif
                        
                        if (Fs > 0.0D0) then
                            DC_s = Fs * (0.375D0 * (vold(i, j) - vold(i, j-1)) - 0.125D0 * (vold(i, j-1) - vold(i, max(j-2, 0))))
                        else
                            DC_s = Fs * (0.375D0 * (vold(i, j-1) - vold(i, j)) - 0.125D0 * (vold(i, j) - vold(i, j+1)))
                        endif
                    endif
                    
                    S_dc = DC_w - DC_e + DC_s - DC_n

                    aw_s = Du + max(Fw, 0.0D0)
                    ae_s = Du + max(-Fe, 0.0D0)
                    as_s = Dv + max(Fs, 0.0D0)
                    an_s = Dv + max(-Fn, 0.0D0)
                    
                    ap_s = aw_s + ae_s + as_s + an_s 
                    ap = (ap_o + theta * ap_s) / alph_v
                    ap_prime = ap_o - (1.0D0 - theta) * ap_s
                    dvv(i, j) = dx / ap
                    
                    w_node = aw_s * (theta * vold(i-1, j) + (1.0D0 - theta) * vn(i-1, j))
                    e_node = ae_s * (theta * vold(i+1, j) + (1.0D0 - theta) * vn(i+1, j))
                    n_node = an_s * (theta * vold(i, j+1) + (1.0D0 - theta) * vn(i, j+1))
                    s_node = as_s * (theta * vold(i, j-1) + (1.0D0 - theta) * vn(i, j-1))
                    p_past = ap_prime * vn(i, j)
                    
                    T_face = 0.5D0 * (T(i, j-1) + T(i, j))
                    buoy = (T_face - T_ref) * volume
                    
                    vstar(i, j) = (w_node + e_node + n_node + s_node + p_past + S_dc + &
                                   (1.0D0 - alph_v) * ap * vold(i, j) + buoy + &
                                   (pstar(i, j-1) - pstar(i, j)) * dx) / ap
                enddo
            enddo
            !$omp end parallel do

            ! Applying boundary conditions for u* and v*
            !$omp parallel do
            do j = 1, ny
                ustar(1, j) = 0.0D0
                ustar(nx+1, j) = 0.0D0
                vstar(0, j) = -vstar(1, j) 
                vstar(nx+1, j) = -vstar(nx, j)
            enddo
            !$omp end parallel do
            !$omp parallel do
            do i = 1, nx
                ustar(i, 0) = -ustar(i, 1)
                ustar(i, ny+1) = -ustar(i, ny)
                vstar(i, 1) = 0.0D0
                vstar(i, ny+1) = 0.0D0
            enddo
            !$omp end parallel do

            ! For mass source
            bmax = 0.0D0
            b_mean = 0.0D0
            !$omp parallel do collapse(2) reduction(max:bmax) reduction(+:b_mean)
            do j = 1, ny
                do i = 1, nx
                    bprime(i, j) = (rho * dy * ustar(i, j)) - (rho * dy * ustar(i+1, j)) &
                                 + (rho * dx * vstar(i, j)) - (rho * dx * vstar(i, j+1))
                    
                    b_mean = b_mean + bprime(i, j)
                    
                    if (abs(bprime(i, j)) > bmax) bmax = abs(bprime(i, j))
                enddo
            enddo
            !$omp end parallel do
            
            ! Mean mass subtraction
            b_mean = b_mean / dble(nx * ny)
            !$omp parallel do collapse(2)
            do j = 1, ny
                do i = 1, nx
                    bprime(i, j) = bprime(i, j) - b_mean
                enddo
            enddo
            !$omp end parallel do

            ! Red-Black cells SOR
            !$omp parallel do collapse(2)
            do j = 0, ny+1
                do i = 0, nx+1
                    pprime(i, j) = 0.0D0
                enddo
            enddo
            !$omp end parallel do

            do p_iter = 1, 500
                p_err = 0.0D0
                !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do j = 1, ny
                    do i = 1 + mod(j+1, 2), nx, 2
                        
                        aw = rho * dy * duu(i, j)
                        ae = rho * dy * duu(i+1, j)
                        as = rho * dx * dvv(i, j)
                        an = rho * dx * dvv(i, j+1)
                        
                        if (i == 1)  aw = 0.0D0
                        if (i == nx) ae = 0.0D0
                        if (j == 1)  as = 0.0D0
                        if (j == ny) an = 0.0D0
                        ap = aw + ae + as + an
                        
                        if (ap > 1.0D-15) then
                            p_new = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + an * pprime(i, j+1) + as * pprime(i, j-1) &
                                    + bprime(i, j)) / ap
                            
                            p_new = (1.0D0 - omega) * pprime(i, j) + omega * p_new
                            
                            if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                            
                            pprime(i, j) = p_new
                        endif
                    enddo
                enddo
                !$omp end parallel do
                
                !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do j = 1, ny
                    do i = 1 + mod(j, 2), nx, 2
                        
                        aw = rho * dy * duu(i, j)
                        ae = rho * dy * duu(i+1, j)
                        as = rho * dx * dvv(i, j)
                        an = rho * dx * dvv(i, j+1)
                        
                        if (i == 1)  aw = 0.0D0
                        if (i == nx) ae = 0.0D0
                        if (j == 1)  as = 0.0D0
                        if (j == ny) an = 0.0D0
                        ap = aw + ae + as + an
                        
                        if (ap > 1.0D-15) then
                            p_new = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + an * pprime(i, j+1) + as * pprime(i, j-1) &
                                   + bprime(i, j)) / ap
                            
                            p_new = (1.0D0 - omega) * pprime(i, j) + omega * p_new
                            
                            if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                            
                            pprime(i, j) = p_new
                        endif
                    enddo
                enddo
                !$omp end parallel do
                
                if (p_err < sor_tol) exit
            enddo
            
            ! Symmetrical matrix
            p_mean = 0.0D0
            !$omp parallel do collapse(2) reduction(+:p_mean)
            do j = 1, ny
                do i = 1, nx
                    p_mean = p_mean + pprime(i, j)
                enddo
            enddo
            !$omp end parallel do
            p_mean = p_mean / dble(nx * ny)
            !$omp parallel do collapse(2)
            do j = 1, ny
                do i = 1, nx
                    pprime(i, j) = pprime(i, j) - p_mean
                enddo
            enddo
            !$omp end parallel do

            ! Correcting velocity and pressure fields
            !$omp parallel do collapse(2) private(i, j)
            do j = 1, ny
                do i = 1, nx
                    p(i, j) = pstar(i, j) + alph_p * pprime(i, j)
                    pstar(i, j) = p(i, j)
                    
                    if (i > 1) then
                        u(i, j) = ustar(i, j) + duu(i, j) * (pprime(i-1, j) - pprime(i, j))
                    endif
        
                    if (j > 1) then
                        v(i, j) = vstar(i, j) + dvv(i, j) * (pprime(i, j-1) - pprime(i, j))
                    endif
                enddo
            enddo
            !$omp end parallel do

            ! Applying boundary conditions for corrected u and v
            !$omp parallel do
            do j = 1, ny
                u(1, j) = 0.0D0 
                u(nx+1, j) = 0.0D0
                v(0, j) = -v(1, j) 
                v(nx+1, j) = -v(nx, j)
            enddo
            !$omp end parallel do
            !$omp parallel do
            do i = 1, nx
                u(i, 0) = -u(i, 1)
                u(i, ny+1) = -u(i, ny)
                v(i, 1) = 0.0D0
                v(i, ny+1) = 0.0D0
            enddo
            !$omp end parallel do
            
            if (bmax <= tol) exit
        enddo

        write(*,'(A,I6,A,I4,A,ES12.4)') ' Step:', time_step, '  Inner iters:', iter, '  Max mass err:', bmax
            
        if (mod(time_step, output_freq) == 0) then
            write(vtk_filename, '(A, I6.6, A)') 'RBC_2D_', time_step, '.vtk'
            call write_vtk_rbc(trim(vtk_filename), u, v, p, T, nx, ny, dx, dy)
        endif
    enddo

! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) 'Wall time (s) : ', end_time - start_time

! Deallocate memory
    deallocate(u, ustar, uold, duu, un)
    deallocate(v, vstar, vold, dvv, vn)
    deallocate(p, pstar, pprime, bprime)
    deallocate(T, Tstar, Told, Tn)
    
    stop
end program Rayleigh_Benard_Convection

    
! This subroutine aims at writing VTK output files
subroutine write_vtk_rbc(filename, u, v, p, T, nx, ny, dx, dy)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny
    real*8,    intent(in) :: dx, dy
    real*8,    intent(in) :: u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), p(0:nx+1, 0:ny+1), T(0:nx+1, 0:ny+1)
    integer*4 :: i, j
    
    open(unit=20, file=filename, status='unknown')
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') '2D Rayleigh-Benard Convection (QUICK)'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET STRUCTURED_POINTS'
    write(20, '(A, I0, 1X, I0, 1X, I0)')          'DIMENSIONS ', nx, ny, 1
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ',     0.0D0, 0.0D0, 0.0D0
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'SPACING ',    dx, dy, 1.0D0
    write(20, '(A, I0)')                          'POINT_DATA ', nx * ny
    write(20, '(A)') 'SCALARS Pressure float 1'
    write(20, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(20, '(5E15.6)') (p(i, j), i = 1, nx)
    enddo
    write(20, '(A)') 'SCALARS Temperature float 1'
    write(20, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(20, '(5E15.6)') (T(i, j), i = 1, nx)
    enddo
    write(20, '(A)') 'VECTORS Velocity float'
    do j = 1, ny
        do i = 1, nx
            write(20, '(3E15.6)') 0.5D0 * (u(i, j) + u(i+1, j)), 0.5D0 * (v(i, j) + v(i, j+1)), 0.0D0
        enddo
    enddo
    close(20)
end subroutine write_vtk_rbc