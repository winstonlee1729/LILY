program Transient_Cylinder_PISO
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration
    integer*4 nx, ny
    integer*4 i, j, k
    integer*4 iter, maxiter
    integer*4 time_step, max_step, output_freq
    parameter(nx = 200, ny = 100)
    parameter(maxiter = 5000, max_step = 4000)
    parameter(output_freq = 50)
    
    ! Tolerance for convergence
    real*8, parameter :: tol = 1.0D-5
    
    ! Parameters for time calculation
    real*8 start_time, end_time
    
    ! Flow field parameters
    real*8 rho, mu, U_inf
    real*8 Lx, Ly, thickness, volume, dx, dy, dt
    real*8 cx, cy, R, rad_sq
    ! Re = rho * U_inf * (2R) / mu 
    parameter(rho = 1.0D0, mu = 0.01D0, U_inf = 1.0D0)
    parameter(Lx = 20.0D0, Ly = 10.0D0, thickness = 1.0D0)
    parameter(dt = 0.1D0)
    
    ! Cylinder geometry 
    parameter(cx = 5.0D0, cy = 5.0D0, R = 0.5D0) 
    
    ! Under-relaxation
    real*8 alph_p, alph_u, alph_v
    parameter(alph_p = 1.0D0, alph_u = 1.0D0, alph_v = 1.0D0)
    
    ! SOR Acceleration
    real*8, parameter :: omega = 1.3D0
    real*8 p_err, p_new
    
    ! Equation coefficients for pressure correction equation
    real*8 aw, ae, an, as, ap, ap_o, delta_F, ap_p
    
    ! Equation coefficients for momentum predictor
    real*8 aw_s, ae_s, an_s, as_s, ap_s, ap_prime
    real*8 w_node, e_node, n_node, s_node, p_past
    real*8 Fw, Fe, Fn, Fs
    real*8 Dw, De, Dn, Ds, bmax, uc, vc
    
    ! For outputting solution
    real*8 xc, yc, pc
    
    ! For mass conservation
    real*8 M_in, M_out, mass_ratio
    
    ! Theta for temporal discretization scheme
    ! 1.0 = Fully-implicit | 0.5 = Crank-Nicolson | 0.0 = Explicit
    real*8, parameter :: theta = 1.0D0
    
    ! Field matrix
    real*8, allocatable :: u(:,:), ustar(:,:), ustarstar(:,:), uold(:,:), un(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vstarstar(:,:), vold(:,:), vn(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pstarstar(:,:)
    
    real*8, allocatable :: pprime(:,:), pprime2(:,:)
    real*8, allocatable :: bprime(:,:), bprime2(:,:)
    real*8, allocatable :: uprime(:,:), vprime(:,:)
    real*8, allocatable :: ucorr(:,:), vcorr(:,:)
    real*8, allocatable :: du(:,:), dv(:,:)
    
! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), ustarstar(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1), un(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vstarstar(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1), vn(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pstarstar(0:nx+1, 0:ny+1))
    allocate(pprime(0:nx+1, 0:ny+1), pprime2(0:nx+1, 0:ny+1))
    allocate(bprime(0:nx+1, 0:ny+1), bprime2(0:nx+1, 0:ny+1))
    allocate(uprime(0:nx+1, 0:ny+1), vprime(0:nx+1, 0:ny+1))
    allocate(ucorr(0:nx+1, 0:ny+1), vcorr(0:nx+1, 0:ny+1))
    allocate(du(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    volume = dx * dy * thickness
    
    Dw = mu * (dy / dx)
    De = mu * (dy / dx)
    Dn = mu * (dx / dy) 
    Ds = mu * (dx / dy) 
    
    ap_o = rho * (volume / dt)
    
! Initialization
    u = 0.0D0
    v = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0 
    pstar = 0.0D0
    ustarstar = 0.0D0
    vstarstar = 0.0D0
    pstarstar = 0.0D0
    un = 0.0D0
    vn = 0.0D0
    ! Applying initial inlet boundary condition
    do j = 0, ny+1
        u(1, j) = U_inf
        un(1, j) = U_inf
    enddo
    
! Open output file
    open(unit=10, file='cylinder_transient_PISO.txt', status='replace')
    
! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step
        
        ! Save the converged fields from the previous time step
        !$omp parallel do collapse(2)
        do i = 0, nx+1
            do j = 0, ny+1
                un(i, j) = u(i, j)
                vn(i, j) = v(i, j)
            enddo
        enddo
        !$omp end parallel do
    
        ! PISO alogrithm starts
        do iter = 1, maxiter
            
            ! Save old velocities and reset correction matrices
            !$omp parallel do collapse(2)
            do i = 0, nx+1
                do j = 0, ny+1
                    uold(i, j) = u(i, j)
                    vold(i, j) = v(i, j)
                    uprime(i, j) = 0.0D0
                    vprime(i, j) = 0.0D0
                    ucorr(i, j) = 0.0D0
                    vcorr(i, j) = 0.0D0
                enddo
            enddo
            !$omp end parallel do
            
            ! Predictor step
            ! u-momentum predictor 
            !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw_s, ae_s, as_s, an_s, &
            !$omp ap_s, ap_prime, ap, w_node, e_node, s_node, n_node, p_past, rad_sq)
            do i = 2, nx
                do j = 1, ny
                    ! Define cylinder location
                    rad_sq = ((dble(i) - 1.0D0) * dx - cx)**2 + ((dble(j) - 0.5D0) * dy - cy)**2
                    
                    if (rad_sq <= R**2) then        ! Inside the cylinder
                        ap = 1.0D30  
                        du(i, j) = 0.0D0
                        ustar(i, j) = 0.0D0
                    else                            ! Inside the flow field
                        Fe = 0.5D0 * rho * dy * (u(i, j) + u(i+1, j))
                        Fw = 0.5D0 * rho * dy * (u(i-1, j) + u(i, j))
                        Fn = 0.5D0 * rho * dx * (v(i-1, j+1) + v(i, j+1))
                        Fs = 0.5D0 * rho * dx * (v(i-1, j) + v(i, j))
                        
                        delta_F = Fe - Fw + Fn - Fs
                        aw_s = Dw + max(Fw, 0.0D0) 
                        ae_s = De + max(-Fe, 0.0D0)
                        as_s = Ds + max(Fs, 0.0D0) 
                        an_s = Dn + max(-Fn, 0.0D0)  
                        
                        ap_s = aw_s + ae_s + as_s + an_s + delta_F
                        ap = ap_o + theta * ap_s
                        ap_prime = ap_o - (1.0D0 - theta) * ap_s
                        
                        du(i, j) = dy / ap
                        
                        w_node = aw_s * (theta * uold(i-1, j) + (1.0D0 - theta) * un(i-1, j))
                        e_node = ae_s * (theta * uold(i+1, j) + (1.0D0 - theta) * un(i+1, j))
                        s_node = as_s * (theta * uold(i, j-1) + (1.0D0 - theta) * un(i, j-1))
                        n_node = an_s * (theta * uold(i, j+1) + (1.0D0 - theta) * un(i, j+1))
                        p_past = ap_prime * un(i, j)
                        
                        ustar(i,j) = (w_node + e_node + s_node + n_node + p_past + &
                                     (pstar(i-1, j) - pstar(i, j)) * dy) / ap
                    endif
                enddo
            enddo
            !$omp end parallel do
            

            ! v-momentum predictor
            !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw_s, ae_s, as_s, an_s, &
            !$omp ap_s, ap_prime, ap, w_node, e_node, s_node, n_node, p_past, rad_sq)
            do i = 1, nx
                do j = 2, ny
                    ! Define cylinder location
                    rad_sq = ((dble(i) - 0.5D0) * dx - cx)**2 + ((dble(j) - 1.0D0) * dy - cy)**2
                    
                    if (rad_sq <= R**2) then            ! Insider the cylinder
                        ap = 1.0D30  
                        dv(i, j) = 0.0D0
                        vstar(i, j) = 0.0D0
                    else                                ! Inside the flow field
                        Fe = 0.5D0 * rho * dy * (u(i+1, j-1) + u(i+1, j))
                        Fw = 0.5D0 * rho * dy * (u(i, j-1) + u(i, j))
                        Fn = 0.5D0 * rho * dx * (v(i, j) + v(i, j+1))
                        Fs = 0.5D0 * rho * dx * (v(i, j-1) + v(i, j))
                        
                        delta_F = Fe - Fw + Fn - Fs
                        aw_s = Dw + max(Fw, 0.0D0)
                        ae_s = De + max(-Fe, 0.0D0)
                        as_s = Ds + max(Fs, 0.0D0)
                        an_s = Dn + max(-Fn, 0.0D0)
                        
                        ap_s = aw_s + ae_s + as_s + an_s + delta_F
                        ap = ap_o + theta * ap_s
                        ap_prime = ap_o - (1.0D0 - theta) * ap_s
                        
                        dv(i, j) = dx / ap
                        
                        w_node = aw_s * (theta * vold(i-1, j) + (1.0D0 - theta) * vn(i-1, j))
                        e_node = ae_s * (theta * vold(i+1, j) + (1.0D0 - theta) * vn(i+1, j))
                        s_node = as_s * (theta * vold(i, j-1) + (1.0D0 - theta) * vn(i, j-1))
                        n_node = an_s * (theta * vold(i, j+1) + (1.0D0 - theta) * vn(i, j+1))
                        p_past = ap_prime * vn(i, j)
                        
                        vstar(i, j) = (w_node + e_node + s_node + n_node + p_past + &
                                     (pstar(i, j-1) - pstar(i, j)) * dx) / ap
                    endif
                enddo
            enddo
            !$omp end parallel do
            

            ! Open boundary conditions and mass conservation for u* and v*
            ! Inlet and outlet
            !$omp parallel do
            do j = 0, ny+1
                ustar(1, j) = U_inf                 ! Inlet
                vstar(0, j) = 0.0D0                 ! v = 0 at inlet
                vstar(nx+1, j) = vstar(nx, j)       ! Zero gradient v at outlet
            enddo
            !$omp end parallel do
            
            ! Mass conservation at outlet for ustar
            M_in = 0.0D0 
            M_out = 0.0D0
            do j = 1, ny
                M_in = M_in + ustar(1, j) * dy
                M_out = M_out + ustar(nx, j) * dy
            enddo
            
            mass_ratio = M_in / max(M_out, 1.0D-10)
            !$omp parallel do
            do j = 1, ny
                ustar(nx+1, j) = ustar(nx, j) * mass_ratio
            enddo
            !$omp end parallel do
            
            ! Top and bottom
            !$omp parallel do
            do i = 0, nx+1
                vstar(i, 1) = 0.0D0                 ! Bottom wall v = 0
                vstar(i, ny+1) = 0.0D0              ! Top wall v = 0
                ustar(i, 0) = -ustar(i, 1)          
                ustar(i, ny+1) = -ustar(i, ny)       
            enddo
            !$omp end parallel do
            
            ! First corrector step
            ! Mass source
            bmax = 0.0D0
            !$omp parallel do private(j) reduction(max:bmax)
            do i = 1, nx
                do j = 1, ny
                    pprime(i, j) = 0.0D0
                    bprime(i, j) = (rho * dy * ustar(i, j)) + (rho * dx * vstar(i, j)) - &
                                  (rho * dy * ustar(i+1, j)) - (rho * dx * vstar(i, j+1))
                    
                    ! Exclude solid cylinder from max error check
                    if (((dble(i) - 0.5) * dx - cx)**2 + ((dble(j) - 0.5) * dy - cy)**2 > R**2) then
                        if (abs(bprime(i, j)) > bmax) bmax = abs(bprime(i, j))
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            ! Red-Black SOR for first pressure correction
            do k = 1, 5000
                p_err = 0.0D0
            
                ! Red nodes
                !$omp parallel do private(j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do i = 1, nx
                    do j = 1 + mod(i, 2), ny, 2
                        aw = rho * dy * du(i, j)
                        ae = rho * dy * du(i+1, j)
                        as = rho * dx * dv(i, j)
                        an = rho * dx * dv(i, j+1)
                        
                        if (i == 1) aw = 0.0D0
                        if (i == nx) ae = 0.0D0
                        if (j == 1) as = 0.0D0
                        if (j == ny) an = 0.0D0
                        ap = aw + ae + an + as
                        
                        if (ap < 1.0D-15) then
                            p_new = 0.0D0
                        else
                            p_new = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + &
                                    an * pprime(i, j+1) + as * pprime(i, j-1) + bprime(i, j)) / ap
                        
                            p_new = (1.0D0 - omega) * pprime(i, j) + omega * p_new
                        endif
                        
                        if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                        pprime(i, j) = p_new
                    enddo
                enddo
                !$omp end parallel do
                
                ! Black nodes
                !$omp parallel do private(j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do i = 1, nx
                    do j = 1 + mod(i+1, 2), ny, 2
                        aw = rho * dy * du(i, j)
                        ae = rho * dy * du(i+1, j)
                        as = rho * dx * dv(i, j)
                        an = rho * dx * dv(i, j+1)
                        
                        if (i == 1) aw = 0.0D0
                        if (i == nx) ae = 0.0D0
                        if (j == 1) as = 0.0D0
                        if (j == ny) an = 0.0D0
                        ap = aw + ae + an + as
                        
                        if (ap < 1.0D-15) then
                            p_new = 0.0D0
                        else
                            p_new = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + &
                                    an * pprime(i, j+1) + as * pprime(i, j-1) + bprime(i, j)) / ap
                        
                            p_new = (1.0D0 - omega) * pprime(i, j) + omega * p_new
                        endif
                        
                        if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                        pprime(i, j) = p_new
                    enddo
                enddo
                !$omp end parallel do
                
                pprime(1, 1) = 0.0D0 ! Pin pressure
                if (p_err < 1.0D-5) exit 
            enddo
            
            ! First field updates
            !$omp parallel do private(j)
            do i = 1, nx
                do j = 1, ny
                    pstarstar(i, j) = pstar(i, j) + pprime(i, j)
                    
                    if (i > 1) then
                        uprime(i, j) = du(i, j) * (pprime(i-1, j) - pprime(i, j))
                        ustarstar(i, j) = ustar(i, j) + uprime(i, j)
                    endif
                    
                    if (j > 1) then
                        vprime(i, j) = dv(i, j) * (pprime(i, j-1) - pprime(i, j))
                        vstarstar(i, j) = vstar(i, j) + vprime(i, j)
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            ! Open boundary conditions and mass conservation for u** and v**
            ! Inlet and outlet
            !$omp parallel do
            do j = 0, ny+1
                ustarstar(1, j) = U_inf                     ! Inlet
                vstarstar(0, j) = 0.0D0                     ! v = 0 at inlet
                vstarstar(nx+1, j) = vstarstar(nx, j)       ! Zero gradient v at outlet
            enddo
            !$omp end parallel do
            
            ! Mass conservation at outlet for ustar
            M_in = 0.0D0 
            M_out = 0.0D0
            do j = 1, ny
                M_in = M_in + ustarstar(1, j) * dy
                M_out = M_out + ustarstar(nx, j) * dy
            enddo
            
            mass_ratio = M_in / max(M_out, 1.0D-10)
            !$omp parallel do
            do j = 1, ny
                ustarstar(nx+1, j) = ustarstar(nx, j) * mass_ratio
            enddo
            !$omp end parallel do
            
            ! Top and bottom
            !$omp parallel do
            do i = 0, nx+1
                vstarstar(i, 1) = 0.0D0                 ! Bottom wall v = 0
                vstarstar(i, ny+1) = 0.0D0              ! Top wall v = 0
                ustarstar(i, 0) = -ustarstar(i, 1)          
                ustarstar(i, ny+1) = -ustarstar(i, ny)       
            enddo
            !$omp end parallel do
            
            ! Second corrector step (ucorr, vcorr, pprime2)
            ! Calculate ucorr
            !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw_s, ae_s, as_s, an_s, ap_s, ap)
            do i = 2, nx
                do j = 1, ny
                    Fe = 0.5D0 * rho * dy * (uold(i, j) + uold(i+1, j))
                    Fw = 0.5D0 * rho * dy * (uold(i-1, j) + uold(i, j))
                    Fn = 0.5D0 * rho * dx * (vold(i-1, j+1) + vold(i, j+1))
                    Fs = 0.5D0 * rho * dx * (vold(i-1, j) + vold(i, j))
                    
                    delta_F = Fe - Fw + Fn - Fs
                    aw_s = Dw + max(Fw, 0.0D0)
                    ae_s = De + max(-Fe, 0.0D0)
                    as_s = Ds + max(Fs, 0.0D0)   
                    an_s = Dn + max(-Fn, 0.0D0)  
                    
                    ap_s = aw_s + ae_s + as_s + an_s + delta_F
                    ap = ap_o + theta * ap_s
                    
                    if (ap < 1.0D-15) then
                        ucorr(i, j) = 0.0D0
                    else
                        ucorr(i, j) = (aw_s * uprime(i-1, j) + ae_s * uprime(i+1, j) + &
                                    as_s * uprime(i, j-1) + an_s * uprime(i, j+1)) / ap
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            ! Calculate vcorr
            !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw_s, ae_s, as_s, an_s, ap_s, ap)
            do i = 1, nx
                do j = 2, ny
                    Fe = 0.5D0 * rho * dy * (uold(i+1, j-1) + uold(i+1, j))
                    Fw = 0.5D0 * rho * dy * (uold(i, j-1) + uold(i, j))
                    Fn = 0.5D0 * rho * dx * (vold(i, j) + vold(i, j+1))
                    Fs = 0.5D0 * rho * dx * (vold(i, j-1) + vold(i, j))
                    
                    delta_F = Fe - Fw + Fn - Fs
                    aw_s = Dw + max(Fw, 0.0D0)
                    ae_s = De + max(-Fe, 0.0D0)
                    as_s = Ds + max(Fs, 0.0D0)
                    an_s = Dn + max(-Fn, 0.0D0)
                    
                    ap_s = aw_s + ae_s + as_s + an_s + delta_F
                    ap = ap_o + theta * ap_s
                    
                    if (ap < 1.0D-15) then
                        vcorr(i, j) = 0.0D0
                    else
                        vcorr(i, j) = (aw_s * vprime(i-1, j) + ae_s * vprime(i+1, j) + &
                                    as_s * vprime(i, j-1) + an_s * vprime(i, j+1)) / ap
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            
            ! Calculate secondary mass source (bprime2)
            !$omp parallel do private(j)
            do i = 1, nx
                do j = 1, ny
                    pprime2(i, j) = 0.0D0
                    bprime2(i, j) = (rho * dy * ucorr(i, j)) + (rho * dx * vcorr(i, j)) &
                                  - (rho * dy * ucorr(i+1, j)) - (rho * dx * vcorr(i, j+1))
                enddo
            enddo
            !$omp end parallel do
            
            ! Red-Black SOR for second pressure correction
            do k = 1, 5000
                p_err = 0.0D0
                
                ! Red nodes
                !$omp parallel do private(j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do i = 1, nx
                    do j = 1 + mod(i, 2), ny, 2
                        aw = rho * dy * du(i, j)
                        ae = rho * dy * du(i+1, j)
                        as = rho * dx * dv(i, j)
                        an = rho * dx * dv(i, j+1)
                        
                        if (i == 1) aw = 0.0D0
                        if (i == nx) ae = 0.0D0
                        if (j == 1) as = 0.0D0
                        if (j == ny) an = 0.0D0
                        ap = aw + ae + an + as
                        
                        if (ap < 1.0D-15) then
                            p_new = 0.0D0
                        else
                            p_new = (ae * pprime2(i+1, j) + aw * pprime2(i-1, j) + &
                                    an * pprime2(i, j+1) + as * pprime2(i, j-1) + bprime2(i, j)) / ap
                        
                            p_new = (1.0D0 - omega) * pprime2(i, j) + omega * p_new
                        endif
                    
                        if (abs(p_new - pprime2(i, j)) > p_err) p_err = abs(p_new - pprime2(i, j))
                        pprime2(i,j) = p_new
                    enddo
                enddo
                !$omp end parallel do
                
                ! Black nodes
                !$omp parallel do private(j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do i = 1, nx
                    do j = 1 + mod(i+1, 2), ny, 2
                        aw = rho * dy * du(i, j)
                        ae = rho * dy * du(i+1, j)
                        as = rho * dx * dv(i, j)
                        an = rho * dx * dv(i, j+1)
                        
                        if (i == 1) aw = 0.0D0
                        if (i == nx) ae = 0.0D0
                        if (j == 1) as = 0.0D0
                        if (j == ny) an = 0.0D0
                        ap = aw + ae + an + as
                        
                        if (ap < 1.0D-15) then
                            p_new = 0.0D0
                        else
                            p_new = (ae * pprime2(i+1, j) + aw * pprime2(i-1, j) + &
                                 an * pprime2(i, j+1) + as * pprime2(i, j-1) + bprime2(i, j)) / ap
                        
                            p_new = (1.0D0 - omega) * pprime2(i, j) + omega * p_new
                        endif
                        
                        if (abs(p_new - pprime2(i, j)) > p_err) p_err = abs(p_new - pprime2(i, j))
                        pprime2(i, j) = p_new
                    enddo
                enddo
                !$omp end parallel do
                
                pprime2(1, 1) = 0.0D0 ! Pin pressure
                if (p_err < 1.0D-5) exit 
            enddo
            
            ! Field updates for the iteration
            !$omp parallel do private(j, uc, vc)
            do i = 1, nx
                do j = 1, ny
                    p(i, j) = pstarstar(i, j) + pprime2(i, j)
                    pstar(i, j) = p(i, j)
                    
                    if (i > 1) then
                        uc = ustarstar(i, j) + ucorr(i, j) + du(i, j) * (pprime2(i-1, j) - pprime2(i, j))
                        u(i, j) = alph_u * uc + (1.0D0 - alph_u) * uold(i, j)
                    endif
                    
                    if (j > 1) then
                        vc = vstarstar(i, j) + vcorr(i, j) + dv(i, j) * (pprime2(i, j-1) - pprime2(i, j))
                        v(i,j) = alph_v * vc + (1.0D0 - alph_v) * vold(i,j)
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            ! Open boundary conditions and mass conservation for u and v
            ! Inlet and outlet
            !$omp parallel do 
            do j = 0, ny+1
                u(1, j) = U_inf                 ! Inlet
                v(0, j) = 0.0D0                 ! v = 0 at inlet
                v(nx+1, j) = v(nx, j)           ! Zero gradient v at outlet
            enddo
            !$omp end parallel do
            
            ! Mass conservation at outlet for ustar
            M_in = 0.0D0 
            M_out = 0.0D0
            do j = 1, ny
                M_in = M_in + u(1, j) * dy
                M_out = M_out + u(nx, j) * dy
            enddo
            
            mass_ratio = M_in / max(M_out, 1.0D-10)
            !$omp parallel do
            do j = 1, ny
                u(nx+1, j) = u(nx, j) * mass_ratio
            enddo
            !$omp end parallel do
            
            ! Top and bottom
            !$omp parallel do
            do i = 0, nx+1
                v(i, 1) = 0.0D0                 ! Bottom wall v = 0
                v(i, ny+1) = 0.0D0              ! Top wall v = 0
                u(i, 0) = -u(i, 1)          
                u(i, ny+1) = -u(i, ny)       
            enddo
            !$omp end parallel do
            
            if (bmax <= tol) exit
            
        enddo  ! End of SIMPLE algo.
        
        write(*,*) 'Time step: ', time_step, '| Inner Iters: ', iter, '| Max mass error: ', bmax
        
        ! Outputting solution
        if (mod(time_step, output_freq) == 0 .or. time_step == max_step) then
            write(10, *) 'ZONE T="Time Step ', time_step, '", I=', nx, ', J=', ny, ', F=POINT'
            do j = 1, ny
                do i = 1, nx
                    xc = (dble(i) - 0.5D0) * dx
                    yc = (dble(j) - 0.5D0) * dy
                    uc = 0.5D0 * (u(i, j) + u(i+1, j))
                    vc = 0.5D0 * (v(i, j) + v(i, j+1)) 
                    pc = p(i, j)
                    write(10, '(5E16.6)') xc, yc, uc, vc, pc
                enddo
            enddo
        endif
        
    enddo ! End of time loop
    
! Stop the clock and calculation running time
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total Time Steps: ', max_step
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'

! Close the output file
    close(10)

! Deallocation memory
    deallocate(u, ustar, ustarstar, uold, un)
    deallocate(v, vstar, vstarstar, vold, vn)
    deallocate(p, pstar, pstarstar)
    deallocate(pprime, pprime2)
    deallocate(bprime, bprime2)
    deallocate(uprime, vprime)
    deallocate(ucorr, vcorr)
    deallocate(du, dv)

    stop
end program Transient_Cylinder_PISO