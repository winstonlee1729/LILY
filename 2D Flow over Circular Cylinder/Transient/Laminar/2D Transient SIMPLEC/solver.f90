program Transient_Cylinder_SIMPLEC
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
    parameter(dt = 0.05D0)
    
    ! Cylinder geometry 
    parameter(cx = 5.0D0, cy = 5.0D0, R = 0.5D0) 
    
    ! Under-relaxation
    real*8 alph_p, alph_u, alph_v
    parameter(alph_p = 0.5D0, alph_u = 0.5D0, alph_v = 0.5D0)
    
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
    real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), du(:,:), un(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vold(:,:), dv(:,:), vn(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:), bprime(:,:)
    real*8, allocatable :: temp_pprime(:,:)
    
! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1), du(0:nx+1, 0:ny+1), un(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1), vn(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1))
    allocate(temp_pprime(0:nx+1, 0:ny+1))
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    volume = dx * dy * thickness
    
    Dw = mu * (dy / dx)
    De = Dw
    Dn = mu * (dx / dy) 
    Ds = Dn
    
    ap_o = rho * (volume / dt)
    
! Initialization 
    u = 0.0D0
    v = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0 
    pstar = 0.0D0
    un = 0.0D0
    vn = 0.0D0
    ! Applying Initial inlet boundary condition
    do j = 0, ny+1
        u(1, j) = U_inf
        un(1, j) = U_inf
    enddo
    
! Open output file
    open(unit=10, file='cylinder_transient_SIMPLEC.txt', status='replace')
    
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
    
        ! SIMPLEC alogrithm starts
        do iter = 1, maxiter
            
            ! Save old velocities
            !$omp parallel do collapse(2)
            do i = 0, nx+1
                do j = 0, ny+1
                    uold(i, j) = u(i, j)
                    vold(i, j) = v(i, j)
                enddo
            enddo
            !$omp end parallel do
            
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
                        
                        du(i, j) = dy / (ap/alph_u - theta * (aw_s + ae_s + an_s + as_s))
                        
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
                        
                        dv(i, j) = dx / (ap/alph_v - theta * (aw_s + ae_s + as_s + an_s))
                        
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
            do j = 0, ny+1
                ustar(1, j) = U_inf                 ! Inlet
                vstar(0, j) = 0.0D0                 ! v = 0 at inlet
                vstar(nx+1, j) = vstar(nx, j)       ! Zero gradient v at outlet
            enddo
            
            ! Mass conservation at outlet for ustar
            M_in = 0.0D0 
            M_out = 0.0D0
            do j = 1, ny
                M_in = M_in + ustar(1, j) * dy
                M_out = M_out + ustar(nx, j) * dy
            enddo
            
            mass_ratio = M_in / max(M_out, 1.0D-10)
            
            do j = 1, ny
                ustar(nx+1, j) = ustar(nx, j) * mass_ratio
            enddo
            
            ! Top and bottom 
            do i = 0, nx+1
                vstar(i, 1) = 0.0D0                 ! Bottom wall v = 0
                vstar(i, ny+1) = 0.0D0              ! Top wall v = 0
                ustar(i, 0) = -ustar(i, 1)
                ustar(i, ny+1) = -ustar(i, ny)       
            enddo

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
            
            
            ! Pressure correction equation
            do k = 1, 100
                !$omp parallel do private(j, aw, ae, as, an, ap, ap_p)
                do i = 1, nx
                    do j = 1, ny
                        aw = rho * dy * du(i, j)
                        ae = rho * dy * du(i+1, j)
                        as = rho * dx * dv(i, j) 
                        an = rho * dx * dv(i, j+1)
                        
                        if (i == 1) aw = 0.0D0
                        if (i == nx) ae = 0.0D0
                        if (j == 1) as = 0.0D0 
                        if (j == ny) an = 0.0D0
                        ap_p = aw + ae + an + as
                        
                        ! Prevent division by zero inside the cylinder
                        if (ap_p < 1.0D-15) then
                            temp_pprime(i, j) = 0.0D0
                        else
                            temp_pprime(i, j) = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + &
                                                an * pprime(i, j+1) + as * pprime(i, j-1) + &
                                                bprime(i, j)) / ap_p
                        endif
                    enddo
                enddo
                !$omp end parallel do
                
                !$omp parallel do collapse(2)
                do i = 1, nx
                    do j = 1, ny
                        pprime(i, j) = temp_pprime(i, j)
                    enddo
                enddo
                !$omp end parallel do
            enddo
            
            
            ! Field updates
            !$omp parallel do private(j, uc, vc)
            do i = 1, nx
                do j = 1, ny
                    p(i, j) = pstar(i, j) + alph_p * pprime(i, j)
                    pstar(i, j)= p(i, j)
                    
                    if (i > 1) then
                        uc = ustar(i, j) + du(i, j) * (pprime(i-1, j) - pprime(i, j))
                        u(i, j) = alph_u * uc + (1.0D0 - alph_u) * uold(i, j)
                    endif
                    
                    if (j > 1) then
                        vc = vstar(i, j) + dv(i, j) * (pprime(i, j-1) - pprime(i, j))
                        v(i, j) = alph_v * vc + (1.0D0 - alph_v) * vold(i, j)
                    endif
                enddo
            enddo
            !$omp end parallel do
    
            ! Open boundary conditions and mass conservation for u and v
            ! Inlet and outlet
            do j = 0, ny+1
                u(1, j) = U_inf                     ! Inlet
                v(0, j) = 0.0D0                     ! v = 0 at inlet
                v(nx+1, j) = v(nx, j)               ! Zero gradient v at outlet
            enddo
            
            ! Mass conservation at outlet for ustar
            M_in = 0.0D0 
            M_out = 0.0D0
            do j = 1, ny
                M_in = M_in + u(1, j) * dy
                M_out = M_out + u(nx, j) * dy
            enddo
            
            mass_ratio = M_in / max(M_out, 1.0D-10)
            
            do j = 1, ny
                u(nx+1, j) = u(nx, j) * mass_ratio
            enddo
            
            ! Top and bottom 
            do i = 0, nx+1
                v(i, 1) = 0.0D0                 ! Bottom wall v = 0
                v(i, ny+1) = 0.0D0              ! Top wall v = 0
                u(i, 0) = -ustar(i, 1)
                u(i, ny+1) = -ustar(i, ny)       
            enddo
            
            if (bmax <= tol) exit
            
        enddo  ! End of SIMPLEC algo.
        
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

! Deallocate memory
    deallocate(u, ustar, uold, du, un)
    deallocate(v, vstar, vold, dv, vn)
    deallocate(p, pstar, pprime, bprime, temp_pprime)

    stop
end program Transient_Cylinder_SIMPLEC