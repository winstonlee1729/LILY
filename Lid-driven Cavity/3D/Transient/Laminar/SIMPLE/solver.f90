program Transient_Cavity_SIMPLE
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration
    integer*4 nx, ny, nz
    integer*4 i, j, k
    integer*4 iter, maxiter, p_iter
    integer*4 time_step, max_step, output_freq
    parameter(nx = 50, ny = 50, nz = 50)
    parameter(maxiter = 5000, max_step = 1000)
    parameter(output_freq = 20)
    
    ! Tolerance for convergence
    real*8, parameter :: tol = 1.0D-5
    
    ! Parameters for time
    real*8 start_time, end_time
    
    ! Flow field parameters
    real*8 rho, mu, ulid, delta_t
    real*8 Lx, Ly, Lz, volume, dx, dy, dz
    parameter(rho = 1.0D0, mu = 0.01D0, ulid = 5.0D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0)
    parameter(delta_t = 0.1D0)
    
    ! Under-relaxation
    real*8 alph_p, alph_u, alph_v, alph_w
    parameter(alph_p = 0.7D0, alph_u = 0.7D0, alph_v = 0.7D0, alph_w = 0.7D0)
    
    ! Equation coefficients
    real*8 aw, ae, an, as, at, ab, ap, ap_o, ap_p
    real*8 Fw, Fe, Fn, Fs, Ft, Fb, delta_F
    real*8 Dw, De, Dn, Ds, Dt, Db
    real*8 bmax, uc, vc, wc
    
    ! Output string for VTK files
    character(len=50) :: vtk_filename
    
    ! Theta for temporal discretization
    real*8 theta
    real*8 aw_s, ae_s, an_s, as_s, at_s, ab_s, ap_s, ap_prime
    real*8 w_node, e_node, n_node, s_node, t_node, b_node, p_past
    ! 1.0 = Implicit | 0.5 = Crank-Nicolson | 0.0 = Explicit
    parameter(theta = 1.0D0)
    
    ! Field matrix
    real*8, allocatable :: u(:,:,:), ustar(:,:,:), uold(:,:,:), du(:,:,:), un(:,:,:)
    real*8, allocatable :: v(:,:,:), vstar(:,:,:), vold(:,:,:), dv(:,:,:), vn(:,:,:)
    real*8, allocatable :: w(:,:,:), wstar(:,:,:), wold(:,:,:), dww(:,:,:), wn(:,:,:)
    real*8, allocatable :: p(:,:,:), pstar(:,:,:), pprime(:,:,:), bprime(:,:,:)
    real*8, allocatable :: temp_pprime(:,:,:)
    
! Allocate memory
    allocate(u(0:nx+1, 0:ny+1, 0:nz+1), ustar(0:nx+1, 0:ny+1, 0:nz+1), uold(0:nx+1, 0:ny+1, 0:nz+1), du(0:nx+1, 0:ny+1, 0:nz+1), un(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(v(0:nx+1, 0:ny+1, 0:nz+1), vstar(0:nx+1, 0:ny+1, 0:nz+1), vold(0:nx+1, 0:ny+1, 0:nz+1), dv(0:nx+1, 0:ny+1, 0:nz+1), vn(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(w(0:nx+1, 0:ny+1, 0:nz+1), wstar(0:nx+1, 0:ny+1, 0:nz+1), wold(0:nx+1, 0:ny+1, 0:nz+1), dww(0:nx+1, 0:ny+1, 0:nz+1), wn(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(p(0:nx+1, 0:ny+1, 0:nz+1), pstar(0:nx+1, 0:ny+1, 0:nz+1), pprime(0:nx+1, 0:ny+1, 0:nz+1), bprime(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(temp_pprime(0:nx+1, 0:ny+1, 0:nz+1))
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    dz = Lz / dble(nz)
    volume = dx * dy * dz
    
    Dw = mu * (dy * dz / dx)
    De = mu * (dy * dz / dx)
    Dn = mu * (dx * dz / dy)
    Ds = mu * (dx * dz / dy)
    Dt = mu * (dx * dy / dz)
    Db = mu * (dx * dy / dz)
    
    ap_o = rho * (volume / delta_t)
    
! Initialization
    u = 0.0D0
    v = 0.0D0
    w = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0
    wstar = 0.0D0
    pstar = 0.0D0
    un = 0.0D0
    vn = 0.0D0
    wn = 0.0D0
    du = 0.0D0
    dv = 0.0D0
    dww = 0.0D0
    pprime = 0.0D0
    bprime = 0.0D0
    temp_pprime = 0.0D0
    ! Applying boundary condition
    !$omp parallel do collapse(2)
    do i = 0, nx+1
        do j = 0, ny+1
            u(i, j, nz+1) = ulid
        enddo
    enddo
    !$omp end parallel do
    
! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step
        
        ! Save the converged fields from the previous time step
        !$omp parallel do collapse(3)
        do i = 0, nx+1
            do j = 0, ny+1
                do k = 0, nz+1
                    un(i, j, k) = u(i, j, k)
                    vn(i, j, k) = v(i, j, k)
                    wn(i, j, k) = w(i, j, k)
                enddo
            enddo
        enddo
        !$omp end parallel do
    

        ! SIMPLE algo. starts here
        do iter = 1, maxiter
        
            ! Save old velocities
            !$omp parallel do collapse(3)
            do i = 0, nx+1
                do j = 0, ny+1
                    do k = 0, nz+1
                        uold(i, j, k) = u(i, j, k)
                        vold(i, j, k) = v(i, j, k)
                        wold(i, j, k) = w(i, j, k)
                    enddo
                enddo
            enddo
            !$omp end parallel do
         
            ! u-momentum predictor
            !$omp parallel do collapse(3) private(i, j, k, Fe, Fw, Fn, Fs, Ft, Fb, delta_F, aw_s, ae_s, as_s, an_s, at_s, ab_s &
            !$omp , ap_s, ap_prime, ap, w_node, e_node, s_node, n_node, t_node, b_node, p_past)
            do i = 2, nx
                do j = 1, ny
                    do k = 1, nz
                        Fe = 0.5D0 * rho * dy * dz * (u(i, j, k) + u(i+1, j, k))
                        Fw = 0.5D0 * rho * dy * dz * (u(i-1, j, k) + u(i, j, k))
                        Fn = 0.5D0 * rho * dx * dz * (v(i-1, j+1, k) + v(i, j+1, k))
                        Fs = 0.5D0 * rho * dx * dz * (v(i-1, j, k) + v(i, j, k))
                        Ft = 0.5D0 * rho * dx * dy * (w(i-1, j, k+1) + w(i, j, k+1))
                        Fb = 0.5D0 * rho * dx * dy * (w(i-1, j, k) + w(i, j, k))
                        delta_F = Fe - Fw + Fn - Fs + Ft - Fb
                  
                        aw_s = Dw + max(Fw, 0.0D0)
                        ae_s = De + max(-Fe, 0.0D0)
                        as_s = Ds + max(Fs, 0.0D0)   
                        an_s = Dn + max(-Fn, 0.0D0)  
                        ab_s = Db + max(Fb, 0.0D0)
                        at_s = Dt + max(-Ft, 0.0D0)
             
                        ap_s = aw_s + ae_s + as_s + an_s + at_s + ab_s + delta_F
                        ap = ap_o + theta * ap_s
                        ap_prime = ap_o - (1.0D0 - theta) * ap_s
                        
                        du(i, j, k) = (dy * dz) / ap
                        
                        w_node = aw_s * (theta * uold(i-1, j, k) + (1.0D0 - theta) * un(i-1, j, k))
                        e_node = ae_s * (theta * uold(i+1, j, k) + (1.0D0 - theta) * un(i+1, j, k))
                        n_node = an_s * (theta * uold(i, j+1, k) + (1.0D0 - theta) * un(i, j+1, k))
                        s_node = as_s * (theta * uold(i, j-1, k) + (1.0D0 - theta) * un(i, j-1, k))
                        b_node = ab_s * (theta * uold(i, j, k-1) + (1.0D0 - theta) * un(i, j, k-1))
                        t_node = at_s * (theta * uold(i, j, k+1) + (1.0D0 - theta) * un(i, j, k+1))
                        p_past = ap_prime * un(i, j, k)
                    
                        ustar(i, j, k) = (w_node + e_node + n_node + s_node + b_node + t_node + p_past + &
                                        (pstar(i-1, j, k) - pstar(i, j, k)) * dy * dz) / ap
                    enddo
                enddo
            enddo
            !$omp end parallel do
         
            ! v-momentum predictor
            !$omp parallel do collapse(3) private(i, j, k, Fe, Fw, Fn, Fs, Ft, Fb, delta_F, aw_s, ae_s, as_s, an_s, at_s, ab_s &
            !$omp , ap_s, ap_prime, ap, w_node, e_node, s_node, n_node, t_node, b_node, p_past)
            do i = 1, nx
                do j = 2, ny
                    do k = 1, nz
                        Fe = 0.5D0 * rho * dy * dz * (u(i+1, j-1, k) + u(i+1, j, k))
                        Fw = 0.5D0 * rho * dy * dz * (u(i, j-1, k) + u(i, j, k))
                        Fn = 0.5D0 * rho * dx * dz * (v(i, j, k) + v(i, j+1, k))
                        Fs = 0.5D0 * rho * dx * dz * (v(i, j-1, k) + v(i, j, k))
                        Ft = 0.5D0 * rho * dx * dy * (w(i, j-1, k+1) + w(i, j, k+1))
                        Fb = 0.5D0 * rho * dx * dy * (w(i, j-1, k) + w(i, j, k))
                        delta_F = Fe - Fw + Fn - Fs + Ft - Fb
                  
                        aw_s = Dw + max(Fw, 0.0D0)
                        ae_s = De + max(-Fe, 0.0D0)
                        as_s = Ds + max(Fs, 0.0D0)
                        an_s = Dn + max(-Fn, 0.0D0)
                        ab_s = Db + max(Fb, 0.0D0)
                        at_s = Dt + max(-Ft, 0.0D0)
                  
                        ap_s = aw_s + ae_s + as_s + an_s + ab_s + at_s + delta_F
                        ap = ap_o + theta * ap_s
                        ap_prime = ap_o - (1.0D0 - theta) * ap_s
                        
                        dv(i, j, k) = (dx * dz) / ap
                        
                        w_node = aw_s * (theta * vold(i-1, j, k) + (1.0D0 - theta) * vn(i-1, j, k))
                        e_node = ae_s * (theta * vold(i+1, j, k) + (1.0D0 - theta) * vn(i+1, j, k))
                        n_node = an_s * (theta * vold(i, j+1, k) + (1.0D0 - theta) * vn(i, j+1, k))
                        s_node = as_s * (theta * vold(i, j-1, k) + (1.0D0 - theta) * vn(i, j-1, k))
                        t_node = at_s * (theta * vold(i, j, k+1) + (1.0D0 - theta) * vn(i, j, k+1))
                        b_node = ab_s * (theta * vold(i, j, k-1) + (1.0D0 - theta) * vn(i, j, k-1))
                        p_past = ap_prime * vn(i, j, k)
                    
                        vstar(i, j, k) = (w_node + e_node + n_node + s_node + t_node + b_node + p_past + &
                                        (pstar(i, j-1, k) - pstar(i, j, k)) * dx * dz) / ap
                    enddo
                enddo
            enddo
            !$omp end parallel do

            ! w-momentum predictor
            !$omp parallel do collapse(3) private(i, j, k, Fe, Fw, Fn, Fs, Ft, Fb, delta_F, aw_s, ae_s, as_s, an_s, at_s, ab_s &
            !$omp , ap_s, ap_prime, ap, w_node, e_node, s_node, n_node, t_node, b_node, p_past)
            do i = 1, nx
                do j = 1, ny
                    do k = 2, nz
                        Fe = 0.5D0 * rho * dy * dz * (u(i+1, j, k-1) + u(i+1, j, k))
                        Fw = 0.5D0 * rho * dy * dz * (u(i, j, k-1) + u(i, j, k))
                        Fn = 0.5D0 * rho * dx * dz * (v(i, j+1, k-1) + v(i, j+1, k))
                        Fs = 0.5D0 * rho * dx * dz * (v(i, j, k-1) + v(i, j, k))
                        Ft = 0.5D0 * rho * dx * dy * (w(i, j, k) + w(i, j, k+1))
                        Fb = 0.5D0 * rho * dx * dy * (w(i, j, k-1) + w(i, j, k))
                        delta_F = Fe - Fw + Fn - Fs + Ft - Fb
                  
                        aw_s = Dw + max(Fw, 0.0D0)
                        ae_s = De + max(-Fe, 0.0D0)
                        as_s = Ds + max(Fs, 0.0D0)
                        an_s = Dn + max(-Fn, 0.0D0)
                        ab_s = Db + max(Fb, 0.0D0)
                        at_s = Dt + max(-Ft, 0.0D0)
                  
                        ap_s = aw_s + ae_s + as_s + an_s + ab_s + at_s + delta_F
                        ap = ap_o + theta * ap_s
                        ap_prime = ap_o - (1.0D0 - theta) * ap_s
                        
                        dww(i, j, k) = (dx * dy) / ap
                        
                        w_node = aw_s * (theta * wold(i-1, j, k) + (1.0D0 - theta) * wn(i-1, j, k))
                        e_node = ae_s * (theta * wold(i+1, j, k) + (1.0D0 - theta) * wn(i+1, j, k))
                        n_node = an_s * (theta * wold(i, j+1, k) + (1.0D0 - theta) * wn(i, j+1, k))
                        s_node = as_s * (theta * wold(i, j-1, k) + (1.0D0 - theta) * wn(i, j-1, k))
                        t_node = at_s * (theta * wold(i, j, k+1) + (1.0D0 - theta) * wn(i, j, k+1))
                        b_node = ab_s * (theta * wold(i, j, k-1) + (1.0D0 - theta) * wn(i, j, k-1))
                        p_past = ap_prime * wn(i, j, k)
                    
                        wstar(i, j, k) = (w_node + e_node + n_node + s_node + t_node + b_node + p_past + &
                                        (pstar(i, j, k-1) - pstar(i, j, k)) * dx * dy) / ap
                    enddo
                enddo
            enddo
            !$omp end parallel do
         
            ! Reapplying boundary conditions for u*, v*, and w*
            !$omp parallel do collapse(2)
            do i = 0, nx+1
                do k = 0, nz+1
                    ustar(i, 0, k)    = 0.0D0    
                    ustar(i, ny+1, k) = 0.0D0   
                    vstar(i, 1, k)    = 0.0D0    
                    vstar(i, ny+1, k) = 0.0D0   
                    wstar(i, 0, k)    = 0.0D0
                    wstar(i, ny+1, k) = 0.0D0 
                enddo
            enddo
            !$omp end parallel do
        
            !$omp parallel do collapse(2)
            do j = 0, ny+1
                do k = 0, nz+1
                    ustar(1, j, k)    = 0.0D0    
                    ustar(nx+1, j, k) = 0.0D0    
                    vstar(0, j, k)    = 0.0D0
                    vstar(nx+1, j, k) = 0.0D0
                    wstar(0, j, k)    = 0.0D0
                    wstar(nx+1, j, k) = 0.0D0
                enddo
            enddo
            !$omp end parallel do
        
            !$omp parallel do collapse(2)
            do i = 0, nx+1
                do j = 0, ny+1
                    ustar(i, j, 0)    = 0.0D0    
                    ustar(i, j, nz+1) = ulid    
                    vstar(i, j, 0)    = 0.0D0
                    vstar(i, j, nz+1) = 0.0D0
                    wstar(i, j, 1)    = 0.0D0
                    wstar(i, j, nz+1) = 0.0D0
                enddo
            enddo
            !$omp end parallel do
        
            ! Calculate mass source 
            bmax = 0.0D0
            !$omp parallel do collapse(3) private(i, j, k) reduction(max:bmax)
            do i = 1, nx
                do j = 1, ny
                    do k = 1, nz
                        pprime(i, j, k) = 0.0D0
                        bprime(i, j, k) = (rho * dy * dz * ustar(i, j, k)) - (rho * dy * dz * ustar(i+1, j, k)) &
                                        + (rho * dx * dz * vstar(i, j, k)) - (rho * dx * dz * vstar(i, j+1, k)) &
                                        + (rho * dx * dy * wstar(i, j, k)) - (rho * dx * dy * wstar(i, j, k+1))
                    
                        if (abs(bprime(i, j, k)) > bmax) bmax = abs(bprime(i, j, k))
                    enddo
                enddo
            enddo
            !$omp end parallel do
         
            ! Pressure correction equation 
            do p_iter = 1, 200
                !$omp parallel do collapse(3) private(i, j, k, aw, ae, as, an, ab, at, ap_p)
                do i = 1, nx
                    do j = 1, ny
                        do k = 1, nz
                            aw = rho * dy * dz * du(i, j, k)
                            ae = rho * dy * dz * du(i+1, j, k)
                            as = rho * dx * dz * dv(i, j, k)
                            an = rho * dx * dz * dv(i, j+1, k)
                            ab = rho * dx * dy * dww(i, j, k)
                            at = rho * dx * dy * dww(i, j, k+1)
                      
                            if (i == 1) aw = 0.0D0
                            if (i == nx) ae = 0.0D0
                            if (j == 1) as = 0.0D0
                            if (j == ny) an = 0.0D0
                            if (k == 1) ab = 0.0D0
                            if (k == nz) at = 0.0D0
        
                            ap_p = aw + ae + an + as + at + ab
                      
                            temp_pprime(i, j, k) = (ae * pprime(i+1, j, k) + aw * pprime(i-1, j, k) + &
                                                    an * pprime(i, j+1, k) + as * pprime(i, j-1, k) + &
                                                    at * pprime(i, j, k+1) + ab * pprime(i, j, k-1) + &
                                                    bprime(i, j, k)) / ap_p
                        enddo
                    enddo
                enddo
                !$omp end parallel do

                !$omp parallel do collapse(3) private(i, j, k)
                do i = 1, nx
                    do j = 1, ny
                        do k = 1, nz
                            pprime(i, j, k) = temp_pprime(i, j, k)
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            enddo
         
            ! Field updates
            !$omp parallel do collapse(3) private(i, j, k, uc, vc, wc)
            do i = 1, nx
                do j = 1, ny
                    do k = 1, nz
                        p(i, j, k) = pstar(i, j, k) + alph_p * pprime(i, j, k)
                        pstar(i, j, k) = p(i, j, k)
                 
                        if (i > 1) then
                            uc = ustar(i, j, k) + du(i, j, k) * (pprime(i-1, j, k) - pprime(i, j, k))
                            u(i, j, k) = alph_u * uc + (1.0D0 - alph_u) * uold(i, j, k)
                        endif
                   
                        if (j > 1) then
                            vc = vstar(i, j, k) + dv(i, j, k) * (pprime(i, j-1, k) - pprime(i, j, k))
                            v(i, j, k) = alph_v * vc + (1.0D0 - alph_v) * vold(i, j, k)
                        endif

                        if (k > 1) then
                            wc = wstar(i, j, k) + dww(i, j, k) * (pprime(i, j, k-1) - pprime(i, j, k))
                            w(i, j, k) = alph_w * wc + (1.0D0 - alph_w) * wold(i, j, k)
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
            
            ! Reapplying boundary conditions for u, v, w
            !$omp parallel do collapse(2)
            do i = 0, nx+1
                do k = 0, nz+1
                    u(i, 0, k)    = 0.0D0    
                    u(i, ny+1, k) = 0.0D0   
                    v(i, 1, k)    = 0.0D0    
                    v(i, ny+1, k) = 0.0D0   
                    w(i, 0, k)    = 0.0D0
                    w(i, ny+1, k) = 0.0D0 
                enddo
            enddo
            !$omp end parallel do
        
            !$omp parallel do collapse(2)
            do j = 0, ny+1
                do k = 0, nz+1
                    u(1, j, k)    = 0.0D0    
                    u(nx+1, j, k) = 0.0D0    
                    v(0, j, k)    = 0.0D0
                    v(nx+1, j, k) = 0.0D0
                    w(0, j, k)    = 0.0D0
                    w(nx+1, j, k) = 0.0D0
                enddo
            enddo
            !$omp end parallel do
        
            !$omp parallel do collapse(2)
            do i = 0, nx+1
                do j = 0, ny+1
                    u(i, j, 0)    = 0.0D0    
                    u(i, j, nz+1) = ulid    
                    v(i, j, 0)    = 0.0D0
                    v(i, j, nz+1) = 0.0D0
                    w(i, j, 1)    = 0.0D0
                    w(i, j, nz+1) = 0.0D0
                enddo
            enddo
            !$omp end parallel do
        
            if (bmax <= tol) exit
        enddo
        
        write(*,*) 'Time step: ', time_step, '| Inner Iters: ', iter, '| Max mass error: ', bmax
        
        ! Outputting solution
        if (mod(time_step, output_freq) == 0) then
            write(vtk_filename, '(A, I6.6, A)') 'transient_3D_', time_step, '.vtk'
            call write_vtk(trim(vtk_filename), u, v, w, p, nx, ny, nz, dx, dy, dz)
        endif

        
    enddo ! End of transient loop
    
! Stop the clock
    end_time = omp_get_wtime()   
    write(*,*) '------------------------------------'
    write(*,*) 'Total Time Steps: ', max_step
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'


! Deallocate memory
    deallocate(u, ustar, uold, du, un)
    deallocate(v, vstar, vold, dv, vn)
    deallocate(w, wstar, wold, dww, wn)
    deallocate(p, pstar, pprime, bprime)
    deallocate(temp_pprime)

    stop
end program Transient_Cavity_SIMPLE
    
    
! This subroutine aims at writing solution to Legacy VTK files.
subroutine write_vtk(filename, u, v, w, p, nx, ny, nz, dx, dy, dz)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny, nz
    real*8, intent(in) :: dx, dy, dz
    real*8, intent(in) :: u(0:nx+1, 0:ny+1, 0:nz+1)
    real*8, intent(in) :: v(0:nx+1, 0:ny+1, 0:nz+1)
    real*8, intent(in) :: w(0:nx+1, 0:ny+1, 0:nz+1)
    real*8, intent(in) :: p(0:nx+1, 0:ny+1, 0:nz+1)
    integer*4 :: i, j, k
    
    ! Open file
    open(unit=20, file=filename, status='replace')
    
    ! Header 
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') '3D CFD Transient Lid-Driven Cavity'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET STRUCTURED_POINTS'
    
    ! Grid properties
    write(20, '(A, I0, 1X, I0, 1X, I0)') 'DIMENSIONS ', nx, ny, nz
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ', 0.0D0, 0.0D0, 0.0D0
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'SPACING ', dx, dy, dz
    
    ! Data assignment
    write(20, '(A, I0)') 'POINT_DATA ', nx * ny * nz
    
    ! Pressure field 
    write(20, '(A)') 'SCALARS Pressure float 1'
    write(20, '(A)') 'LOOKUP_TABLE default'
    write(20, '(5E15.6)') (((p(i, j, k), i=1,nx), j=1,ny), k=1,nz)

    ! Velocity field 
    write(20, '(A)') 'VECTORS Velocity float'
    write(20, '(3E15.6)') ((( &
        0.5D0 * (u(i, j, k) + u(i+1, j, k)), &
        0.5D0 * (v(i, j, k) + v(i, j+1, k)), &
        0.5D0 * (w(i, j, k) + w(i, j, k+1)), &
        i=1,nx), j=1,ny), k=1,nz)
    
    close(20)
end subroutine write_vtk