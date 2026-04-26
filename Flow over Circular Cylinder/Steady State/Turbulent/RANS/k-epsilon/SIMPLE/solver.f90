program SS_Cylinder_RANS_k_eps
    use omp_lib
    implicit none
      
! Declaration of variables
    integer*4 nx, ny
    integer*4 iter, maxiter
    integer*4 i, j, k
    parameter(nx = 200, ny = 100)
    parameter(maxiter = 500000)
    
    ! Tolerance for convergence
    real*8, parameter :: tol = 1.0D-6
    
    ! Program time
    real*8 start_time, end_time
    
    ! Flow parameters 
    real*8 rho, mu, U_inf
    real*8 Lx, Ly, dx, dy, volume
    real*8 cx, cy, R, rad_sq
    parameter(rho = 1.0D0, mu = 1.0D-4)
    parameter(U_inf = 1.0D0)
    parameter(Lx = 20.0D0, Ly = 10.0D0)
    parameter(cx = 5.0D0, cy = 5.0D0, R = 0.5D0)
    
    ! Under-relaxation 
    real*8 alph_p, alph_u, alph_v, alph_k, alph_eps
    parameter(alph_p = 0.3D0, alph_u = 0.3D0, alph_v = 0.3D0)
    parameter(alph_k = 0.3D0, alph_eps = 0.3D0)
    
    ! Standard k-epsilon constants
    real*8 C_mu, C1_eps, C2_eps, sigma_k, sigma_eps
    parameter(C_mu = 0.09D0, C1_eps = 1.44D0, C2_eps = 1.92D0)
    parameter(sigma_k = 1.0D0, sigma_eps = 1.3D0)
    
    ! Turbulence boundary variables
    real*8 I_turb, L_turb, k_in, eps_in
    
    ! Flow fields
    real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), du(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vold(:,:), dv(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:), bprime(:,:), temp_pprime(:,:)
    
    ! Turbulence fields
    real*8, allocatable :: k_tke(:,:), k_old(:,:)
    real*8, allocatable :: eps(:,:), eps_old(:,:)
    real*8, allocatable :: mut(:,:), Pk(:,:)
    
    ! Equation coefficients
    real*8 ae, aw, an, as, ap, delta_F
    real*8 Fe, Fw, Fn, Fs
    real*8 De, Dw, Dn, Ds, bmax, uc, vc
    real*8 mu_e, mu_w, mu_n, mu_s
    real*8 dudx, dudy, dvdx, dvdy
    real*8 xc, yc, pc
    
    ! For mass conservation
    real*8 M_in, M_out, mass_ratio

! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1), du(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1))
    allocate(temp_pprime(0:nx+1, 0:ny+1))
    
    allocate(k_tke(0:nx+1, 0:ny+1), k_old(0:nx+1, 0:ny+1))
    allocate(eps(0:nx+1, 0:ny+1), eps_old(0:nx+1, 0:ny+1))
    allocate(mut(0:nx+1, 0:ny+1), Pk(0:nx+1, 0:ny+1))

! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    volume = dx * dy
    
    ! Calculate Inlet Turbulence (5% Intensity)
    I_turb = 0.05D0
    L_turb = 0.07D0 * (2.0D0 * R)
    k_in = 1.5D0 * (U_inf * I_turb)**2
    eps_in = (C_mu**0.75D0 * k_in**1.5D0) / L_turb

! Initialization
    u = 0.0D0
    v = 0.0D0
    p = 0.0D0
    pstar = 0.0D0
    
    ! Initialize fluid domain with freestream turbulence
    k_tke = k_in
    eps = eps_in
    mut = rho * C_mu * (k_in**2) / eps_in
    ! Applying inlet boundary condition
    !$omp parallel do
    do j = 0, ny+1
        u(1, j) = U_inf
        k_tke(1, j) = k_in
        eps(1, j) = eps_in
    enddo
    !$omp end parallel do
      
! Start the clock
    start_time = omp_get_wtime()
      
! RANS SIMPLE algorithm starts here
    do iter = 1, maxiter
        
        ! Save old fields
        !$omp parallel do collapse(2)
        do i = 0, nx+1
            do j = 0, ny+1
                uold(i, j) = u(i, j)
                vold(i, j) = v(i, j)
                k_old(i, j) = k_tke(i, j)
                eps_old(i, j) = eps(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        ! Update eddy viscosity (mut) and production (Pk)
        !$omp parallel do private(j, dudx, dudy, dvdx, dvdy, rad_sq)
        do i = 1, nx
            do j = 1, ny
                rad_sq = ((dble(i) - 0.5D0) * dx - cx)**2 + ((dble(j) - 0.5D0) * dy - cy)**2
                
                if (rad_sq <= R**2) then
                    Pk(i, j) = 0.0D0
                    mut(i, j) = 0.0D0
                else
                    dudx = (uold(i+1,j) - uold(i-1,j)) / (2.0D0 * dx)
                    dudy = (uold(i,j+1) - uold(i,j-1)) / (2.0D0 * dy)
                    dvdx = (vold(i+1,j) - vold(i-1,j)) / (2.0D0 * dx)
                    dvdy = (vold(i,j+1) - vold(i,j-1)) / (2.0D0 * dy)
                    
                    Pk(i, j) = mut(i,j) * (2.0D0*(dudx**2 + dvdy**2) + (dudy + dvdx)**2)
                    mut(i, j) = rho * C_mu * (k_old(i,j)**2) / max(eps_old(i,j), 1.0D-12)
                endif
            enddo
        enddo
        !$omp end parallel do
        
        ! u-momentum predictor
        !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw, ae, as, an, ap, rad_sq, mu_e, mu_w, mu_n, mu_s, De, Dw, Dn, Ds)
        do i = 2, nx
            do j = 1, ny
                rad_sq = ((dble(i) - 1.0D0) * dx - cx)**2 + ((dble(j) - 0.5D0) * dy - cy) **2
                
                if (rad_sq <= R**2) then
                    ap = 1.0D30
                    du(i, j) = 0.0D0
                    ustar(i, j) = 0.0D0
                else
                    Fe = 0.5D0 * rho * dy * (uold(i, j) + uold(i+1, j))
                    Fw = 0.5D0 * rho * dy * (uold(i-1, j) + uold(i, j))
                    Fn = 0.5D0 * rho * dx * (vold(i-1, j+1) + vold(i, j+1))
                    Fs = 0.5D0 * rho * dx * (vold(i-1, j) + vold(i, j))
                    delta_F = Fe - Fw + Fn - Fs
                    
                    mu_e = mu + 0.5D0 * (mut(i, j) + mut(i+1, j))
                    mu_w = mu + 0.5D0 * (mut(i-1, j) + mut(i, j))
                    mu_n = mu + 0.5D0 * (mut(i, j) + mut(i, j+1))
                    mu_s = mu + 0.5D0 * (mut(i, j-1) + mut(i, j))
                    
                    De = mu_e * (dy / dx)
                    Dw = mu_w * (dy / dx)
                    Dn = mu_n * (dx / dy)
                    Ds = mu_s * (dx / dy)
                
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)   
                    an = Dn + max(-Fn, 0.0D0)  
                    
                    ap = aw + ae + as + an + max(delta_F, 0.0D0)
                    ap = ap / alph_u
                
                    du(i, j) = dy / ap 
                
                    ustar(i, j) = (aw * uold(i-1, j) + ae * uold(i+1, j) + as * uold(i, j-1) + an * uold(i, j+1) &
                                + (1.0D0 - alph_u) * ap * uold(i, j) &
                                + (pstar(i-1, j) - pstar(i, j)) * dy) / ap
                endif
            enddo
        enddo
        !$omp end parallel do
         
        ! v-momentum predictor
        !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw, ae, as, an, ap, rad_sq, mu_e, mu_w, mu_n, mu_s, De, Dw, Dn, Ds)
        do i = 1, nx
            do j = 2, ny
                rad_sq = ((dble(i) - 0.5D0) * dx - cx)**2 + ((dble(j) - 1.0D0) * dy - cy)**2
                
                if (rad_sq <= R**2) then
                    ap = 1.0D30
                    dv(i, j) = 0.0D0
                    vstar(i, j) = 0.0D0
                else
                    Fe = 0.5D0 * rho * dy * (uold(i+1, j-1) + uold(i+1, j))
                    Fw = 0.5D0 * rho * dy * (uold(i, j-1) + uold(i, j))
                    Fn = 0.5D0 * rho * dx * (vold(i, j) + vold(i, j+1))
                    Fs = 0.5D0 * rho * dx * (vold(i, j-1) + vold(i, j))
                    delta_F = Fe - Fw + Fn - Fs
                    
                    mu_e = mu + 0.5D0 * (mut(i+1, j) + mut(i, j))
                    mu_w = mu + 0.5D0 * (mut(i, j) + mut(i-1, j))
                    mu_n = mu + 0.5D0 * (mut(i, j+1) + mut(i, j))
                    mu_s = mu + 0.5D0 * (mut(i, j) + mut(i, j-1))
                    
                    De = mu_e * (dy / dx)
                    Dw = mu_w * (dy / dx)
                    Dn = mu_n * (dx / dy)
                    Ds = mu_s * (dx / dy)
                
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)
                    an = Dn + max(-Fn, 0.0D0)
                    
                    ap = aw + ae + as + an + max(delta_F, 0.0D0)
                    ap = ap / alph_v
                
                    dv(i, j) = dx / ap 
                
                    vstar(i, j) = (aw * vold(i-1, j) + ae * vold(i+1, j) + as * vold(i, j-1) + an * vold(i, j+1) &
                                + (1.0D0 - alph_v) * ap * vold(i, j) &
                                + (pstar(i, j-1) - pstar(i, j)) * dx) / ap
                endif
            enddo
        enddo
        !$omp end parallel do
         
        ! Reapply boundary conditions for u* and v*
        do j = 0, ny+1  
            ustar(1, j) = U_inf                  
            vstar(0, j) = 0.0D0                  
            vstar(nx+1, j) = vstar(nx, j)         
        enddo
        
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
        
        do i = 0, nx+1
            vstar(i, 1) = 0.0D0                  
            vstar(i, ny+1) = 0.0D0               
            ustar(i, 0) = -ustar(i, 1)             
            ustar(i, ny+1) = -ustar(i, ny)         
        enddo
        
        ! Pressure correction equation
        bmax = 0.0D0
        !$omp parallel do private(j) reduction(max:bmax)
        do i = 1, nx
            do j = 1, ny
                pprime(i, j) = 0.0D0
                bprime(i, j) = (rho * dy * ustar(i, j)) - (rho * dy * ustar(i+1, j)) &
                             + (rho * dx * vstar(i, j)) - (rho * dx * vstar(i, j+1))
                
                if (((dble(i) - 0.5) * dx - cx)**2 + ((dble(j) - 0.5) * dy - cy)**2 > R**2) then
                    if (abs(bprime(i, j)) > bmax) bmax = abs(bprime(i, j))
                endif
            enddo
        enddo
        !$omp end parallel do
        
        do k = 1, 100
            !$omp parallel do private(j, aw, ae, as, an, ap)
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
                    ap = aw + ae + an + as
                    
                    if (ap < 1.0D-15) then
                        temp_pprime(i, j) = 0.0D0
                    else
                        temp_pprime(i, j) = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + &
                                            an * pprime(i, j+1) + as * pprime(i, j-1) + &
                                            bprime(i, j)) / ap
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            temp_pprime(1, 1) = 0.0D0
            
            !$omp parallel do collapse(2)
            do i = 1, nx
                do j = 1, ny
                    pprime(i, j) = temp_pprime(i, j)
                enddo
            enddo
            !$omp end parallel do
        enddo
        
        ! Correct fields
        !$omp parallel do private(j)
        do i = 1, nx
            do j = 1, ny
                p(i, j) = pstar(i, j) + alph_p * pprime(i, j)
                pstar(i, j) = p(i, j)
                
                if (i > 1) then
                    u(i, j) = ustar(i, j) + du(i, j) * (pprime(i-1, j) - pprime(i, j))
                endif
                if (j > 1) then
                    v(i, j) = vstar(i, j) + dv(i, j) * (pprime(i, j-1) - pprime(i, j))
                endif
            enddo
        enddo
        !$omp end parallel do
        
        ! Reapply boundary conditions for final u and v
        do j = 0, ny+1
            u(1, j) = U_inf                     
            v(0, j) = 0.0D0                     
            v(nx+1, j) = v(nx, j)               
        enddo
        
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
        
        do i = 0, nx+1
            v(i, 1) = 0.0D0                 
            v(i, ny+1) = 0.0D0              
            u(i, 0) = u(i, 1)               
            u(i, ny+1) = u(i, ny)           
        enddo


        ! Solve turbulent kinetic energy (k)
        !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw, ae, as, an, ap, rad_sq, mu_e, mu_w, mu_n, mu_s, De, Dw, Dn, Ds)
        do i = 1, nx
            do j = 1, ny
                rad_sq = ((dble(i) - 0.5D0) * dx - cx)**2 + ((dble(j) - 0.5D0) * dy - cy)**2
                
                if (rad_sq <= R**2) then
                    k_tke(i, j) = 1.0D-12
                else
                    Fe = 0.5D0 * rho * dy * (u(i, j) + u(i+1, j))
                    Fw = 0.5D0 * rho * dy * (u(i-1, j) + u(i, j))
                    Fn = 0.5D0 * rho * dx * (v(i, j) + v(i, j+1))
                    Fs = 0.5D0 * rho * dx * (v(i, j-1) + v(i, j))
                    delta_F = Fe - Fw + Fn - Fs
                    
                    mu_e = mu + 0.5D0 * (mut(i+1, j) + mut(i, j)) / sigma_k
                    mu_w = mu + 0.5D0 * (mut(i, j) + mut(i-1, j)) / sigma_k
                    mu_n = mu + 0.5D0 * (mut(i, j+1) + mut(i, j)) / sigma_k
                    mu_s = mu + 0.5D0 * (mut(i, j) + mut(i, j-1)) / sigma_k
                    
                    De = mu_e * (dy / dx)
                    Dw = mu_w * (dy / dx)
                    Dn = mu_n * (dx / dy)
                    Ds = mu_s * (dx / dy)
                
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)
                    an = Dn + max(-Fn, 0.0D0)
                    
                    ap = aw + ae + as + an + max(delta_F, 0.0D0) + (rho * eps_old(i,j) / max(k_old(i,j), 1.0D-12)) * volume
                    ap = ap / alph_k
                
                    k_tke(i, j) = (aw * k_old(i-1, j) + ae * k_old(i+1, j) + as * k_old(i, j-1) + an * k_old(i, j+1) &
                                + Pk(i, j) * volume &
                                + (1.0D0 - alph_k) * ap * k_old(i, j)) / ap
                                
                    k_tke(i, j) = max(k_tke(i, j), 1.0D-12)
                endif
            enddo
        enddo
        !$omp end parallel do


        ! Solve turbulent dissipation (eps)
        !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw, ae, as, an, ap, rad_sq, mu_e, mu_w, mu_n, mu_s, De, Dw, Dn, Ds)
        do i = 1, nx
            do j = 1, ny
                rad_sq = ((dble(i) - 0.5D0) * dx - cx)**2 + ((dble(j) - 0.5D0) * dy - cy)**2
                
                if (rad_sq <= R**2) then
                    eps(i, j) = 1.0D-12
                else
                    Fe = 0.5D0 * rho * dy * (u(i, j) + u(i+1, j))
                    Fw = 0.5D0 * rho * dy * (u(i-1, j) + u(i, j))
                    Fn = 0.5D0 * rho * dx * (v(i, j) + v(i, j+1))
                    Fs = 0.5D0 * rho * dx * (v(i, j-1) + v(i, j))
                    delta_F = Fe - Fw + Fn - Fs
                    
                    mu_e = mu + 0.5D0 * (mut(i+1, j) + mut(i, j)) / sigma_eps
                    mu_w = mu + 0.5D0 * (mut(i, j) + mut(i-1, j)) / sigma_eps
                    mu_n = mu + 0.5D0 * (mut(i, j+1) + mut(i, j)) / sigma_eps
                    mu_s = mu + 0.5D0 * (mut(i, j) + mut(i, j-1)) / sigma_eps
                    
                    De = mu_e * (dy / dx)
                    Dw = mu_w * (dy / dx)
                    Dn = mu_n * (dx / dy)
                    Ds = mu_s * (dx / dy)
                
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)
                    an = Dn + max(-Fn, 0.0D0)
                    
                    ap = aw + ae + as + an + max(delta_F, 0.0D0) + (C2_eps * rho * eps_old(i,j) / max(k_old(i,j), 1.0D-12)) * volume
                    ap = ap / alph_eps
                
                    eps(i, j) = (aw * eps_old(i-1, j) + ae * eps_old(i+1, j) + as * eps_old(i, j-1) + an * eps_old(i, j+1) &
                                + C1_eps * (eps_old(i,j) / max(k_old(i,j), 1.0D-12)) * Pk(i, j) * volume &
                                + (1.0D0 - alph_eps) * ap * eps_old(i, j)) / ap
                                
                    eps(i, j) = max(eps(i, j), 1.0D-12) 
                endif
            enddo
        enddo
        !$omp end parallel do

        ! Apply boundary conditions for k and eps
        do j = 0, ny+1  
            k_tke(1, j) = k_in
            eps(1, j) = eps_in
            k_tke(0, j) = k_in
            eps(0, j) = eps_in
            
            k_tke(nx+1, j) = k_tke(nx, j)
            eps(nx+1, j) = eps(nx, j)
        enddo
        
        do i = 0, nx+1
            k_tke(i, 1) = k_tke(i, 2)
            eps(i, 1) = eps(i, 2)
            k_tke(i, ny+1) = k_tke(i, ny)
            eps(i, ny+1) = eps(i, ny)
        enddo

        if (mod(iter, 100) == 0) then
            write(*,*) 'Iteration: ', iter, '|Max mass error: ', bmax
        endif
         
        if (bmax <= tol) exit
         
    enddo  ! End of SIMPLE algo.
      
! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total iterations: ', iter
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'

! Outputting solution (Added mut so you can plot Eddy Viscosity)
    open(unit=10, file='SS_Cylinder_RANS.txt', status='replace')
    do j = 1, ny
        do i = 1, nx
            xc = (dble(i) - 0.5D0) * dx
            yc = (dble(j) - 0.5D0) * dy
            uc = 0.5D0 * (u(i, j) + u(i+1, j))
            vc = 0.5D0 * (v(i, j) + v(i, j+1)) 
            pc = p(i, j)
            
            write(10, '(6E16.6)') xc, yc, uc, vc, pc, mut(i, j)
        enddo
    enddo

! Close the output file
    close(10)

! Memory Deallocation
    deallocate(u, ustar, uold, du)
    deallocate(v, vstar, vold, dv)
    deallocate(p, pstar, pprime, bprime, temp_pprime)
    deallocate(k_tke, k_old, eps, eps_old, mut, Pk)

    stop
end program SS_Cylinder_RANS_k_eps