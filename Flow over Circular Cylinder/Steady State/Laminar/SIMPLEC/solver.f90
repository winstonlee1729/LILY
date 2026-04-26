program SS_Cylinder_SIMPLEC
    use omp_lib
    implicit none
      
! Declaration of variables
    ! Grid size and iteration number
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
    real*8 Lx, Ly, dx, dy
    real*8 cx, cy, R, rad_sq
    parameter(rho = 1.0D0, mu = 0.05D0)
    parameter(U_inf = 1.0D0)
    parameter(Lx = 20.0D0, Ly = 10.0D0)
    parameter(cx = 5.0D0, cy = 5.0D0, R = 0.5D0)
    
    ! Under-relaxation
    real*8 alph_p, alph_u, alph_v
    parameter(alph_p = 0.8D0, alph_u = 0.3D0, alph_v = 0.3D0)
    
    ! Flow field
    real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), du(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vold(:,:), dv(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:),bprime(:,:)
    real*8, allocatable :: temp_pprime(:,:)
    
    ! Equation coefficients
    real*8 ae, aw, an, as, ap, delta_F
    real*8 Fe, Fw, Fn, Fs
    real*8 De, Dw, Dn, Ds, bmax, uc, vc
    real*8 xc, yc, pc
    
    ! For mass conservation
    real*8 M_in, M_out, mass_ratio

! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1))
    allocate(uold(0:nx+1, 0:ny+1), du(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1))
    allocate(vold(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1))
    allocate(pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1))
    allocate(temp_pprime(0:nx+1, 0:ny+1))

! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    
    De = mu * (dy / dx)
    Dw = mu * (dy / dx)
    Dn = mu * (dx / dy)
    Ds = mu * (dx / dy)

! Initialization
    u = 0.0D0
    v = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0
    pstar = 0.0D0
    !$omp parallel do
    ! Applying boundary condition
    do j = 0, ny+1
        u(1, j) = U_inf
        uold(1, j) = U_inf
    enddo
    !$omp end parallel do
      
! Start the clock
    start_time = omp_get_wtime()
      
! SIMPLEC algorithm starts here
    do iter = 1, maxiter
        
        ! Save old velocities
        !$omp parallel do private(j)
        do i = 0, nx+1
            do j = 0, ny+1
                uold(i, j) = u(i, j)
                vold(i, j) = v(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        ! u-momentum predictor
        !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw, ae, as, an, ap, rad_sq)
        do i = 2, nx
            do j = 1, ny
                ! Define cylinder position
                rad_sq = ((dble(i) - 1.0D0) * dx - cx)**2 + ((dble(j) - 0.5D0) * dy - cy) **2
                
                if (rad_sq <= R**2) then
                    ap = 1.0D30
                    du(i, j) = 0.0D0
                    ustar(i, j) = 0.0D0
                else
                    Fe = 0.5D0 * rho * dy * (u(i, j) + u(i+1, j))
                    Fw = 0.5D0 * rho * dy * (u(i-1, j) + u(i, j))
                    Fn = 0.5D0 * rho * dx * (v(i-1, j+1) + v(i, j+1))
                    Fs = 0.5D0 * rho * dx * (v(i-1, j) + v(i, j))
                    delta_F = Fe - Fw + Fn - Fs
                
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)   
                    an = Dn + max(-Fn, 0.0D0)  
                    ap = aw + ae + as + an + delta_F
                
                    du(i, j) = dy / (ap / alph_u - (aw + ae + an + as)) 
                
                    ustar(i, j) = (aw * uold(i-1, j) + ae * uold(i+1, j) + as * uold(i, j-1) + an * uold(i, j+1) &
                                + (pstar(i-1, j) - pstar(i, j)) * dy) / ap
                endif
            enddo
        enddo
        !$omp end parallel do
         
        ! v-momentum predictor
        !$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, aw, ae, as, an, ap, rad_sq)
        do i = 1, nx
            do j = 2, ny
                ! Define cylinder location
                rad_sq = ((dble(i) - 0.5D0) * dx - cx)**2 + ((dble(j) - 1.0D0) * dy - cy)**2
                
                if (rad_sq <= R**2) then
                    ap = 1.0D30
                    dv(i, j) = 0.0D0
                    vstar(i, j) = 0.0D0
                else
                    Fe = 0.5D0 * rho * dy * (u(i+1, j-1) + u(i+1, j))
                    Fw = 0.5D0 * rho * dy * (u(i, j-1) + u(i, j))
                    Fn = 0.5D0 * rho * dx * (v(i, j) + v(i, j+1))
                    Fs = 0.5D0 * rho * dx * (v(i, j-1) + v(i, j))
                    delta_F = Fe - Fw + Fn - Fs
                
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)
                    an = Dn + max(-Fn, 0.0D0)
                    ap = aw + ae + as + an + delta_F
                
                    dv(i, j) = dx / (ap / alph_v - (aw + ae + an + as))
                
                    vstar(i, j) = (aw * vold(i-1, j) + ae * vold(i+1, j) + as * vold(i, j-1) + an * vold(i, j+1) &
                                + (pstar(i, j-1) - pstar(i, j)) * dx) / ap
                endif
            enddo
        enddo
        !$omp end parallel do
         
        ! Reapplying boundary conditions for u* and v*
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
            ! If assuming flow is fully developed
            ! ustar(i, 0) = ustar(i, 1)           ! Bottom slip (zero shear)
            ! ustar(i, ny+1) = ustar(i, ny)       ! Top slip (zero shear)
        enddo
        
        ! Calculate mass source 
        bmax = 0.0D0
        !$omp parallel do private(j) reduction(max:bmax)
        do i = 1, nx
            do j = 1, ny
                pprime(i, j) = 0.0D0
                bprime(i, j) = (rho * dy * ustar(i, j)) - (rho * dy * ustar(i+1, j)) &
                             + (rho * dx * vstar(i, j)) - (rho * dx * vstar(i, j+1))
                
                ! Exclude solid cylinder from max error check
                if (((dble(i) - 0.5) * dx - cx)**2 + ((dble(j) - 0.5) * dy - cy)**2 > R**2) then
                    if (abs(bprime(i, j)) > bmax) bmax = abs(bprime(i, j))
                endif
            enddo
        enddo
        !$omp end parallel do
        
        ! Pressure correction equation
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
                        temp_pprime(i, j) = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + an * pprime(i, j+1) &
                                          + as * pprime(i, j-1) + bprime(i, j)) / ap
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
        

        ! Field updates
        !$omp parallel do private(j, uc, vc)
        do i = 1, nx
            do j = 1, ny
                p(i, j) = pstar(i, j) + alph_p * pprime(i, j)
                pstar(i, j) = p(i, j)
               
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
        
        !if (mod(iter, 100) == 0) then
            write(*,*) 'Iteration: ', iter, '|Max mass error: ', bmax
        !endif
            
            
        ! Reapplying boundary conditions for u and v
        ! Inlet and outlet
        do j = 0, ny+1
            u(1, j) = U_inf                 ! Inlet
            v(0, j) = 0.0D0                 ! v = 0 at inlet
            v(nx+1, j) = v(nx, j)           ! Zero gradient v at outlet
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
            u(i, 0) = -u(i, 1)              
            u(i, ny+1) = -u(i, ny)           
            ! If assuming flow is fully developed
            ! u(i, 0) = u(i, 1)               ! Bottom slip (zero shear)
            ! u(i, ny+1) = u(i, ny)           ! Top slip (zero shear)
        enddo
         
        if (bmax <= tol) exit
         
    enddo  ! SIMPLEC algorithm stops here
      
! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total iterations: ', iter
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'

! Outputting solution
    open(unit=10, file='SS_Cylinder_SIMPLEC.txt', status='unknown')
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

! Close the output file
    close(10)

! Memory Deallocation
    deallocate(u, ustar, uold, du)
    deallocate(v, vstar, vold, dv)
    deallocate(p, pstar, pprime, bprime)

    stop
end program SS_Cylinder_SIMPLEC