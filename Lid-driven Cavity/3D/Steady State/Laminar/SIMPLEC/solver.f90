program SIMPLEC_Lid_Driven_Cavity
    use omp_lib
    implicit none
    
! Declaration of variables  
    ! Grid size and iteration number 
    integer*4 nx, ny, nz
    integer*4 iter, maxiter, p_iter
    integer*4 i, j, k
    parameter(nx = 50, ny = 50, nz = 50)
    parameter(maxiter = 50000)
    
    ! Tolerance for convergence
    real*8, parameter :: tol = 1.0D-6
    
    ! Program time
    real*8 start_time, end_time
    
    ! Flow parameters
    real*8 rho, mu, ulid
    real*8 Lx, Ly, Lz, dx, dy, dz
    parameter(rho = 1.0D0, mu = 0.01D0)
    parameter(ulid = 5.0D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0)
    
    ! Under-relaxation
    real*8 alph_p, alph_u, alph_v, alph_w
    parameter(alph_p = 0.5D0, alph_u = 0.5D0, alph_v = 0.5D0, alph_w = 0.5D0)
    
    ! Flow field arrays
    real*8, allocatable :: u(:,:,:), ustar(:,:,:), uold(:,:,:), du(:,:,:)
    real*8, allocatable :: v(:,:,:), vstar(:,:,:), vold(:,:,:), dv(:,:,:)
    real*8, allocatable :: w(:,:,:), wstar(:,:,:), wold(:,:,:), dww(:,:,:)
    real*8, allocatable :: p(:,:,:), pstar(:,:,:), pprime(:,:,:), bprime(:,:,:)
    real*8, allocatable :: temp_pprime(:,:,:) 
    
    ! Output string for VTK files
    character(len=50) :: vtk_filename
    
    ! Equation coefficients
    real*8 aw, ae, an, as, at, ab, ap, delta_F
    real*8 Fe, Fw, Fn, Fs, Ft, Fb
    real*8 De, Dw, Dn, Ds, Dt, Db, bmax, uc, vc, wc
    
    ! Allocate memory
    allocate(u(0:nx+1, 0:ny+1, 0:nz+1), ustar(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(uold(0:nx+1, 0:ny+1, 0:nz+1), du(0:nx+1, 0:ny+1, 0:nz+1))
    
    allocate(v(0:nx+1, 0:ny+1, 0:nz+1), vstar(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(vold(0:nx+1, 0:ny+1, 0:nz+1), dv(0:nx+1, 0:ny+1, 0:nz+1))
    
    allocate(w(0:nx+1, 0:ny+1, 0:nz+1), wstar(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(wold(0:nx+1, 0:ny+1, 0:nz+1), dww(0:nx+1, 0:ny+1, 0:nz+1))
    
    allocate(p(0:nx+1, 0:ny+1, 0:nz+1), pstar(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(pprime(0:nx+1, 0:ny+1, 0:nz+1), bprime(0:nx+1, 0:ny+1, 0:nz+1))
    allocate(temp_pprime(0:nx+1, 0:ny+1, 0:nz+1))

! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    dz = Lz / dble(nz)
    
    De = mu * (dy * dz / dx)
    Dw = mu * (dy * dz / dx)
    Dn = mu * (dx * dz / dy)
    Ds = mu * (dx * dz / dy)
    Dt = mu * (dx * dy / dz)
    Db = mu * (dx * dy / dz)

! Initialization 
    u = 0.0D0
    v = 0.0D0
    w = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0
    wstar = 0.0D0
    pstar = 0.0D0 
    du = 0.0D0
    dv = 0.0D0
    dww = 0.0D0
    pprime = 0.0D0
    bprime = 0.0D0
    temp_pprime = 0.0D0
    ! Applying boundary condition
    do i = 0, nx+1
        do j = 0, ny+1
            u(i, j, nz+1) = ulid
        enddo
    enddo
      
! Start the clock
    start_time = omp_get_wtime()
      
! SIMPLEC algo. starts here
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
        !$omp parallel do collapse(3) private(i, j, k, Fe, Fw, Fn, Fs, Ft, Fb, delta_F, aw, ae, as, an, at, ab, ap)
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
                  
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)   
                    an = Dn + max(-Fn, 0.0D0)  
                    ab = Db + max(Fb, 0.0D0)
                    at = Dt + max(-Ft, 0.0D0)
             
                    ap = aw + ae + as + an + at + ab + delta_F
                    du(i, j, k) = (dy * dz) / (ap / alph_u - ap)
                    
                    ustar(i, j, k) = (aw * uold(i-1, j, k) + ae * uold(i+1, j, k) + &
                                      as * uold(i, j-1, k) + an * uold(i, j+1, k) + &
                                      ab * uold(i, j, k-1) + at * uold(i, j, k+1) + &
                                      (pstar(i-1, j, k) - pstar(i, j, k)) * dy * dz) / ap
                enddo
            enddo
        enddo
        !$omp end parallel do
         
        ! v-momentum predictor
        !$omp parallel do collapse(3) private(i, j, k, Fe, Fw, Fn, Fs, Ft, Fb, delta_F, aw, ae, as, an, at, ab, ap)
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
                  
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)
                    an = Dn + max(-Fn, 0.0D0)
                    ab = Db + max(Fb, 0.0D0)
                    at = Dt + max(-Ft, 0.0D0)
                  
                    ap = aw + ae + as + an + at + ab + delta_F
                    dv(i, j, k) = (dx * dz) / (ap / alph_v - ap)
                    
                    vstar(i, j, k) = (aw * vold(i-1, j, k) + ae * vold(i+1, j, k) + &
                                      as * vold(i, j-1, k) + an * vold(i, j+1, k) + &
                                      ab * vold(i, j, k-1) + at * vold(i, j, k+1) + &
                                      (pstar(i, j-1, k) - pstar(i, j, k)) * dx * dz) / ap
                enddo
            enddo
        enddo
        !$omp end parallel do

        ! w-momentum predictor
        !$omp parallel do collapse(3) private(i, j, k, Fe, Fw, Fn, Fs, Ft, Fb, delta_F, aw, ae, as, an, at, ab, ap)
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
                  
                    aw = Dw + max(Fw, 0.0D0)
                    ae = De + max(-Fe, 0.0D0)
                    as = Ds + max(Fs, 0.0D0)
                    an = Dn + max(-Fn, 0.0D0)
                    ab = Db + max(Fb, 0.0D0)
                    at = Dt + max(-Ft, 0.0D0)
                  
                    ap = aw + ae + as + an + at + ab + delta_F
                    dww(i, j, k) = (dx * dy) / (ap / alph_w - ap)
                    
                    wstar(i, j, k) = (aw * wold(i-1, j, k) + ae * wold(i+1, j, k) + &
                                      as * wold(i, j-1, k) + an * wold(i, j+1, k) + &
                                      ab * wold(i, j, k-1) + at * wold(i, j, k+1) + &
                                      (pstar(i, j, k-1) - pstar(i, j, k)) * dx * dy) / ap
                enddo
            enddo
        enddo
        !$omp end parallel do
         
        ! Reapplying boundary conditions
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
            !$omp parallel do collapse(3) private(i, j, k, aw, ae, as, an, ab, at, ap)
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
        
                        ap = aw + ae + an + as + at + ab
                      
                        temp_pprime(i, j, k) = (ae * pprime(i+1, j, k) + aw * pprime(i-1, j, k) + &
                                                an * pprime(i, j+1, k) + as * pprime(i, j-1, k) + &
                                                at * pprime(i, j, k+1) + ab * pprime(i, j, k-1) + &
                                                bprime(i, j, k)) / ap
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
        
        write(*,*) 'Iteration: ', iter, '|Max mass error: ', bmax
        
        ! Outputting solution
        if (mod(iter, 500) == 0) then
            write(vtk_filename, '(A, I6.6, A)') 'cavity_3D_', iter, '.vtk'
            call write_vtk(trim(vtk_filename), u, v, w, p, nx, ny, nz, dx, dy, dz)
        endif
         
        if (bmax <= tol) exit
         
    enddo  
      
! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total iterations: ', iter
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'

! Final output
    call write_vtk('cavity_SIMPLE_3D_final.vtk', u, v, w, p, nx, ny, nz, dx, dy, dz)

! Deallocate memory
    deallocate(u, ustar, uold, du)
    deallocate(v, vstar, vold, dv)
    deallocate(w, wstar, wold, dww)
    deallocate(p, pstar, pprime, bprime)
    deallocate(temp_pprime)

    stop
end program SIMPLEC_Lid_Driven_Cavity



! This subroutine aims at outputting velocity vectors and pressure as Legacy VTK format. 
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
    real*8 :: uc, vc, wc
    
    ! Open file
    open(unit=20, file=filename, status='replace')
    
    ! Header
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') '3D CFD Lid-Driven Cavity Output'
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
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                write(20, '(E14.6)') p(i, j, k)
            enddo
        enddo
    enddo

    ! Velocity field 
    write(20, '(A)') 'VECTORS Velocity float'
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                uc = 0.5D0 * (u(i, j, k) + u(i+1, j, k))
                vc = 0.5D0 * (v(i, j, k) + v(i, j+1, k))
                wc = 0.5D0 * (w(i, j, k) + w(i, j, k+1))
                
                write(20, '(E14.6, 1X, E14.6, 1X, E14.6)') uc, vc, wc
            enddo
        enddo
    enddo
    
    close(20)
end subroutine write_vtk