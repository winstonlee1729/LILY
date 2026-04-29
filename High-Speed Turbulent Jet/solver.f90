program LES_Turbulent_Jet
    use omp_lib
    implicit none

! Declaration of variables
    integer*4 nx, ny
    integer*4 i, j, p_iter, time_step
    integer*4 max_step, output_freq
    
    ! Long domain to watch the jet break down into vortices
    parameter(nx = 250, ny = 100)
    parameter(max_step = 8000, output_freq = 50)

    ! Convergence and physical parameters
    real*8, parameter :: tol = 1.0D-5, sor_tol = 1.0D-6
    real*8, parameter :: rho = 1.0D0, mu = 0.0001D0
    real*8, parameter :: delta_t = 0.002D0
    real*8, parameter :: omega = 1.7D0, theta = 1.0D0
    
    ! Smagorinsky LES Model Parameter
    real*8, parameter :: C_s = 0.1D0

    real*8 :: Lx, Ly, dx, dy, volume, ap_o, filter_delta
    real*8 :: Fe, Fw, Fn, Fs
    real*8 :: aw, ae, as, an, ap
    real*8 :: aw_s, ae_s, an_s, as_s, ap_s
    real*8 :: w_node, e_node, n_node, s_node
    real*8 :: p_err, p_new, bmax1, bmax2
    real*8 :: DC_e, DC_w, DC_n, DC_s, S_dc
    
    ! Inlet conditions and LES variables
    real*8 :: U_jet, y_c, rnd_val, strain_rate, mu_eff_face
    real*8 :: start_time, end_time
    character(len=60) :: vtk_filename

    ! Flow field matrices
    real*8, allocatable :: u(:,:), ustar(:,:), ustar2(:,:), uold(:,:), du(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vstar2(:,:), vold(:,:), dv(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:), bprime(:,:)
    real*8, allocatable :: pprime2(:,:), bprime2(:,:)
    real*8, allocatable :: mu_t(:,:) 

! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), ustar2(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1), du(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vstar2(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1))
    allocate(pprime2(0:nx+1, 0:ny+1), bprime2(0:nx+1, 0:ny+1))
    allocate(mu_t(0:nx+1, 0:ny+1))

! Defining variables
    Lx = 5.0D0; Ly = 2.0D0
    dx = Lx / dble(nx); dy = Ly / dble(ny)
    volume = dx * dy
    ap_o = rho * volume / delta_t
    
    ! LES Filter width (grid size)
    filter_delta = sqrt(dx * dy)

    U_jet = 2.0D0

! Initialization
    u = 0.0D0
    v = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0
    pstar = 0.0D0
    uold = 0.0D0
    vold = 0.0D0
    du = 0.0D0
    dv = 0.0D0
    pprime = 0.0D0
    bprime = 0.0D0
    pprime2 = 0.0D0
    bprime2 = 0.0D0
    mu_t = 0.0D0
    
    call random_seed()

! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step
        
        ! Save old fields
        !$omp parallel do collapse(2)
        do j = 0, ny+1
            do i = 0, nx+1
                uold(i, j) = u(i, j)
                vold(i, j) = v(i, j)
            enddo
        enddo
        !$omp end parallel do

        ! SMAGORINSKY LES MODEL: Calculate Turbulent Eddy Viscosity 
        !$omp parallel do collapse(2) private(i, j, strain_rate)
        do j = 1, ny
            do i = 1, nx
                ! Magnitude of the strain rate tensor 
                strain_rate = sqrt(2.0D0 * ((uold(i+1, j) - uold(i, j)) / dx)**2 + &
                                   2.0D0 * ((vold(i, j+1) - vold(i, j)) / dy)**2 + &
                                   ((uold(i, j+1) - uold(i, j-1)) / (2.0D0*dy) + &
                                    (vold(i+1, j) - vold(i-1, j)) / (2.0D0*dx))**2)
                                    
                mu_t(i, j) = rho * (C_s * filter_delta)**2 * strain_rate
            enddo
        enddo
        !$omp end parallel do
        
        mu_t(0, :) = mu_t(1, :)
        mu_t(nx+1, :) = mu_t(nx, :)
        mu_t(:, 0) = mu_t(:, 1)
        mu_t(:, ny+1) = mu_t(:, ny)

        ! u-momentum predictor
        !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, aw_s, ae_s, as_s, an_s, ap_s, ap, &
        !$omp w_node, e_node, s_node, n_node, DC_e, DC_w, DC_n, DC_s, S_dc, mu_eff_face)
        do j = 1, ny
            do i = 2, nx
                Fe = 0.5D0 * rho * dy * (uold(i, j) + uold(i+1, j))
                Fw = 0.5D0 * rho * dy * (uold(i-1, j) + uold(i, j))
                Fn = 0.5D0 * rho * dx * (vold(i-1, j+1) + vold(i, j+1))
                Fs = 0.5D0 * rho * dx * (vold(i-1, j) + vold(i, j))
                
                DC_e = 0.0D0; DC_w = 0.0D0; DC_n = 0.0D0; DC_s = 0.0D0
                
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
                
                mu_eff_face = mu + 0.5D0 * (mu_t(i-1, j) + mu_t(i, j))
                
                aw_s = mu_eff_face * dy/dx + max( Fw, 0.0D0)
                ae_s = mu_eff_face * dy/dx + max(-Fe, 0.0D0)
                as_s = mu_eff_face * dx/dy + max( Fs, 0.0D0)
                an_s = mu_eff_face * dx/dy + max(-Fn, 0.0D0)
                
                ap_s = aw_s + ae_s + as_s + an_s 
                ap = ap_o + ap_s
                du(i, j) = dy / ap
                
                w_node = aw_s * uold(i-1, j)
                e_node = ae_s * uold(i+1, j)
                n_node = an_s * uold(i, j+1)
                s_node = as_s * uold(i, j-1)
                
                ustar(i, j) = (w_node + e_node + n_node + s_node + S_dc + &
                               ap_o * uold(i, j) + (p(i-1, j) - p(i, j)) * dy) / ap
            enddo
        enddo
        !$omp end parallel do

        ! v-momentum Predictor (QUICK + LES)
        !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, aw_s, ae_s, as_s, an_s, ap_s, ap, &
        !$omp w_node, e_node, s_node, n_node, DC_e, DC_w, DC_n, DC_s, S_dc, mu_eff_face)
        do j = 2, ny
            do i = 1, nx
                Fe = 0.5D0 * rho * dy * (uold(i+1, j-1) + uold(i+1, j))
                Fw = 0.5D0 * rho * dy * (uold(i, j-1)   + uold(i, j))
                Fn = 0.5D0 * rho * dx * (vold(i, j)     + vold(i, j+1))
                Fs = 0.5D0 * rho * dx * (vold(i, j-1)   + vold(i, j))
                
                DC_e = 0.0D0; DC_w = 0.0D0; DC_n = 0.0D0; DC_s = 0.0D0
                
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

                mu_eff_face = mu + 0.5D0 * (mu_t(i, j-1) + mu_t(i, j))

                aw_s = mu_eff_face * dy/dx + max( Fw, 0.0D0)
                ae_s = mu_eff_face * dy/dx + max(-Fe, 0.0D0)
                as_s = mu_eff_face * dx/dy + max( Fs, 0.0D0)
                an_s = mu_eff_face * dx/dy + max(-Fn, 0.0D0)
                
                ap_s = aw_s + ae_s + as_s + an_s 
                ap = ap_o + ap_s
                dv(i, j) = dx / ap
                
                w_node = aw_s * vold(i-1, j)
                e_node = ae_s * vold(i+1, j)
                n_node = an_s * vold(i, j+1) 
                s_node = as_s * vold(i, j-1)
                
                vstar(i, j) = (w_node + e_node + n_node + s_node + S_dc + &
                               ap_o * vold(i, j) + (p(i, j-1) - p(i, j)) * dx) / ap
            enddo
        enddo
        !$omp end parallel do

        ! Jet boundary conditions
        do j = 1, ny
            y_c = dble(j) * dy
            if (y_c > 0.4D0*Ly .and. y_c < 0.6D0*Ly) then
                call random_number(rnd_val)
                ustar(1, j) = U_jet + 0.05D0 * U_jet * (rnd_val - 0.5D0) 
            else
                ustar(1, j) = 0.0D0 
            endif
            
            vstar(0, j) = -vstar(1, j)
            ustar(nx+1, j) = ustar(nx, j) 
            vstar(nx+1, j) = vstar(nx, j)
        enddo
        do i = 1, nx
            ustar(i, 0) = -ustar(i, 1)
            ustar(i, ny+1) = -ustar(i, ny) ! Slip walls
            vstar(i, 1) = 0.0D0   
            vstar(i, ny+1) = 0.0D0
        enddo
        
        du(nx+1, :) = du(nx, :)

        ! First corrector
        bmax1 = 0.0D0
        !$omp parallel do collapse(2) reduction(max:bmax1)
        do j = 1, ny
            do i = 1, nx
                pprime(i, j) = 0.0D0
                bprime(i, j) = rho*dy*ustar(i,j) - rho*dy*ustar(i+1,j) + rho*dx*vstar(i,j) - rho*dx*vstar(i,j+1)
                if (abs(bprime(i, j)) > bmax1) bmax1 = abs(bprime(i, j))
            enddo
        enddo
        !$omp end parallel do

        do p_iter = 1, 500
            p_err = 0.0D0
            !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j+1, 2), nx, 2
                    
                    aw = rho*dy*du(i,j)
                    ae = rho*dy*du(i+1,j)
                    as = rho*dx*dv(i,j)
                    an = rho*dx*dv(i,j+1)
                    
                    if (i == 1)  aw = 0.0D0 
                    if (j == 1)  as = 0.0D0
                    if (j == ny) an = 0.0D0
                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae*pprime(min(i+1,nx),j) + aw*pprime(max(i-1,1),j) + an*pprime(i,j+1) + as*pprime(i,j-1) + bprime(i,j)) / ap
                        p_new = (1.0D0 - omega) * pprime(i,j) + omega * p_new
                        
                        if (abs(p_new - pprime(i,j)) > p_err) p_err = abs(p_new - pprime(i,j))
                        
                        pprime(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j, 2), nx, 2
                    
                    aw = rho*dy*du(i,j)
                    ae = rho*dy*du(i+1,j)
                    as = rho*dx*dv(i,j)
                    an = rho*dx*dv(i,j+1)
                    
                    if (i == 1)  aw = 0.0D0 
                    if (j == 1)  as = 0.0D0
                    if (j == ny) an = 0.0D0
                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae*pprime(min(i+1,nx),j) + aw*pprime(max(i-1,1),j) + an*pprime(i,j+1) + as*pprime(i,j-1) + bprime(i,j)) / ap
                        p_new = (1.0D0 - omega) * pprime(i,j) + omega * p_new
                        
                        if (abs(p_new - pprime(i,j)) > p_err) p_err = abs(p_new - pprime(i,j))
                        
                        pprime(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            pprime(nx+1, :) = 0.0D0 
            if (p_err < sor_tol) exit
        enddo

        ! Apply first corrector
        !$omp parallel do collapse(2)
        do j = 1, ny
            do i = 1, nx
                pstar(i, j) = p(i, j) + pprime(i, j)
                if (i > 1) ustar2(i, j) = ustar(i, j) + du(i, j)*(pprime(i-1, j) - pprime(i, j))
                if (j > 1) vstar2(i, j) = vstar(i, j) + dv(i, j)*(pprime(i, j-1) - pprime(i, j))
            enddo
        enddo
        !$omp end parallel do
        
        do j = 1, ny
            ustar2(1, j) = ustar(1, j)
            ustar2(nx+1, j) = ustar2(nx, j)
            vstar2(0, j) = -vstar2(1, j) 
            vstar2(nx+1, j) = vstar2(nx, j)
        enddo
        do i = 1, nx
            ustar2(i, 0) = -ustar2(i, 1)
            ustar2(i, ny+1) = -ustar2(i, ny)
            vstar2(i, 1) = 0.0D0  
            vstar2(i, ny+1) = 0.0D0
        enddo

        ! Second corrector step
        bmax2 = 0.0D0
        !$omp parallel do collapse(2) reduction(max:bmax2)
        do j = 1, ny
            do i = 1, nx
                pprime2(i, j) = 0.0D0
                bprime2(i, j) = rho*dy*ustar2(i,j) - rho*dy*ustar2(i+1,j) + rho*dx*vstar2(i,j) - rho*dx*vstar2(i,j+1)
                if (abs(bprime2(i, j)) > bmax2) bmax2 = abs(bprime2(i, j))
            enddo
        enddo
        !$omp end parallel do

        do p_iter = 1, 500
            p_err = 0.0D0
            !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j+1, 2), nx, 2
                    
                    aw = rho*dy*du(i,j)
                    ae = rho*dy*du(i+1,j)
                    as = rho*dx*dv(i,j)
                    an = rho*dx*dv(i,j+1)
                    
                    if (i == 1)  aw = 0.0D0
                    if (j == 1)  as = 0.0D0 
                    if (j == ny) an = 0.0D0
                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae*pprime2(min(i+1,nx),j) + aw*pprime2(max(i-1,1),j) + an*pprime2(i,j+1) + as*pprime2(i,j-1) + bprime2(i,j)) / ap
                        p_new = (1.0D0 - omega) * pprime2(i,j) + omega * p_new
                        
                        if (abs(p_new - pprime2(i,j)) > p_err) p_err = abs(p_new - pprime2(i,j))
                        
                        pprime2(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j, 2), nx, 2
                    aw = rho*dy*du(i,j)
                    ae = rho*dy*du(i+1,j)
                    as = rho*dx*dv(i,j)
                    an = rho*dx*dv(i,j+1)
                    
                    if (i == 1)  aw = 0.0D0
                    if (j == 1)  as = 0.0D0
                    if (j == ny) an = 0.0D0
                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae*pprime2(min(i+1,nx),j) + aw*pprime2(max(i-1,1),j) + an*pprime2(i,j+1) + as*pprime2(i,j-1) + bprime2(i,j)) / ap
                        p_new = (1.0D0 - omega) * pprime2(i,j) + omega * p_new
                        
                        if (abs(p_new - pprime2(i,j)) > p_err) p_err = abs(p_new - pprime2(i,j))
                        
                        pprime2(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do
            
            pprime2(nx+1, :) = 0.0D0 
            if (p_err < sor_tol) exit
        enddo

        ! Final Velocity Update
        !$omp parallel do collapse(2)
        do j = 1, ny
            do i = 1, nx
                p(i, j) = pstar(i, j) + pprime2(i, j)
                if (i > 1) u(i, j) = ustar2(i, j) + du(i, j)*(pprime2(i-1, j) - pprime2(i, j))
                if (j > 1) v(i, j) = vstar2(i, j) + dv(i, j)*(pprime2(i, j-1) - pprime2(i, j))
            enddo
        enddo
        !$omp end parallel do

        do j = 1, ny
            u(1, j) = ustar(1, j)
            u(nx+1, j) = u(nx, j)
            v(0, j) = -v(1, j)
            v(nx+1, j) = v(nx, j)
        enddo
        do i = 1, nx
            u(i, 0) = -u(i, 1)
            u(i, ny+1) = -u(i, ny)
            v(i, 1) = 0.0D0  
            v(i, ny+1) = 0.0D0
        enddo

        if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
            write(*,'(A,I6,A,ES12.4,A,ES12.4)') ' Step:', time_step, '  Err1:', bmax1, '  Err2:', bmax2
        endif

        if (mod(time_step, output_freq) == 0) then
            write(vtk_filename, '(A, I6.6, A)') 'LES_JET_', time_step, '.vtk'
            call write_vtk_les(trim(vtk_filename), u, v, p, mu_t, nx, ny, dx, dy)
        endif

    enddo

    end_time = omp_get_wtime()
    write(*,*) 'Wall time (s) : ', end_time - start_time

    deallocate(u, ustar, ustar2, uold, du)
    deallocate(v, vstar, vstar2, vold, dv)
    deallocate(p, pstar, pprime, pprime2, bprime, bprime2)
    deallocate(mu_t)
    stop
end program LES_Turbulent_Jet

subroutine write_vtk_les(filename, u, v, p, mu_t, nx, ny, dx, dy)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny
    real*8,    intent(in) :: dx, dy
    real*8,    intent(in) :: u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), p(0:nx+1, 0:ny+1), mu_t(0:nx+1, 0:ny+1)
    integer*4 :: i, j
    open(unit=21, file=filename, status='unknown')
    write(21, '(A)') '# vtk DataFile Version 3.0'
    write(21, '(A)') '2D LES Turbulent Jet'
    write(21, '(A)') 'ASCII'
    write(21, '(A)') 'DATASET STRUCTURED_POINTS'
    write(21, '(A, I0, 1X, I0, 1X, I0)')          'DIMENSIONS ', nx, ny, 1
    write(21, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ',     0.0D0, 0.0D0, 0.0D0
    write(21, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'SPACING ',    dx, dy, 1.0D0
    write(21, '(A, I0)')                          'POINT_DATA ', nx * ny
    write(21, '(A)') 'SCALARS Pressure float 1'
    write(21, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(21, '(5E15.6)') (p(i, j), i = 1, nx)
    enddo
    write(21, '(A)') 'SCALARS Eddy_Viscosity float 1'
    write(21, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(21, '(5E15.6)') (mu_t(i, j), i = 1, nx)
    enddo
    write(21, '(A)') 'VECTORS Velocity float'
    do j = 1, ny
        do i = 1, nx
            write(21, '(3E15.6)') 0.5D0 * (u(i, j) + u(i+1, j)), 0.5D0 * (v(i, j) + v(i, j+1)), 0.0D0
        enddo
    enddo
    close(21)
end subroutine write_vtk_les