program Kelvin_Helmholtz
    use omp_lib
    implicit none

! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny
    integer*4 time_step, max_step, output_freq
    integer*4 i, j, k
    integer*4 p_iter, max_p_iter
    integer*4 im1, ip1
    integer*4 pert_mode
    parameter(nx = 400, ny = 200)
    parameter(max_step = 9000, output_freq = 100)
    parameter(max_p_iter = 5000)
    parameter(pert_mode = 4)
    
    ! For program convergence
    real*8, parameter :: p_tol = 1.0D-8
    
    ! Flow field parameters
    real*8 rho, mu, U0, dt
    real*8 Lx, Ly, dx, dy, pi
    real*8 ya, yb, shear_thick, pert_amp
    parameter(rho = 1.0D0, mu = 1.0D-4, U0 = 0.5D0)
    parameter(Lx = 2.0D0, Ly = 1.0D0)
    parameter(dt = 5.0D-4)
    parameter(shear_thick = 0.025D0, pert_amp = 1.0D-2)
    
    ! SOR factor
    reakl*8, parameter :: omega = 1.7D0
    
    ! Temporary variables and coefficients
    real*8 x_face, x_cell, y_cell, y_face
    real*8 Fe, Fw, Fn, Fs
    real*8 aw, ae, as, an, ap, ap_o
    real*8 p_new, p_err
    real*8 max_div, t_now, p_mean, b_mean
    real*8 tr_flux
    
    ! For program time
    real*8 start_time, end_time
    
    ! For outputting solution
    character(len=60) :: vtk_filename

    ! Flow field matrices
    real*8, allocatable :: u(:,:), ustar(:,:), ustarstar(:,:), uold(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vstarstar(:,:), vold(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pstarstar(:,:)
    real*8, allocatable :: pprime(:,:), pprime2(:,:)
    real*8, allocatable :: bprime(:,:), bprime2(:,:)
    real*8, allocatable :: uprime(:,:), vprime(:,:)
    real*8, allocatable :: ucorr(:,:), vcorr(:,:)
    real*8, allocatable :: du(:,:), dv(:,:)
    real*8, allocatable :: tr(:,:), tr_old(:,:)

! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), ustarstar(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vstarstar(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pstarstar(0:nx+1, 0:ny+1))
    allocate(pprime(0:nx+1, 0:ny+1), pprime2(0:nx+1, 0:ny+1))
    allocate(bprime(0:nx+1, 0:ny+1), bprime2(0:nx+1, 0:ny+1))
    allocate(uprime(0:nx+1, 0:ny+1), vprime(0:nx+1, 0:ny+1))
    allocate(ucorr(0:nx+1, 0:ny+1), vcorr(0:nx+1, 0:ny+1))
    allocate(du(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
    allocate(tr(0:nx+1, 0:ny+1), tr_old(0:nx+1, 0:ny+1))

! Defining variable
    pi = acos(-1.0D0)
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    ap_o = rho * dx * dy / dt
    
    ya = 0.25D0 * Ly         ! Lower interface
    yb = 0.75D0 * Ly         ! Upper interface

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
    uold = 0.0D0
    vold = 0.0D0
    uprime = 0.0D0
    vprime = 0.0D0
    ucorr = 0.0D0
    vcorr = 0.0D0
    pprime = 0.0D0
    pprime2 = 0.0D0
    bprime = 0.0D0
    bprime2 = 0.0D0
    du = 0.0D0
    dv = 0.0D0
    tr = 0.0D0
    tr_old = 0.0D0
    ! Applying ICs and BCs
    ! Initial conditions: (1) 3-band tanh shear (2) Sine perturbation (3) Tracer
    do j = 1, ny
        y_cell = (dble(j) - 0.5D0) * dy
        y_face = dble(j-1) * dy
        do i = 1, nx
            x_cell = (dble(i) - 0.5D0) * dx

            u(i, j) = U0 * (tanh((y_cell - ya) / shear_thick) - tanh((y_cell - yb) / shear_thick) - 1.0D0)

            v(i, j) = pert_amp * sin(2.0D0 * pi * dble(pert_mode) * x_cell / Lx) * &
                      (exp(-((y_face - ya) / (2.0D0 * shear_thick))**2) + &
                       exp(-((y_face - yb) / (2.0D0 * shear_thick))**2))

            tr(i, j) = 0.5D0 * (tanh((y_cell - ya) / shear_thick) - tanh((y_cell - yb) / shear_thick))
        enddo
    enddo
    
    !$omp parallel do
    do j = 0, ny+1
        u(0, j) = u(nx, j)
        u(nx+1, j) = u(1, j)
        v(0, j) = v(nx, j)
        v(nx+1, j) = v(1, j)
        tr(0, j) = tr(nx, j)
        tr(nx+1, j) = tr(1, j)
    enddo
    !$omp end parallel do
    
    !$omp parallel do
    do i = 0, nx+1
        u(i, 0) = u(i, 1)
        u(i, ny+1) = u(i, ny)
        v(i, 1) = 0.0D0   
        v(i, ny+1) = 0.0D0
        v(i, 0) = -v(i, 1)
        tr(i, 0) = tr(i, 1) 
        tr(i, ny+1) = tr(i, ny)
    enddo
    !$omp end parallel do

! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step
        t_now = dble(time_step) * dt

        ! Save fields for this time step
        !$omp parallel do collapse(2)
        do i = 0, nx+1
            do j = 0, ny+1
                tr_old(i, j) = tr(i, j)
                uold(i, j) = u(i, j)
                vold(i, j) = v(i, j)
                pstar(i, j) = p(i, j)
            enddo
        enddo
        !$omp end parallel do

        ! For tracer advection
        !$omp parallel do collapse(2) private(i, j, im1, ip1, tr_flux)
        do j = 1, ny
            do i = 1, nx
                
                im1 = i - 1
                if (im1 < 1) im1 = nx
                
                ip1 = i + 1
                if (ip1 > nx) ip1 = 1

                tr_flux = 0.0D0
                if (u(i, j) > 0.0D0) then
                    tr_flux = tr_flux + u(i, j) * tr_old(im1, j) * dy
                else
                    tr_flux = tr_flux + u(i, j) * tr_old(i, j) * dy
                endif
                
                if (u(ip1, j) > 0.0D0) then
                    tr_flux = tr_flux - u(ip1, j) * tr_old(i, j) * dy
                else
                    tr_flux = tr_flux - u(ip1, j) * tr_old(ip1, j) * dy
                endif
                
                if (v(i, j) > 0.0D0) then
                    tr_flux = tr_flux + v(i, j) * tr_old(i, j-1) * dx
                else
                    tr_flux = tr_flux + v(i, j) * tr_old(i, j) * dx
                endif
                
                if (v(i, j+1) > 0.0D0) then
                    tr_flux = tr_flux - v(i, j+1) * tr_old(i, j) * dx
                else
                    tr_flux = tr_flux - v(i, j+1) * tr_old(i, j+1) * dx
                endif

                tr(i, j) = tr_old(i, j) + (dt / (dx * dy)) * tr_flux
            enddo
        enddo
        !$omp end parallel do
        
        ! Apply tracer BCs
        !$omp parallel do
        do j = 0, ny+1
            tr(0, j) = tr(nx, j)
            tr(nx+1, j) = tr(1, j)
        enddo
        !$omp end parallel do
        
        !$omp parallel do
        do i = 0, nx+1
            tr(i, 0) = tr(i, 1)
            tr(i, ny+1) = tr(i, ny)
        enddo
        !$omp end parallel do

        ! Predictor step 
        !$omp parallel do collapse(2) private(i, j, im1, ip1, Fw, Fe, Fn, Fs, aw, ae, as, an, ap)
        do j = 1, ny
            do i = 1, nx
                
                im1 = i - 1
                if (im1 < 1) im1 = nx
                
                ip1 = i + 1
                if (ip1 > nx) ip1 = 1

                Fe = 0.5D0 * rho * dy * (uold(i, j) + uold(ip1, j))
                Fw = 0.5D0 * rho * dy * (uold(im1, j) + uold(i, j))
                Fn = 0.5D0 * rho * dx * (vold(im1, j+1) + vold(i, j+1))
                Fs = 0.5D0 * rho * dx * (vold(im1, j) + vold(i, j))

                aw = mu * dy / dx + max(Fw, 0.0D0)
                ae = mu * dy / dx + max(-Fe, 0.0D0)
                as = mu * dx / dy + max(Fs, 0.0D0)
                an = mu * dx / dy + max(-Fn, 0.0D0)
                ap = ap_o + aw + ae + as + an + (Fe - Fw + Fn - Fs)
                
                du(i, j) = dy / ap

                ustar(i, j) = (aw * uold(im1, j) + ae * uold(ip1, j) + as * uold(i, j-1) + an * uold(i, j+1) &
                             + ap_o * uold(i, j) + (pstar(im1, j) - pstar(i, j)) * dy) / ap
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(2) private(i, j, im1, ip1, Fw, Fe, Fn, Fs, aw, ae, as, an, ap)
        do j = 2, ny
            do i = 1, nx
                
                im1 = i - 1
                if (im1 < 1) im1 = nx
                
                ip1 = i + 1
                if (ip1 > nx) ip1 = 1

                Fe = 0.5D0 * rho * dy * (uold(ip1, j-1) + uold(ip1, j))
                Fw = 0.5D0 * rho * dy * (uold(i, j-1) + uold(i, j))
                Fn = 0.5D0 * rho * dx * (vold(i, j) + vold(i, j+1))
                Fs = 0.5D0 * rho * dx * (vold(i, j-1) + vold(i, j))

                aw = mu * dy / dx + max(Fw, 0.0D0)
                ae = mu * dy / dx + max(-Fe, 0.0D0)
                as = mu * dx / dy + max(Fs, 0.0D0)
                an = mu * dx / dy + max(-Fn, 0.0D0)
                ap = ap_o + aw + ae + as + an + (Fe - Fw + Fn - Fs)
                
                dv(i, j) = dx / ap

                vstar(i, j) = (aw * vold(im1, j) + ae * vold(ip1, j) + as * vold(i, j-1) + an * vold(i, j+1) &
                             + ap_o * vold(i, j) + (pstar(i, j-1) - pstar(i, j)) * dx) / ap
            enddo
        enddo
        !$omp end parallel do

        ! Apply BCs to predicted velocities
        !$omp parallel do
        do j = 0, ny+1
            ustar(0, j) = ustar(nx, j)
            ustar(nx+1, j) = ustar(1, j)
            vstar(0, j) = vstar(nx, j)
            vstar(nx+1, j) = vstar(1, j)
        enddo
        !$omp end parallel do
        
        !$omp parallel do
        do i = 0, nx+1
            ustar(i, 0) = ustar(i, 1)
            ustar(i, ny+1) = ustar(i, ny)
            vstar(i, 1) = 0.0D0
            vstar(i, ny+1) = 0.0D0
        enddo
        !$omp end parallel do
        
        ! First corrector 
        b_mean = 0.0D0
        !$omp parallel do collapse(2) private(ip1) reduction(+:b_mean)
        do j = 1, ny
            do i = 1, nx
                ip1 = i + 1
                if (ip1 > nx) ip1 = 1
                
                pprime(i, j) = 0.0D0
                
                bprime(i, j) = rho * dy * ustar(i, j) - rho * dy * ustar(ip1, j) + rho * dx * vstar(i, j) - rho * dx * vstar(i, j+1)
                b_mean = b_mean + bprime(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        b_mean = b_mean / dble(nx*ny)
        !$omp parallel do collapse(2)
        do j = 1, ny
            do i = 1, nx
                bprime(i, j) = bprime(i, j) - b_mean
            enddo
        enddo
        !$omp end parallel do

        do p_iter = 1, max_p_iter
            p_err = 0.0D0
            !$omp parallel do private(i, j, im1, ip1, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j+1, 2), nx, 2
                    
                    im1 = i - 1
                    if (im1 < 1) im1 = nx
                    
                    ip1 = i + 1 
                    if (ip1 > nx) ip1 = 1
                    
                    aw = rho * dy * du(i, j)
                    ae = rho * dy * du(ip1, j)
                    
                    if (j == 1) then
                        as = 0.0D0
                    else
                        as = rho * dx * dv(i, j)
                    endif
                    
                    if (j == ny) then
                        an = 0.0D0
                    else
                        an = rho * dx * dv(i, j+1)
                    endif
                    
                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae * pprime(ip1, j) + aw * pprime(im1, j) + an * pprime(i, j+1) + as * pprime(i, j-1) + bprime(i, j)) / ap
                        p_new = (1.0D0 - omega)*pprime(i, j) + omega*p_new
                        
                        if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                        
                        pprime(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do

            !$omp parallel do private(i, j, im1, ip1, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j, 2), nx, 2
                    
                    im1 = i - 1
                    if (im1 < 1) im1 = nx
                
                    ip1 = i + 1
                    if (ip1 > nx) ip1 = 1
                    
                    aw = rho * dy * du(i, j)
                    ae = rho * dy * du(ip1, j)
                    
                    if (j == 1) then
                        as = 0.0D0
                    else
                        as = rho * dx * dv(i, j)
                    endif
                    
                    if (j == ny) then
                        an = 0.0D0
                    else
                        an = rho * dx * dv(i, j+1)
                    endif

                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae * pprime(ip1, j) + aw * pprime(im1, j) + an * pprime(i, j+1) + as * pprime(i, j-1) + bprime(i, j)) / ap
                        p_new = (1.0D0 - omega) * pprime(i, j) + omega * p_new
                        
                        if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                        
                        pprime(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do

            if (p_err < p_tol) exit
        enddo

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
        do j = 0, ny+1
            do i = 0, nx+1
                pprime(i, j) = pprime(i, j) - p_mean
            enddo
        enddo
        !$omp end parallel do

        ! Applying first correction
        !$omp parallel do collapse(2) private(i, j, im1)
        do j = 1, ny
            do i = 1, nx
                
                im1 = i - 1
                if (im1 < 1) im1 = nx
                
                pstarstar(i, j) = pstar(i, j) + pprime(i, j)
                
                uprime(i, j) = du(i, j) * (pprime(im1, j) - pprime(i, j))
                ustarstar(i, j) = ustar(i, j) + uprime(i, j)
                
                if (j > 1) then
                    vprime(i, j) = dv(i, j) * (pprime(i, j-1) - pprime(i, j))
                    vstarstar(i, j) = vstar(i, j) + vprime(i, j)
                endif
            enddo
        enddo
        !$omp end parallel do

        ! Apply BCs to first corrected velocities
        do j = 0, ny+1
            ustarstar(0, j) = ustarstar(nx, j)
            ustarstar(nx+1, j) = ustarstar(1, j)
            vstarstar(0, j) = vstarstar(nx, j)
            vstarstar(nx+1, j) = vstarstar(1, j)
            uprime(0, j) = uprime(nx, j)
            uprime(nx+1, j) = uprime(1, j)
            vprime(0, j) = vprime(nx, j)
            vprime(nx+1, j) = vprime(1, j)
        enddo
        do i = 0, nx+1
            ustarstar(i, 0) = ustarstar(i, 1)
            ustarstar(i, ny+1) = ustarstar(i, ny)
            vstarstar(i, 1) = 0.0D0
            vstarstar(i, ny+1) = 0.0D0
            uprime(i, 0) = uprime(i, 1)
            uprime(i, ny+1) = uprime(i, ny)
            vprime(i, 1) = 0.0D0
            vprime(i, ny+1) = 0.0D0
        enddo

        ! Second corrector
        !$omp parallel do collapse(2) private(i, j, im1, ip1, Fw, Fe, Fn, Fs, aw, ae, as, an, ap)
        do j = 1, ny
            do i = 1, nx
                
                im1 = i - 1
                if (im1 < 1) im1 = nx
                
                ip1 = i + 1
                if (ip1 > nx) ip1 = 1

                Fe = 0.5D0 * rho * dy * (uold(i, j) + uold(ip1, j))
                Fw = 0.5D0 * rho * dy * (uold(im1, j) + uold(i, j))
                Fn = 0.5D0 * rho * dx * (vold(im1, j+1) + vold(i, j+1))
                Fs = 0.5D0 * rho * dx * (vold(im1, j) + vold(i, j))

                aw = mu * dy / dx + max(Fw, 0.0D0)
                ae = mu * dy / dx + max(-Fe, 0.0D0)
                as = mu * dx / dy + max(Fs, 0.0D0)
                an = mu * dx / dy + max(-Fn, 0.0D0)
                ap = ap_o + aw + ae + as + an + (Fe - Fw + Fn - Fs)
                
                ucorr(i, j) = (aw * uprime(im1, j) + ae * uprime(ip1, j) + as * uprime(i, j-1) + an * uprime(i, j+1)) / ap
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(2) private(i, j, im1, ip1, Fw, Fe, Fn, Fs, aw, ae, as, an, ap)
        do j = 2, ny
            do i = 1, nx
                
                im1 = i - 1
                if (im1 < 1) im1 = nx
                
                ip1 = i + 1
                if (ip1 > nx) ip1 = 1

                Fe = 0.5D0 * rho * dy * (uold(ip1, j-1) + uold(ip1, j))
                Fw = 0.5D0 * rho * dy * (uold(i, j-1) + uold(i, j))
                Fn = 0.5D0 * rho * dx * (vold(i, j) + vold(i, j+1))
                Fs = 0.5D0 * rho * dx * (vold(i, j-1) + vold(i, j))

                aw = mu * dy / dx + max(Fw, 0.0D0)
                ae = mu * dy / dx + max(-Fe, 0.0D0)
                as = mu * dx / dy + max(Fs, 0.0D0)
                an = mu * dx / dy + max(-Fn, 0.0D0)
                ap = ap_o + aw + ae + as + an + (Fe - Fw + Fn - Fs)
                
                vcorr(i, j) = (aw * vprime(im1, j) + ae * vprime(ip1, j) + as * vprime(i, j-1) + an * vprime(i, j+1)) / ap
            enddo
        enddo
        !$omp end parallel do

        ! Apply BCs to second corrected velocities
        !$omp parallel do
        do j = 0, ny+1
            ucorr(0, j) = ucorr(nx, j)
            ucorr(nx+1, j) = ucorr(1, j)
            vcorr(0, j) = vcorr(nx, j) 
            vcorr(nx+1, j) = vcorr(1, j)
        enddo
        !$omp end parallel do
        
        !$omp parallel do
        do i = 0, nx+1
            ucorr(i, 0) = ucorr(i, 1)
            ucorr(i, ny+1) = ucorr(i, ny)
            vcorr(i, 1) = 0.0D0
            vcorr(i, ny+1) = 0.0D0
        enddo
        !$omp end parallel do
        
        b_mean = 0.0D0
        !$omp parallel do collapse(2) private(ip1) reduction(+:b_mean)
        do j = 1, ny
            do i = 1, nx
                
                ip1 = i + 1 
                if (ip1 > nx) ip1 = 1
                
                pprime2(i, j) = 0.0D0
                bprime2(i, j) = rho * dy * ucorr(i, j) - rho * dy * ucorr(ip1, j) + rho * dx * vcorr(i, j) - rho * dx * vcorr(i, j+1)
                b_mean = b_mean + bprime2(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        b_mean = b_mean / dble(nx*ny)
        !$omp parallel do collapse(2)
        do j = 1, ny
            do i = 1, nx
                bprime2(i, j) = bprime2(i, j) - b_mean
            enddo
        enddo
        !$omp end parallel do

        do k = 1, max_p_iter
            p_err = 0.0D0
            
            ! Red cells
            !$omp parallel do private(i, j, im1, ip1, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j+1, 2), nx, 2
                    im1 = i - 1
                    if (im1 < 1) im1 = nx
                    
                    ip1 = i + 1
                    if (ip1 > nx) ip1 = 1
                    
                    aw = rho * dy * du(i, j)
                    ae = rho * dy * du(ip1, j)
                    
                    if (j == 1) then
                        as = 0.0D0
                    else
                        as = rho * dx * dv(i, j)
                    endif
                    
                    if (j == ny) then
                        an = 0.0D0
                    else
                        an = rho * dx * dv(i, j+1)
                    endif
                    
                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae * pprime2(ip1, j) + aw * pprime2(im1, j) + an * pprime2(i, j+1) + as * pprime2(i, j-1) + bprime2(i, j)) / ap
                        p_new = (1.0D0 - omega) * pprime2(i, j) + omega * p_new
                        
                        if (abs(p_new - pprime2(i, j)) > p_err) p_err = abs(p_new - pprime2(i, j))
                        
                        pprime2(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do

            ! Black cells
            !$omp parallel do private(i, j, im1, ip1, aw, ae, as, an, ap, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j, 2), nx, 2
                    
                    im1 = i - 1
                    if (im1 < 1) im1 = nx
                    
                    ip1 = i + 1
                    if (ip1 > nx) ip1 = 1
                    
                    aw = rho * dy * du(i, j)
                    ae = rho * dy * du(ip1, j)
                    
                    if (j == 1) then
                        as = 0.0D0
                    else
                        as = rho * dx * dv(i, j)
                    endif
                    
                    if (j == ny) then
                        an = 0.0D0
                    else
                        an = rho * dx * dv(i, j+1)
                    endif

                    ap = aw + ae + as + an
                    
                    if (ap > 1.0D-15) then
                        p_new = (ae * pprime2(ip1, j) + aw * pprime2(im1, j) + an * pprime2(i, j+1) + as * pprime2(i, j-1) + bprime2(i, j)) / ap
                        p_new = (1.0D0 - omega) * pprime2(i, j) + omega * p_new
                        
                        if (abs(p_new - pprime2(i, j)) > p_err) p_err = abs(p_new - pprime2(i, j))
                        
                        pprime2(i, j) = p_new
                    endif
                enddo
            enddo
            !$omp end parallel do

            if (p_err < p_tol) exit
        enddo

        p_mean = 0.0D0
        !$omp parallel do collapse(2) reduction(+:p_mean)
        do j = 1, ny
            do i = 1, nx
                p_mean = p_mean + pprime2(i, j)
            enddo
        enddo
        !$omp end parallel do
        p_mean = p_mean / dble(nx*ny)
        
        !$omp parallel do collapse(2)
        do j = 0, ny+1
            do i = 0, nx+1
                pprime2(i, j) = pprime2(i, j) - p_mean
            enddo
        enddo
        !$omp end parallel do

        ! Final update of velocities and pressure
        !$omp parallel do collapse(2) private(i, j, im1)
        do j = 1, ny
            do i = 1, nx
                
                im1 = i - 1
                if (im1 < 1) im1 = nx
                
                p(i, j) = pstarstar(i, j) + pprime2(i, j)
                u(i, j) = ustarstar(i, j) + ucorr(i, j) + du(i, j) * (pprime2(im1, j) - pprime2(i, j))
                
                if (j > 1) then
                    v(i, j) = vstarstar(i, j) + vcorr(i, j) + dv(i, j) * (pprime2(i, j-1) - pprime2(i, j))
                endif
            enddo
        enddo
        !$omp end parallel do
        
        ! Apply BCs to final velocities and pressure
        !$omp parallel do
        do j = 0, ny+1
            u(0, j) = u(nx, j)
            u(nx+1, j) = u(1, j)
            v(0, j) = v(nx, j)
            v(nx+1, j) = v(1, j)
            p(0, j) = p(nx, j)
            p(nx+1, j) = p(1, j)
        enddo
        !$omp end parallel do
        
        !$omp parallel do
        do i = 0, nx+1
            u(i, 0) = u(i, 1)
            u(i, ny+1) = u(i, ny)
            v(i, 1) = 0.0D0
            v(i, ny+1) = 0.0D0
            v(i, 0) = -v(i, 1)
            p(i, 0) = p(i, 1)
            p(i, ny+1) = p(i, ny)
        enddo
        !$omp end parallel do
        
        ! Calculate final divergence for output
        max_div = 0.0D0
        !$omp parallel do collapse(2) private(i, j, ip1) reduction(max:max_div)
        do j = 1, ny
            do i = 1, nx
                
                ip1 = i + 1 
                if (ip1 > nx) ip1 = 1
                
                bprime(i, j) = rho * dy * u(i, j) - rho * dy * u(ip1, j) + rho * dx * v(i, j) - rho * dx * v(i, j+1)
                if (abs(bprime(i, j)) > max_div) max_div = abs(bprime(i, j))
            enddo
        enddo
        !$omp end parallel do
        
        ! Outputting solutions
        if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
            write(*,'(A,I6,A,F8.4,A,I5,A,ES12.4)') ' Step:', time_step, '  t:', t_now, '  P-iters:', p_iter, '  div(u):', max_div
            write(vtk_filename, '(A, I6.6, A)') 'KH_PISO_', time_step, '.vtk'
            call write_vtk_kh(trim(vtk_filename), u, v, p, tr, nx, ny, dx, dy)
        endif
    enddo

! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) 'Wall time (s): ', end_time - start_time
    
! Deallocate memory
    deallocate(u, ustar, ustarstar, uold)
    deallocate(v, vstar, vstarstar, vold)
    deallocate(p, pstar, pstarstar)
    deallocate(pprime, pprime2, bprime, bprime2)
    deallocate(uprime, vprime, ucorr, vcorr)
    deallocate(du, dv)
    deallocate(tr, tr_old)
    
    stop
end program Kelvin_Helmholtz
    
    
    
! This subroutine aims at outputting VTK files. 
subroutine write_vtk_kh(filename, u, v, p, tr, nx, ny, dx, dy)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny
    real*8,    intent(in) :: dx, dy
    real*8,    intent(in) :: u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), p(0:nx+1, 0:ny+1), tr(0:nx+1, 0:ny+1)
    integer*4 :: i, j, ip1

    open(unit=22, file=filename, status='replace')
    write(22, '(A)') '# vtk DataFile Version 3.0'
    write(22, '(A)') '2D Kelvin-Helmholtz Instability (True PISO)'
    write(22, '(A)') 'ASCII'
    write(22, '(A)') 'DATASET STRUCTURED_POINTS'
    write(22, '(A, I0, 1X, I0, 1X, I0)')          'DIMENSIONS ', nx, ny, 1
    write(22, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ',     0.0D0, 0.0D0, 0.0D0
    write(22, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'SPACING ',    dx, dy, 1.0D0
    write(22, '(A, I0)')                          'POINT_DATA ', nx * ny
    write(22, '(A)') 'SCALARS Tracer float 1'
    write(22, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(22, '(5E15.6)') (tr(i, j), i = 1, nx)
    enddo
    write(22, '(A)') 'SCALARS Pressure float 1'
    write(22, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(22, '(5E15.6)') (p(i, j), i = 1, nx)
    enddo
    write(22, '(A)') 'VECTORS Velocity float'
    do j = 1, ny
        do i = 1, nx
            ip1 = i + 1; if (ip1 > nx) ip1 = 1
            write(22, '(3E15.6)') 0.5D0 * (u(i, j) + u(ip1, j)), 0.5D0 * (v(i, j) + v(i, j+1)), 0.0D0
        enddo
    enddo
    close(22)
end subroutine write_vtk_kh