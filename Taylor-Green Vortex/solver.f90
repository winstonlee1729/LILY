program TGV_SIMPLE
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration numbers
    integer*4 nx, ny
    integer*4 i, j
    integer*4 iter, maxiter, p_iter, max_p_iter
    integer*4 time_step, max_step, output_freq
    parameter(nx = 64, ny = 64)
    parameter(maxiter = 50, max_step = 1000)
    parameter(output_freq = 20)
    
    ! For program convergence and time
    real*8, parameter :: tol = 1.0D-7
    real*8, parameter :: sor_tol = 1.0D-8
    real*8 start_time, end_time
    
    ! Flow field parameters 
    real*8 rho, mu, U_e, delta_t
    real*8 Lx, Ly, volume, dx, dy, pi
    parameter(U_e = 1.0D0, rho = 1.0D0, mu = 0.01D0) 
    parameter(delta_t = 5.0D-3)
    
    ! Under-relaxation 
    real*8 alph_u, alph_v, alph_p
    parameter(alph_u = 0.7D0, alph_v = 0.7D0, alph_p = 0.3D0)
    
    ! SOR Acceleration
    real*8, parameter :: omega = 1.6D0
    real*8 p_err, p_new
    
    ! Equation coefficients
    real*8 aw, ae, an, as, ap, ap_o, ap_s
    real*8 Fw, Fe, Fn, Fs
    real*8 Dw, De, Dn, Ds
    real*8 bp, bmax, b_mean, sor_res
    
    ! Validation variables
    real*8 decay, err_u, err_v
    real*8 x_u, y_u, x_v, y_v, x_p, y_p, u_ex, v_ex
    integer*4 count_n
    
    ! For solution output
    character(len=50) :: vtk_filename
    
    ! Theta for temporal discretization scheme
    ! 1.0 = Fully-implicit | 0.5 = Crank-Nicolson | 0.0 = Explicit
    real*8, parameter :: theta = 1.0D0
    
    ! Flow matrices
    real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), du(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vold(:,:), dv(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:), bprime(:,:)
    
! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1), du(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1)) 
    
! Defining variables
    pi = acos(-1.0D0)
    Lx = 2.0D0 * pi
    Ly = 2.0D0 * pi
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    volume = dx * dy
    
    Dw = mu * (dy / dx)
    De = mu * (dy / dx)
    Dn = mu * (dx / dy)
    Ds = mu * (dx / dy)
    ap_o = rho * (volume / delta_t)
    
! Initialization 
    u = 0.0D0
    v = 0.0D0
    p = 0.0D0
    ustar = 0.0D0
    vstar = 0.0D0
    pstar = 0.0D0
    du = 0.0D0
    dv = 0.0D0
    pprime = 0.0D0 
    bprime = 0.0D0
    
    do j = 0, ny+1
        do i = 0, nx+1
            x_u = dble(i) * dx
            y_u = (dble(j) - 0.5D0) * dy
            u(i, j) = U_e * cos(x_u) * sin(y_u)

            x_v = (dble(i) - 0.5D0) * dx
            y_v = dble(j) * dy
            v(i, j) = -U_e * sin(x_v) * cos(y_v)

            x_p = (dble(i) - 0.5D0) * dx
            y_p = (dble(j) - 0.5D0) * dy
            p(i, j) = -0.25D0 * rho * U_e**2 * (cos(2.0D0*x_p) + cos(2.0D0*y_p))
        enddo
    enddo
    ! Applying initial periodic boundaries
    u(0, :) = u(nx, :); u(nx+1, :) = u(1, :)
    u(:, 0) = u(:, ny); u(:, ny+1) = u(:, 1)
    v(0, :) = v(nx, :); v(nx+1, :) = v(1, :)
    v(:, 0) = v(:, ny); v(:, ny+1) = v(:, 1)
    p(0, :) = p(nx, :); p(nx+1, :) = p(1, :)
    p(:, 0) = p(:, ny); p(:, ny+1) = p(:, 1)
    pstar = p

! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step
        
        !$omp parallel do collapse(2)
        do i = 0, nx+1
            do j = 0, ny+1
                uold(i, j) = u(i, j)
                vold(i, j) = v(i, j)
            enddo
        enddo
        !$omp end parallel do
    
        do iter = 1, maxiter
            
            ! u-momentum predictor
            !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, aw, ae, as, an, ap_s, ap)
            do i = 1, nx
                do j = 1, ny
                    Fw = rho * dy * 0.5D0 * (u(i-1, j) + u(i, j))
                    Fe = rho * dy * 0.5D0 * (u(i, j) + u(i+1, j))
                    Fs = rho * dx * 0.5D0 * (v(i, j-1) + v(i+1, j-1))
                    Fn = rho * dx * 0.5D0 * (v(i, j) + v(i+1, j))

                    aw = max( Fw, 0.0D0) + Dw
                    ae = max(-Fe, 0.0D0) + De
                    as = max( Fs, 0.0D0) + Ds
                    an = max(-Fn, 0.0D0) + Dn

                    ap_s = aw + ae + as + an + (Fe - Fw) + (Fn - Fs)
                    ap   = (ap_o + theta * ap_s) / alph_u

                    du(i, j) = dy / ap
                    ustar(i, j) = (aw * u(i-1, j) + ae * u(i+1, j) &
                                 + as * u(i, j-1) + an * u(i, j+1) &
                                 + ap_o * uold(i, j) &
                                 + (1.0D0 - alph_u) * ap * u(i, j) &
                                 + (pstar(i, j) - pstar(i+1, j)) * dy) / ap
                enddo
            enddo
            !$omp end parallel do
            
            ! Periodic wrap
            ustar(0, :) = ustar(nx, :)
            ustar(nx+1, :) = ustar(1, :)
            ustar(:, 0) = ustar(:, ny)
            ustar(:, ny+1) = ustar(:, 1)
            du(0, :) = du(nx, :)
            du(nx+1, :) = du(1, :)
            du(:, 0) = du(:, ny)
            du(:, ny+1) = du(:, 1)

            ! v-momentum predictor
            !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, aw, ae, as, an, ap_s, ap)
            do i = 1, nx
                do j = 1, ny
                    Fw = rho * dy * 0.5D0 * (u(i-1, j) + u(i-1, j+1))
                    Fe = rho * dy * 0.5D0 * (u(i, j) + u(i, j+1))
                    Fs = rho * dx * 0.5D0 * (v(i, j-1) + v(i, j))
                    Fn = rho * dx * 0.5D0 * (v(i, j) + v(i, j+1))

                    aw = max( Fw, 0.0D0) + Dw
                    ae = max(-Fe, 0.0D0) + De
                    as = max( Fs, 0.0D0) + Ds
                    an = max(-Fn, 0.0D0) + Dn

                    ap_s = aw + ae + as + an + (Fe - Fw) + (Fn - Fs)
                    ap   = (ap_o + theta * ap_s) / alph_v

                    dv(i, j) = dx / ap
                    vstar(i, j) = (aw * v(i-1, j) + ae * v(i+1, j) &
                                 + as * v(i, j-1) + an * v(i, j+1) &
                                 + ap_o * vold(i, j) &
                                 + (1.0D0 - alph_v) * ap * v(i, j) &
                                 + (pstar(i, j) - pstar(i, j+1)) * dx) / ap
                enddo
            enddo
            !$omp end parallel do
            
            ! Periodic wrap
            vstar(0, :) = vstar(nx, :)
            vstar(nx+1, :) = vstar(1, :)
            vstar(:, 0) = vstar(:, ny)
            vstar(:, ny+1) = vstar(:, 1)
            dv(0, :) = dv(nx, :)
            dv(nx+1, :) = dv(1, :)
            dv(:, 0) = dv(:, ny)
            dv(:, ny+1) = dv(:, 1)


            ! Mass source and mean subtraction
            pprime = 0.0D0
            b_mean = 0.0D0
            bmax   = 0.0D0
            
            !$omp parallel do collapse(2) private(i, j, bp) reduction(+:b_mean) reduction(max:bmax)
            do i = 1, nx
                do j = 1, ny
                    bp = -rho * ((ustar(i, j) - ustar(i-1, j)) * dy &
                               + (vstar(i, j) - vstar(i, j-1)) * dx)
                    b_mean = b_mean + bp
                    if (abs(bp) > bmax) bmax = abs(bp)
                enddo
            enddo
            !$omp end parallel do
            
            b_mean = b_mean / dble(nx * ny)

            ! Red-Black SOR for pressure correction equation
            do p_iter = 1, max_p_iter
                sor_res = 0.0D0

                ! Red Cells
                !$omp parallel do private(i, j, aw, ae, as, an, ap, bp) reduction(max:sor_res)
                do i = 1, nx
                    do j = 1 + mod(i+1, 2), ny, 2
                        aw = rho * dy * du(i-1, j)
                        ae = rho * dy * du(i, j)
                        as = rho * dx * dv(i, j-1)
                        an = rho * dx * dv(i, j)
                        ap = aw + ae + as + an
                        
                        bp = -rho * ((ustar(i, j) - ustar(i-1, j)) * dy &
                                   + (vstar(i, j) - vstar(i, j-1)) * dx) - b_mean
                                   
                        bprime(i, j) = (aw * pprime(i-1, j) + ae * pprime(i+1, j) &
                                      + as * pprime(i, j-1) + an * pprime(i, j+1) &
                                      + bp - ap * pprime(i, j)) / ap
                                      
                        if (abs(bprime(i, j)) > sor_res) sor_res = abs(bprime(i, j))
                        pprime(i, j) = pprime(i, j) + omega * bprime(i, j)
                    enddo
                enddo
                !$omp end parallel do
                
                ! Wrap red cell update
                pprime(0, :) = pprime(nx, :)
                pprime(nx+1, :) = pprime(1, :)
                pprime(:, 0) = pprime(:, ny)
                pprime(:, ny+1) = pprime(:, 1)

                ! Black Cells
                !$omp parallel do private(i, j, aw, ae, as, an, ap, bp) reduction(max:sor_res)
                do i = 1, nx
                    do j = 1 + mod(i, 2), ny, 2
                        aw = rho * dy * du(i-1, j)
                        ae = rho * dy * du(i, j)
                        as = rho * dx * dv(i, j-1)
                        an = rho * dx * dv(i, j)
                        ap = aw + ae + as + an
                        
                        bp = -rho * ((ustar(i, j) - ustar(i-1, j)) * dy &
                                   + (vstar(i, j) - vstar(i, j-1)) * dx) - b_mean
                                   
                        bprime(i, j) = (aw * pprime(i-1, j) + ae * pprime(i+1, j) &
                                      + as * pprime(i, j-1) + an * pprime(i, j+1) &
                                      + bp - ap * pprime(i, j)) / ap
                                      
                        if (abs(bprime(i, j)) > sor_res) sor_res = abs(bprime(i, j))
                        pprime(i, j) = pprime(i, j) + omega * bprime(i, j)
                    enddo
                enddo
                !$omp end parallel do
                
                ! Wrap black update
                pprime(0, :) = pprime(nx, :)
                pprime(nx+1, :) = pprime(1, :)
                pprime(:, 0) = pprime(:, ny)
                pprime(:, ny+1) = pprime(:, 1)

                ! Anchor poisson matrix 
                pprime(1:nx, 1:ny) = pprime(1:nx, 1:ny) - pprime(1, 1)
                
                ! Wrap anchor
                pprime(0, :) = pprime(nx, :)
                pprime(nx+1, :) = pprime(1, :)
                pprime(:, 0) = pprime(:, ny)
                pprime(:, ny+1) = pprime(:, 1)

                if (sor_res < sor_tol) exit
            enddo

            ! Correcting velocities and pressure directly
            !$omp parallel do collapse(2)
            do i = 1, nx
                do j = 1, ny
                    p(i, j) = pstar(i, j) + alph_p * pprime(i, j)
                    u(i, j) = ustar(i, j) - du(i, j) * (pprime(i+1, j) - pprime(i, j))
                    v(i, j) = vstar(i, j) - dv(i, j) * (pprime(i, j+1) - pprime(i, j))
                enddo
            enddo
            !$omp end parallel do
            
            p(0, :) = p(nx, :); p(nx+1, :) = p(1, :)
            p(:, 0) = p(:, ny); p(:, ny+1) = p(:, 1)
            u(0, :) = u(nx, :); u(nx+1, :) = u(1, :)
            u(:, 0) = u(:, ny); u(:, ny+1) = u(:, 1)
            v(0, :) = v(nx, :); v(nx+1, :) = v(1, :)
            v(:, 0) = v(:, ny); v(:, ny+1) = v(:, 1)

            pstar = p

            if (bmax < tol) exit
        enddo  
        
        ! Error validation 
        decay = exp(-2.0D0 * (mu / rho) * dble(time_step) * delta_t)
        err_u = 0.0D0
        err_v = 0.0D0
        count_n = 0
        
        !$omp parallel do collapse(2) private(i, j, x_u, y_u, x_v, y_v, u_ex, v_ex) reduction(+:err_u, err_v, count_n)
        do i = 1, nx
            do j = 1, ny
                x_u = dble(i) * dx
                y_u = (dble(j) - 0.5D0) * dy
                u_ex = U_e * cos(x_u) * sin(y_u) * decay
                err_u = err_u + (u(i, j) - u_ex)**2

                x_v = (dble(i) - 0.5D0) * dx
                y_v = dble(j) * dy
                v_ex = -U_e * sin(x_v) * cos(y_v) * decay
                err_v = err_v + (v(i, j) - v_ex)**2

                count_n = count_n + 1
            enddo
        enddo
        !$omp end parallel do
        
        err_u = sqrt(err_u / dble(count_n))
        err_v = sqrt(err_v / dble(count_n))

        ! Preventing crash
        if (err_u /= err_u .or. err_v /= err_v) then
            write(*,*) 'NaN detected at step ', time_step, ' Exploded.'
            stop 1
        endif

        ! Outputting solution
        if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
            write(*,'(A,I6,A,I3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
                'Step: ', time_step, ' | Iters: ', iter, &
                ' | Max Mass Err: ', bmax, ' | L2_u: ', err_u, ' | L2_v: ', err_v, &
                ' | Decay: ', decay
                
            write(vtk_filename, '(A, I6.6, A)') 'TGV_2D_', time_step, '.vtk'
            call write_vtk(trim(vtk_filename), u, v, p, nx, ny, dx, dy)
        endif

    enddo 
    
! Stop the clock
    end_time = omp_get_wtime()   
    write(*,*) '------------------------------------'
    write(*,'(A,F10.3,A)') 'Total Compute Time: ', end_time - start_time, ' s'
    write(*,*) '------------------------------------'

! Outputting final state
    write(vtk_filename, '(A)') 'TGV_2D_final.vtk'
    call write_vtk(trim(vtk_filename), u, v, p, nx, ny, dx, dy)

! Deallocate memory
    deallocate(u, ustar, uold, du)
    deallocate(v, vstar, vold, dv)
    deallocate(p, pstar, pprime, bprime)

    stop
end program TGV_SIMPLE
    
    
    
! This subroutine aims at outputting VTK solution files.
subroutine write_vtk(filename, u, v, p, nx, ny, dx, dy)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny
    real*8, intent(in) :: dx, dy
    real*8, intent(in) :: u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), p(0:nx+1, 0:ny+1)
    integer*4 :: i, j
    
    open(unit=20, file=filename, status='replace')
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') '2D CFD Transient Taylor-Green Vortex'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET STRUCTURED_POINTS'
    
    write(20, '(A, I0, 1X, I0, 1X, I0)') 'DIMENSIONS ', nx, ny, 1
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ', 0.5D0*dx, 0.5D0*dy, 0.0D0
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'SPACING ', dx, dy, 1.0D0
    write(20, '(A, I0)') 'POINT_DATA ', nx * ny
    
    write(20, '(A)') 'SCALARS Pressure float 1'
    write(20, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(20, '(5E15.6)') (p(i, j), i=1,nx)
    enddo

    write(20, '(A)') 'VECTORS Velocity float'
    do j = 1, ny
        do i = 1, nx
            write(20, '(3E15.6)') 0.5D0 * (u(i-1, j) + u(i, j)), &
                                  0.5D0 * (v(i, j-1) + v(i, j)), &
                                  0.0D0
        enddo
    enddo
    close(20)
end subroutine write_vtk