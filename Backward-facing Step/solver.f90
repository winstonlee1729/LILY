program BFS_SIMPLE
    use omp_lib
    implicit none

! Declaration of variables
    ! Grid size and iteration numbers
    integer*4 nx, ny
    integer*4 i, j
    integer*4 iter, maxiter, p_iter
    integer*4 time_step, max_step, output_freq
    parameter(nx = 150, ny = 60)
    parameter(maxiter = 500, max_step = 5000)
    parameter(output_freq = 50)

    ! For program convergence
    real*8, parameter :: tol = 1.0D-5
    real*8 start_time, end_time

    ! Flow field parameters
    real*8 rho, mu, U_e, delta_t
    real*8 Lx, Ly, volume, dx, dy
    parameter(U_e = 1.0D0, rho = 1.0D0, mu = 0.005D0)
    parameter(delta_t = 0.02D0)

    ! Under-relaxation
    real*8 alph_u, alph_v, alph_p
    parameter(alph_u = 0.7D0, alph_v = 0.7D0, alph_p = 0.3D0)

    ! SOR Acceleration
    real*8, parameter :: omega = 1.3D0
    real*8 p_err, p_new

    ! Equation coefficients
    real*8 aw, ae, an, as, ap, ap_o
    real*8 Fw, Fe, Fn, Fs, delta_F
    real*8 Dw, De, Dn, Ds, bmax

    ! Mass conservation variables
    real*8 mass_in, mass_out, U_corr

    ! For outputting solution
    character(len=50) :: vtk_filename
    
    ! Temporal coefficients
    real*8 aw_s, ae_s, an_s, as_s, ap_s, ap_prime
    real*8 w_node, e_node, n_node, s_node, p_past
    ! Theta for temporal discretization scheme
    ! 1.0 = Fully-implicit | 0.5 = Crank-Nicolson | 0.0 = Explicit
    real*8, parameter :: theta = 1.0D0

    ! Field matrices
    real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), du(:,:), un(:,:)
    real*8, allocatable :: v(:,:), vstar(:,:), vold(:,:), dv(:,:), vn(:,:)
    real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:), bprime(:,:)

    ! Geometry
    integer*4, allocatable :: is_fluid(:,:)
    integer*4 step_length, step_height
    real*8 h_channel, y_coord, y_center

! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), uold(0:nx+1, 0:ny+1), du(0:nx+1, 0:ny+1), un(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1), vn(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1))
    allocate(is_fluid(0:nx+1, 0:ny+1))

! Defining variables
    dx = 0.05D0
    dy = 0.05D0
    Lx = dble(nx) * dx
    Ly = dble(ny) * dy
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
    un = 0.0D0 
    vn = 0.0D0
    du = 0.0D0
    dv = 0.0D0
    pprime = 0.0D0
    bprime = 0.0D0

    ! Block geometry
    is_fluid = 0
    step_height = ny / 2
    step_length = nx / 5

    ! Open the main fluid domain
    do i = 1, nx
        do j = 1, ny
            is_fluid(i, j) = 1
        enddo
    enddo

    ! For outlet
    do j = 1, ny
        is_fluid(nx+1, j) = 1
    enddo

    ! Block the physical step
    do j = 1, step_height
        do i = 1, step_length
            is_fluid(i, j) = 0
        enddo
    enddo

    ! Applying initial inlet profile
    do j = step_height + 1, ny
        h_channel = dble(ny - step_height) * dy
        y_center  = dble(ny + step_height) / 2.0D0 * dy
        y_coord   = (dble(j) - 0.5D0) * dy - y_center
        u(1, j) = U_e * (1.0D0 - (2.0D0 * y_coord / h_channel)**2)
    enddo

! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step

        !$omp parallel do collapse(2)
        do i = 0, nx+1
            do j = 0, ny+1
                un(i, j) = u(i, j)
                vn(i, j) = v(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        ! SIMPLE algo. starts here.
        do iter = 1, maxiter

            !$omp parallel do collapse(2)
            do i = 0, nx+1
                do j = 0, ny+1
                    uold(i, j) = u(i, j)
                    vold(i, j) = v(i, j)
                enddo
            enddo
            !$omp end parallel do


            ! u-momentum predictor
            !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, delta_F, aw_s, ae_s, as_s, an_s, ap_s, ap_prime, ap, &
            !$omp w_node, e_node, s_node, n_node, p_past)
            do i = 2, nx
                do j = 1, ny
                    if (is_fluid(i, j) == 1 .and. is_fluid(i-1, j) == 1) then
                        Fe = 0.5D0 * rho * dy * (u(i, j) + u(i+1, j))
                        Fw = 0.5D0 * rho * dy * (u(i-1, j) + u(i, j))
                        Fn = 0.5D0 * rho * dx * (v(i-1, j+1) + v(i, j+1))
                        Fs = 0.5D0 * rho * dx * (v(i-1, j) + v(i, j))
                        delta_F = Fe - Fw + Fn - Fs

                        aw_s = Dw + max(Fw, 0.0D0); ae_s = De + max(-Fe, 0.0D0)
                        as_s = Ds + max(Fs, 0.0D0); an_s = Dn + max(-Fn, 0.0D0)

                        ap_s = aw_s + ae_s + as_s + an_s + delta_F
                        ap = (ap_o + theta * ap_s) / alph_u
                        ap_prime = ap_o - (1.0D0 - theta) * ap_s

                        du(i, j) = dy / ap

                        w_node = aw_s * (theta * uold(i-1, j) + (1.0D0 - theta) * un(i-1, j))
                        e_node = ae_s * (theta * uold(i+1, j) + (1.0D0 - theta) * un(i+1, j))
                        n_node = an_s * (theta * uold(i, j+1) + (1.0D0 - theta) * un(i, j+1))
                        s_node = as_s * (theta * uold(i, j-1) + (1.0D0 - theta) * un(i, j-1))
                        p_past = ap_prime * un(i, j)

                        ustar(i, j) = (w_node + e_node + n_node + s_node + p_past + &
                                        (1.0D0 - alph_u) * ap * uold(i, j) + &
                                        (pstar(i-1, j) - pstar(i, j)) * dy) / ap
                    else
                        ustar(i, j) = 0.0D0
                        du(i, j) = 0.0D0
                    endif
                enddo
            enddo
            !$omp end parallel do

            ! v-momentum predictor (1st-Order Upwind)
            !$omp parallel do collapse(2) private(i, j, Fe, Fw, Fn, Fs, delta_F, aw_s, ae_s, as_s, an_s, ap_s, ap_prime, ap, &
            !$omp w_node, e_node, s_node, n_node, p_past)
            do i = 1, nx
                do j = 2, ny
                    if (is_fluid(i, j) == 1 .and. is_fluid(i, j-1) == 1) then
                        Fe = 0.5D0 * rho * dy * (u(i+1, j-1) + u(i+1, j))
                        Fw = 0.5D0 * rho * dy * (u(i, j-1) + u(i, j))
                        Fn = 0.5D0 * rho * dx * (v(i, j) + v(i, j+1))
                        Fs = 0.5D0 * rho * dx * (v(i, j-1) + v(i, j))
                        delta_F = Fe - Fw + Fn - Fs

                        aw_s = Dw + max(Fw, 0.0D0); ae_s = De + max(-Fe, 0.0D0)
                        as_s = Ds + max(Fs, 0.0D0); an_s = Dn + max(-Fn, 0.0D0)

                        ap_s = aw_s + ae_s + as_s + an_s + delta_F
                        ap = (ap_o + theta * ap_s) / alph_v
                        ap_prime = ap_o - (1.0D0 - theta) * ap_s

                        dv(i, j) = dx / ap

                        w_node = aw_s * (theta * vold(i-1, j) + (1.0D0 - theta) * vn(i-1, j))
                        e_node = ae_s * (theta * vold(i+1, j) + (1.0D0 - theta) * vn(i+1, j))
                        n_node = an_s * (theta * vold(i, j+1) + (1.0D0 - theta) * vn(i, j+1))
                        s_node = as_s * (theta * vold(i, j-1) + (1.0D0 - theta) * vn(i, j-1))
                        p_past = ap_prime * vn(i, j)

                        vstar(i, j) = (w_node + e_node + n_node + s_node + p_past + &
                                        (1.0D0 - alph_v) * ap * vold(i, j) + &
                                        (pstar(i, j-1) - pstar(i, j)) * dx) / ap
                    else
                        vstar(i, j) = 0.0D0
                        dv(i, j) = 0.0D0
                    endif
                enddo
            enddo
            !$omp end parallel do


            ! Applying boundary conditions and mass conservation
            !$omp parallel do private(j, y_coord, h_channel, y_center)
            do j = 1, ny
                ! Inlet
                if (j > step_height) then
                    h_channel = dble(ny - step_height) * dy
                    y_center  = dble(ny + step_height) / 2.0D0 * dy
                    y_coord   = (dble(j) - 0.5D0) * dy - y_center
                    ustar(1, j) = U_e * (1.0D0 - (2.0D0 * y_coord / h_channel)**2)
                else
                    ustar(1, j) = 0.0D0
                endif
                vstar(0, j) = -vstar(1, j)

                ! Right wall 
                ustar(nx+1, j) = ustar(nx, j)
                vstar(nx+1, j) = vstar(nx, j)
            enddo
            !$omp end parallel do

            ! Top and bottom wall
            !$omp parallel do private(i)
            do i = 1, nx
                ustar(i, 0) = -ustar(i, 1)
                ustar(i, ny+1) = -ustar(i, ny)
                vstar(i, 1) = 0.0D0        
                vstar(i, ny+1) = 0.0D0

                ! Step top wall 
                if (i <= step_length) then
                    ustar(i, step_height) = -ustar(i, step_height+1)
                endif
            enddo
            !$omp end parallel do

            ! Step right wall
            do j = 1, step_height
                vstar(step_length, j) = -vstar(step_length+1, j)
            enddo

            ! For mass conservation
            mass_in = 0.0D0
            mass_out = 0.0D0
            do j = step_height + 1, ny
                mass_in = mass_in + ustar(1, j)
            enddo
            do j = 1, ny
                mass_out = mass_out + ustar(nx+1, j)
            enddo

            U_corr = (mass_in - mass_out) / dble(ny)
            do j = 1, ny
                ustar(nx+1, j) = ustar(nx+1, j) + U_corr
            enddo

            do j = 1, ny
                du(nx+1, j) = du(nx, j)
            enddo

            ! For mass source bprime
            bmax = 0.0D0
            !$omp parallel do collapse(2) reduction(max:bmax)
            do i = 1, nx
                do j = 1, ny
                    if (is_fluid(i, j) == 1) then
                        bprime(i, j) = (rho * dy * ustar(i, j)) - (rho * dy * ustar(i+1, j)) &
                                     + (rho * dx * vstar(i, j)) - (rho * dx * vstar(i, j+1))
                        if (abs(bprime(i, j)) > bmax) bmax = abs(bprime(i, j))
                    else
                        bprime(i, j) = 0.0D0
                        pprime(i, j) = 0.0D0
                    endif
                enddo
            enddo
            !$omp end parallel do

            ! Red-Black SOR for pressure correction equation
            do p_iter = 1, 500
                p_err = 0.0D0

                ! Red Cells
                !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do i = 1, nx
                    do j = 1 + mod(i+1, 2), ny, 2
                        if (is_fluid(i, j) == 1) then
                            aw = rho * dy * du(i, j)
                            ae = rho * dy * du(i+1, j)
                            as = rho * dx * dv(i, j)
                            an = rho * dx * dv(i, j+1)
                            
                            if (is_fluid(i-1, j) == 0) aw = 0.0D0
                            if (is_fluid(i+1, j) == 0) ae = 0.0D0
                            if (is_fluid(i, j-1) == 0) as = 0.0D0
                            if (is_fluid(i, j+1) == 0) an = 0.0D0
                            
                            ap = aw + ae + an + as
                            
                            if (ap > 1.0D-15) then
                                p_new = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + an * pprime(i, j+1) + &
                                        as * pprime(i, j-1) + bprime(i, j)) / ap
                                p_new = (1.0D0 - omega) * pprime(i, j) + omega * p_new
                                
                                if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                                
                                pprime(i, j) = p_new
                                
                            endif
                        endif
                    enddo
                enddo
                !$omp end parallel do

                ! Black Cells
                !$omp parallel do private(i, j, aw, ae, as, an, ap, p_new) reduction(max:p_err)
                do i = 1, nx
                    do j = 1 + mod(i, 2), ny, 2
                        if (is_fluid(i, j) == 1) then
                            aw = rho * dy * du(i, j)
                            ae = rho * dy * du(i+1, j)
                            as = rho * dx * dv(i, j)
                            an = rho * dx * dv(i, j+1)
                            
                            if (is_fluid(i-1, j) == 0) aw = 0.0D0
                            if (is_fluid(i+1, j) == 0) ae = 0.0D0
                            if (is_fluid(i, j-1) == 0) as = 0.0D0
                            if (is_fluid(i, j+1) == 0) an = 0.0D0
                            
                            ap = aw + ae + an + as
                            
                            if (ap > 1.0D-15) then
                                p_new = (ae * pprime(i+1, j) + aw * pprime(i-1, j) + an * pprime(i, j+1) + &
                                        as * pprime(i, j-1) + bprime(i, j)) / ap
                                p_new = (1.0D0 - omega) * pprime(i, j) + omega * p_new
                                
                                if (abs(p_new - pprime(i, j)) > p_err) p_err = abs(p_new - pprime(i, j))
                                
                                pprime(i, j) = p_new
                                
                            endif
                        endif
                    enddo
                enddo
                !$omp end parallel do

                pprime(nx+1, :) = 0.0D0
                
                if (p_err < tol) exit
            enddo

            ! Correcting velocities
            !$omp parallel do collapse(2) private(i, j)
            do i = 1, nx
                do j = 1, ny
                    if (is_fluid(i, j) == 1) then
                        p(i, j) = pstar(i, j) + alph_p * pprime(i, j)
                        pstar(i, j) = p(i, j)
                    endif

                    if (i > 1) then
                        if (is_fluid(i, j) == 1 .and. is_fluid(i-1, j) == 1) then
                            u(i, j) = ustar(i, j) + du(i, j) * (pprime(i-1, j) - pprime(i, j))
                        endif
                    endif

                    if (j > 1) then
                        if (is_fluid(i, j) == 1 .and. is_fluid(i, j-1) == 1) then
                            v(i, j) = vstar(i, j) + dv(i, j) * (pprime(i, j-1) - pprime(i, j))
                        endif
                    endif
                enddo
            enddo
            !$omp end parallel do

            ! Reapplying boundary conditions
            !$omp parallel do private(j, y_coord, h_channel, y_center)
            do j = 1, ny
                
                if (j > step_height) then
                    h_channel = dble(ny - step_height) * dy
                    y_center  = dble(ny + step_height) / 2.0D0 * dy
                    y_coord   = (dble(j) - 0.5D0) * dy - y_center
                    u(1, j) = U_e * (1.0D0 - (2.0D0 * y_coord / h_channel)**2)
                else
                    u(1, j) = 0.0D0
                endif
                
                v(0, j) = -v(1, j)

                u(nx+1, j) = u(nx, j)
                v(nx+1, j) = v(nx, j)

                p(nx+1, j) = 0.0D0
                pstar(nx+1, j) = 0.0D0
            enddo
            !$omp end parallel do

            !$omp parallel do private(i)
            do i = 1, nx
                u(i, 0) = -u(i, 1)
                u(i, ny+1) = -u(i, ny)
                v(i, 1) = 0.0D0    
                v(i, ny+1) = 0.0D0
                
                if (i <= step_length) then
                    u(i, step_height) = -u(i, step_height+1)
                endif
                
            enddo
            !$omp end parallel do

            do j = 1, step_height
                v(step_length, j) = -v(step_length+1, j)
            enddo

            mass_in = 0.0D0
            mass_out = 0.0D0
            do j = step_height + 1, ny
                mass_in = mass_in + u(1, j)
            enddo
            do j = 1, ny
                mass_out = mass_out + u(nx+1, j)
            enddo
            U_corr = (mass_in - mass_out) / dble(ny)
            do j = 1, ny
                u(nx+1, j) = u(nx+1, j) + U_corr
            enddo

            if (bmax <= tol) exit
            
        enddo       ! End of SIMPLE algo.
        

        write(*,*) 'Time step: ', time_step, '| Inner Iters: ', iter, '| Max mass error: ', bmax

        ! Outputting solution
        if (mod(time_step, output_freq) == 0) then
            write(vtk_filename, '(A, I6.6, A)') 'BFS_2D_', time_step, '.vtk'
            call write_vtk(trim(vtk_filename), u, v, p, nx, ny, dx, dy)
        endif

    enddo   ! End of time loop

! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total Time Steps: ', max_step
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'
    
! Deallocate memory
    deallocate(u, ustar, uold, du, un)
    deallocate(v, vstar, vold, dv, vn)
    deallocate(p, pstar, pprime, bprime, is_fluid)
    
    stop
end program BFS_SIMPLE

    
    
! This subroutine aims at outputting VTK files
subroutine write_vtk(filename, u, v, p, nx, ny, dx, dy)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny
    real*8, intent(in) :: dx, dy
    real*8, intent(in) :: u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), p(0:nx+1, 0:ny+1)
    integer*4 :: i, j

    open(unit=20, file=filename, status='replace')
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') '2D CFD Transient Backward Facing Step'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET STRUCTURED_POINTS'

    write(20, '(A, I0, 1X, I0, 1X, I0)') 'DIMENSIONS ', nx, ny, 1
    write(20, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ', 0.0D0, 0.0D0, 0.0D0
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
            write(20, '(3E15.6)') 0.5D0 * (u(i, j) + u(i+1, j)), &
                                  0.5D0 * (v(i, j) + v(i, j+1)), &
                                  0.0D0
        enddo
    enddo
    close(20)
end subroutine write_vtk