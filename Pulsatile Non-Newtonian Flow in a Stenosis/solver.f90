program Blood_Stenosis
    use omp_lib
    implicit none
    
! Declaration of variables
    integer*4, parameter :: nx = 200, ny = 40
    integer*4, parameter :: max_step = 10000, output_freq = 50
    integer*4, parameter :: max_p_iter = 5000

    real*8, parameter :: p_tol = 1.0D-7
    real*8, parameter :: dt = 1.0D-4

    ! Blood properties
    real*8, parameter :: rho = 1060.0D0
    real*8, parameter :: mu_inf = 0.0035D0
    real*8, parameter :: mu_zero = 0.056D0
    real*8, parameter :: lambda_C = 3.313D0
    real*8, parameter :: n_C = 0.3568D0

    ! Geometry 
    real*8, parameter :: Lx = 0.05D0
    real*8, parameter :: Ly = 0.01D0

    ! Cosine-bump stenosis
    real*8, parameter :: x_st = 0.5D0 * Lx
    real*8, parameter :: L_st = 0.015D0
    real*8, parameter :: delta_s = 0.5D0      ! This is 50% diameter reduction at throat

    ! Pulsatile inlet
    real*8, parameter :: U_cl_mean = 0.2D0
    real*8, parameter :: U_pulse_amp = 0.5D0
    real*8, parameter :: T_period = 0.5D0

    real*8, parameter :: omega = 1.5D0

    integer*4 :: i, j, p_iter, time_step
    real*8 :: dx, dy, pi
    real*8 :: x_pos, y_pos, occlusion
    real*8 :: U_inlet, t_now
    real*8 :: u_face, v_face
    real*8 :: dudx, dvdy, dudy, dvdx, gamma_dot
    real*8 :: mu_e, mu_w_loc, mu_n, mu_s
    real*8 :: adv, diff_x, diff_y, p_grad
    real*8 :: aw, ae, as, an, ap, b_src, p_new, p_err
    real*8 :: max_div, mass_in, mass_out, U_corr
    real*8 :: start_time, end_time
    integer*4 :: nfluid_in, nfluid_out
    character(len=60) :: vtk_filename

    real*8, allocatable :: u(:,:), us(:,:)
    real*8, allocatable :: v(:,:), vs(:,:)
    real*8, allocatable :: p(:,:), phi(:,:), divus(:,:)
    real*8, allocatable :: mu_cell(:,:)
    integer*4, allocatable :: is_fluid(:,:)
    
! Allocate memory
    allocate(u(0:nx+1, 0:ny+1), us(0:nx+1, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+1), vs(0:nx+1, 0:ny+1))
    allocate(p(0:nx+1, 0:ny+1), phi(0:nx+1, 0:ny+1), divus(0:nx+1, 0:ny+1))
    allocate(mu_cell(0:nx+1, 0:ny+1))
    allocate(is_fluid(0:nx+1, 0:ny+1))

! Defining memory
    pi = acos(-1.0D0)
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)

! Initialization
    u = 0.0D0
    v = 0.0D0
    us = 0.0D0
    vs = 0.0D0
    p = 0.0D0
    phi = 0.0D0
    divus = 0.0D0
    mu_cell = mu_inf
    is_fluid = 0
    ! Applying BCs and ICs
    ! Building the smooth cosine-bump stenosis
    do j = 1, ny
        y_pos = (dble(j) - 0.5D0) * dy
        do i = 1, nx
            x_pos = (dble(i) - 0.5D0) * dx
            occlusion = 0.0D0
            if (abs(x_pos - x_st) < 0.5D0 * L_st) then
                occlusion = (Ly * 0.5D0) * delta_s * 0.5D0 * (1.0D0 + cos(2.0D0 * pi * (x_pos - x_st) / L_st))
            endif
            if (y_pos > occlusion .and. y_pos < Ly - occlusion) is_fluid(i, j) = 1
        enddo
    enddo
    
    !$omp parallel do
    do j = 1, ny
        is_fluid(nx+1, j) = is_fluid(nx, j)
    enddo
    !$omp end parallel do
    
    do j = 1, ny
        y_pos = (dble(j) - 0.5D0) * dy
        if (is_fluid(1, j) == 1) then
            u(1, j) = U_cl_mean * (1.0D0 - ((y_pos - Ly/2.0D0)/(Ly/2.0D0))**2)
        endif
    enddo

! Start the clock
    start_time = omp_get_wtime()

! Start time loop
    do time_step = 1, max_step
        t_now = dble(time_step) * dt

        ! Pulsatile inlet
        U_inlet = U_cl_mean * (1.0D0 + U_pulse_amp * sin(2.0D0*pi*t_now/T_period))
        do j = 1, ny
            y_pos = (dble(j) - 0.5D0) * dy
            if (is_fluid(1, j) == 1) then
                u(1, j) = U_inlet * (1.0D0 - ((y_pos - Ly/2.0D0)/(Ly/2.0D0))**2)
            else
                u(1, j) = 0.0D0
            endif
        enddo

        ! Carreau viscosity
        !$omp parallel do collapse(2) private(i, j, dudx, dvdy, dudy, dvdx, gamma_dot)
        do j = 1, ny
            do i = 1, nx
                if (is_fluid(i, j) == 1) then
                    dudx = (u(i+1, j) - u(i, j)) / dx
                    dvdy = (v(i, j+1) - v(i, j)) / dy
                    
                    dudy = (0.5D0*(u(i, j+1) + u(i+1, j+1)) - 0.5D0*(u(i, j-1) + u(i+1, j-1))) / (2.0D0*dy)
                    dvdx = (0.5D0*(v(i+1, j) + v(i+1, j+1)) - 0.5D0*(v(i-1, j) + v(i-1, j+1))) / (2.0D0*dx)
                    
                    gamma_dot = sqrt(2.0D0*dudx**2 + 2.0D0*dvdy**2 + (dudy + dvdx)**2)
                    
                    mu_cell(i, j) = mu_inf + (mu_zero - mu_inf) * &
                                    (1.0D0 + (lambda_C * gamma_dot)**2)**((n_C - 1.0D0)*0.5D0)
                else
                    mu_cell(i, j) = mu_inf
                endif
            enddo
        enddo
        !$omp end parallel do

        ! u-momentum predictor
        !$omp parallel do collapse(2) private(i, j, mu_e, mu_w_loc, mu_n, mu_s, v_face, adv, diff_x, diff_y, p_grad)
        do j = 1, ny
            do i = 2, nx
                if (is_fluid(i, j) == 1 .and. is_fluid(i-1, j) == 1) then
                    
                    v_face = 0.25D0 * (v(i-1, j) + v(i, j) + v(i-1, j+1) + v(i, j+1))
                    
                    if (u(i, j) > 0.0D0) then
                        adv = u(i, j) * (u(i, j) - u(i-1, j)) / dx
                    else
                        adv = u(i, j) * (u(i+1, j) - u(i, j)) / dx
                    endif
                    
                    if (v_face > 0.0D0) then
                        adv = adv + v_face * (u(i, j) - u(i, j-1)) / dy
                    else
                        adv = adv + v_face * (u(i, j+1) - u(i, j)) / dy
                    endif
                    
                    mu_e = mu_cell(i, j) 
                    mu_w_loc = mu_cell(i-1, j)
                    diff_x = (mu_e * (u(i+1, j) - u(i, j)) - mu_w_loc * (u(i, j) - u(i-1, j))) / (dx*dx)
                    
                    mu_n = 0.25D0 * (mu_cell(i-1, j) + mu_cell(i, j) + mu_cell(i-1, j+1) + mu_cell(i, j+1))
                    mu_s = 0.25D0 * (mu_cell(i-1, j-1) + mu_cell(i, j-1) + mu_cell(i-1, j) + mu_cell(i, j))
                    diff_y = (mu_n * (u(i, j+1) - u(i, j)) - mu_s * (u(i, j) - u(i, j-1))) / (dy*dy)
                    
                    p_grad = (p(i, j) - p(i-1, j)) / dx
                    
                    us(i, j) = u(i, j) + dt * (-adv + (diff_x + diff_y - p_grad)/rho)
                else
                    us(i, j) = 0.0D0
                endif
            enddo
        enddo
        !$omp end parallel do

        ! v-momentum predictor
        !$omp parallel do collapse(2) private(i, j, mu_e, mu_w_loc, mu_n, mu_s, u_face, adv, diff_x, diff_y, p_grad)
        do j = 2, ny
            do i = 1, nx
                if (is_fluid(i, j) == 1 .and. is_fluid(i, j-1) == 1) then
                    
                    u_face = 0.25D0 * (u(i, j-1) + u(i+1, j-1) + u(i, j) + u(i+1, j))
                    
                    if (u_face > 0.0D0) then
                        adv = u_face * (v(i, j) - v(i-1, j)) / dx
                    else
                        adv = u_face * (v(i+1, j) - v(i, j)) / dx
                    endif
                    
                    if (v(i, j) > 0.0D0) then
                        adv = adv + v(i, j) * (v(i, j) - v(i, j-1)) / dy
                    else
                        adv = adv + v(i, j) * (v(i, j+1) - v(i, j)) / dy
                    endif
                    
                    mu_n = mu_cell(i, j)
                    mu_s = mu_cell(i, j-1)           
                    diff_y = (mu_n * (v(i, j+1) - v(i, j)) - mu_s * (v(i, j) - v(i, j-1))) / (dy*dy)
                    
                    mu_e = 0.25D0 * (mu_cell(i, j-1) + mu_cell(i+1, j-1) + mu_cell(i, j) + mu_cell(i+1, j))
                    mu_w_loc = 0.25D0 * (mu_cell(i-1, j-1) + mu_cell(i, j-1) + mu_cell(i-1, j) + mu_cell(i, j))
                    diff_x = (mu_e * (v(i+1, j) - v(i, j)) - mu_w_loc * (v(i, j) - v(i-1, j))) / (dx*dx)
                    
                    p_grad = (p(i, j) - p(i, j-1)) / dy
                    
                    vs(i, j) = v(i, j) + dt * (-adv + (diff_x + diff_y - p_grad)/rho)
                else
                    vs(i, j) = 0.0D0
                endif
            enddo
        enddo
        !$omp end parallel do

        ! BCs on u*, v*
        do j = 1, ny
            y_pos = (dble(j) - 0.5D0) * dy
            if (is_fluid(1, j) == 1) then
                us(1, j) = U_inlet * (1.0D0 - ((y_pos - Ly/2.0D0)/(Ly/2.0D0))**2)
            else
                us(1, j) = 0.0D0
            endif
            vs(0, j) = -vs(1, j); us(nx+1, j) = us(nx, j); vs(nx+1, j) = vs(nx, j)
        enddo
        
        do i = 1, nx
            us(i, 0) = -us(i, 1); us(i, ny+1) = -us(i, ny)
            vs(i, 1) = 0.0D0; vs(i, ny+1) = 0.0D0
        enddo

        ! Mass balance correction at outlet
        mass_in = 0.0D0
        mass_out = 0.0D0
        nfluid_in = 0
        nfluid_out = 0
        
        do j = 1, ny
            if (is_fluid(1, j) == 1) then; mass_in = mass_in + u(1, j); nfluid_in = nfluid_in + 1; endif
            if (is_fluid(1, j) == 1) then; mass_in = mass_in + us(1, j); nfluid_in = nfluid_in + 1; endif
            if (is_fluid(nx+1, j) == 1) then; mass_out = mass_out + us(nx+1, j); nfluid_out = nfluid_out + 1; endif
        enddo
            
        if (nfluid_out > 0) then
            U_corr = (mass_in - mass_out) / dble(nfluid_out)
            do j = 1, ny
                if (is_fluid(nx+1, j) == 1) us(nx+1, j) = us(nx+1, j) + U_corr
            enddo
        endif

        ! Divergence
        max_div = 0.0D0
        !$omp parallel do collapse(2) reduction(max:max_div)
        do j = 1, ny
            do i = 1, nx
                if (is_fluid(i, j) == 1) then
                    
                    divus(i, j) = (us(i+1, j) - us(i, j))/dx + (vs(i, j+1) - vs(i, j))/dy
                    if (abs(divus(i, j)) > max_div) max_div = abs(divus(i, j))
                    
                else
                    divus(i, j) = 0.0D0
                endif
            enddo
        enddo
        !$omp end parallel do

        ! Poisson SOR for phi
        !$omp parallel do collapse(2)
        do i = 0, nx+1
            do j = 0, ny+1
                phi(i, j) = 0.0D0
            enddo
        enddo
        !$omp end parallel do

        do p_iter = 1, max_p_iter
            p_err = 0.0D0
            !$omp parallel do private(i, j, aw, ae, as, an, ap, b_src, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j+1, 2), nx, 2
                    if (is_fluid(i, j) == 1) then
                        
                        aw = dy/dx
                        ae = dy/dx
                        as = dx/dy
                        an = dx/dy
                        
                        if (is_fluid(i-1, j) == 0) aw = 0.0D0
                        if (is_fluid(i+1, j) == 0) ae = 0.0D0
                        if (is_fluid(i, j-1) == 0) as = 0.0D0
                        if (is_fluid(i, j+1) == 0) an = 0.0D0
                        ap = aw + ae + as + an
                        
                        b_src = (rho*dx*dy/dt) * divus(i, j)
                        
                        if (ap > 1.0D-15) then
                            p_new = (ae*phi(i+1, j) + aw*phi(i-1, j) + an*phi(i, j+1) + as*phi(i, j-1) - b_src) / ap
                            p_new = (1.0D0 - omega)*phi(i, j) + omega*p_new
                            
                            if (abs(p_new - phi(i, j)) > p_err) p_err = abs(p_new - phi(i, j))
                            
                            phi(i, j) = p_new
                        endif
                    endif
                enddo
            enddo
            !$omp end parallel do
            !$omp parallel do private(i, j, aw, ae, as, an, ap, b_src, p_new) reduction(max:p_err)
            do j = 1, ny
                do i = 1 + mod(j, 2), nx, 2
                    if (is_fluid(i, j) == 1) then
                        
                        aw = dy/dx
                        ae = dy/dx
                        as = dx/dy
                        an = dx/dy
                        
                        if (is_fluid(i-1, j) == 0) aw = 0.0D0
                        if (is_fluid(i+1, j) == 0) ae = 0.0D0
                        if (is_fluid(i, j-1) == 0) as = 0.0D0
                        if (is_fluid(i, j+1) == 0) an = 0.0D0
                        ap = aw + ae + as + an
                        
                        b_src = (rho*dx*dy/dt) * divus(i, j)
                        
                        if (ap > 1.0D-15) then
                            p_new = (ae*phi(i+1, j) + aw*phi(i-1, j) + an*phi(i, j+1) + as*phi(i, j-1) - b_src) / ap
                            p_new = (1.0D0 - omega)*phi(i, j) + omega*p_new
                            
                            if (abs(p_new - phi(i, j)) > p_err) p_err = abs(p_new - phi(i, j))
                            
                            phi(i, j) = p_new
                        endif
                    endif
                enddo
            enddo
            !$omp end parallel do
            phi(nx+1, :) = 0.0D0
            if (p_err < p_tol) exit
        enddo

        ! Corrector step for u, v
        !$omp parallel do collapse(2)
        do j = 1, ny
            do i = 2, nx
                if (is_fluid(i, j) == 1 .and. is_fluid(i-1, j) == 1) then
                    u(i, j) = us(i, j) - (dt/rho) * (phi(i, j) - phi(i-1, j)) / dx
                else
                    u(i, j) = 0.0D0
                endif
            enddo
        enddo
        !$omp end parallel do
        !$omp parallel do collapse(2)
        do j = 2, ny
            do i = 1, nx
                if (is_fluid(i, j) == 1 .and. is_fluid(i, j-1) == 1) then
                    v(i, j) = vs(i, j) - (dt/rho) * (phi(i, j) - phi(i, j-1)) / dy
                else
                    v(i, j) = 0.0D0
                endif
            enddo
        enddo
        !$omp end parallel do

        ! Update pressure field
        !$omp parallel do collapse(2)
        do j = 1, ny
            do i = 1, nx
                if (is_fluid(i, j) == 1) p(i, j) = p(i, j) + phi(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        ! Applying BCs on u, v, p
        do j = 1, ny
            y_pos = (dble(j) - 0.5D0) * dy
            if (is_fluid(1, j) == 1) then
                u(1, j) = U_inlet * (1.0D0 - ((y_pos - Ly/2.0D0)/(Ly/2.0D0))**2)
            else
                u(1, j) = 0.0D0
            endif
            v(0, j) = -v(1, j)
            u(nx+1, j) = u(nx, j)
            v(nx+1, j) = v(nx, j)
            p(nx+1, j) = 0.0D0
        enddo
        do i = 1, nx
            u(i, 0) = -u(i, 1)
            u(i, ny+1) = -u(i, ny)
            v(i, 1) = 0.0D0
            v(i, ny+1) = 0.0D0
        enddo

        ! Final outlet mass correction
        mass_in = 0.0D0
        mass_out = 0.0D0
        nfluid_out = 0
        do j = 1, ny
            if (is_fluid(1, j) == 1) mass_in = mass_in + u(1, j)
            if (is_fluid(nx+1, j) == 1) then; mass_out = mass_out + u(nx+1, j); nfluid_out = nfluid_out + 1; endif
        enddo
            
        if (nfluid_out > 0) then
            U_corr = (mass_in - mass_out) / dble(nfluid_out)
            do j = 1, ny
                if (is_fluid(nx+1, j) == 1) u(nx+1, j) = u(nx+1, j) + U_corr
            enddo
        endif

        write(*,'(A,I6,A,F8.5,A,F8.4,A,I5,A,ES12.4)') ' Step:', time_step, '  t:', t_now, '  U_in:', U_inlet, &
                '  P-iters:', p_iter, '  div(u*):', max_div
        
        if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
            write(vtk_filename, '(A, I6.6, A)') 'STENO_2D_', time_step, '.vtk'
            call write_vtk_steno(trim(vtk_filename), u, v, p, mu_cell, is_fluid, nx, ny, dx, dy)
        endif
    enddo

! Stop the clock   
    end_time = omp_get_wtime()
    write(*,*) 'Wall time (s) : ', end_time - start_time

! Deallocate memory
    deallocate(u, us, v, vs, p, phi, divus, mu_cell, is_fluid)
    
    stop
end program Blood_Stenosis

            
! This subroutine aims at outputting VTK files.   
subroutine write_vtk_steno(filename, u, v, p, mu_cell, is_fluid, nx, ny, dx, dy)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny
    real*8,    intent(in) :: dx, dy
    real*8,    intent(in) :: u(0:nx+1, 0:ny+1), v(0:nx+1, 0:ny+1), p(0:nx+1, 0:ny+1), mu_cell(0:nx+1, 0:ny+1)
    integer*4, intent(in) :: is_fluid(0:nx+1, 0:ny+1)
    integer*4 :: i, j
    open(unit=23, file=filename, status='replace')
    write(23, '(A)') '# vtk DataFile Version 3.0'
    write(23, '(A)') '2D Pulsatile Non-Newtonian Stenosis'
    write(23, '(A)') 'ASCII'
    write(23, '(A)') 'DATASET STRUCTURED_POINTS'
    write(23, '(A, I0, 1X, I0, 1X, I0)')          'DIMENSIONS ', nx, ny, 1
    write(23, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'ORIGIN ',     0.0D0, 0.0D0, 0.0D0
    write(23, '(A, F10.6, 1X, F10.6, 1X, F10.6)') 'SPACING ',    dx, dy, 1.0D0
    write(23, '(A, I0)')                          'POINT_DATA ', nx * ny
    write(23, '(A)') 'SCALARS Pressure float 1'
    write(23, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(23, '(5E15.6)') (p(i, j), i = 1, nx)
    enddo
    write(23, '(A)') 'SCALARS Viscosity float 1'
    write(23, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(23, '(5E15.6)') (mu_cell(i, j), i = 1, nx)
    enddo
    write(23, '(A)') 'SCALARS Solid float 1'
    write(23, '(A)') 'LOOKUP_TABLE default'
    do j = 1, ny
        write(23, '(5E15.6)') (dble(1 - is_fluid(i, j)), i = 1, nx)
    enddo
    write(23, '(A)') 'VECTORS Velocity float'
    do j = 1, ny
        do i = 1, nx
            write(23, '(3E15.6)') 0.5D0 * (u(i, j) + u(i+1, j)), 0.5D0 * (v(i, j) + v(i, j+1)), 0.0D0
        enddo
    enddo
    close(23)
end subroutine write_vtk_steno