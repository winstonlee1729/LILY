program RTI_Solver
    use omp_lib
    implicit none

! Declaration of variables
    ! Grid size and iteration number
    integer :: nx, ny
    integer :: i, j
    integer :: time_step, max_step, output_freq, p_max_it
    parameter(nx = 128, ny = 512)
    parameter(max_step = 100000, p_max_it = 4000)
    parameter(output_freq = 300)
    
    ! Flow parameters
    real(8) :: Lx, Ly, dx, dy
    real(8) :: rho_H, rho_L, mu_H, mu_L, g_y, pi
    real(8) :: dt, amp0
    real(8) :: x, y, t
    parameter(Lx = 1.0D0, Ly = 4.0D0)               ! Aspect ratio 1:4
    parameter(rho_H = 3.0D0, rho_L = 1.0D0)         ! Heavy fluid (top), Light fluid (bottom)
    parameter(mu_H = 0.005D0, mu_L = 0.005D0)
    parameter(g_y = -1.0D0)
    parameter(dt = 2.0D-4, amp0 = 0.05D0)           ! Initial perturbation amplitude

    ! For SOR factor
    real(8), parameter :: omega_sor = 1.7d0
    
    ! For pressure Poisson solver convergence
    real(8), parameter :: ptol      = 1.0d-6

    ! Flow field matrics
    real(8), allocatable :: u(:,:),  v(:,:),  p(:,:)
    real(8), allocatable :: us(:,:), vs(:,:)
    real(8), allocatable :: rho(:,:), rho_new(:,:), mu_c(:,:)

    ! For monitoring and output
    real(8) :: divmax, umax, vmax
    
    ! For outputting solution
    character(len=80) :: fname

! Defining variables
    pi = 4.0D0 * atan(1.0D0)
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)

! Allocate memory
    allocate(u(0:nx+2, 0:ny+1))
    allocate(v(0:nx+1, 0:ny+2))
    allocate(p(0:nx+1, 0:ny+1))
    allocate(us(0:nx+2, 0:ny+1))
    allocate(vs(0:nx+1, 0:ny+2))
    allocate(rho(0:nx+1, 0:ny+1))
    allocate(rho_new(0:nx+1, 0:ny+1))
    allocate(mu_c(0:nx+1, 0:ny+1))

! Initialization
    u = 0.0D0 
    v = 0.0D0 
    p = 0.0D0
    us = 0.0D0
    vs = 0.0D0
    rho = rho_L
    rho_new = rho_L 
    mu_c = mu_L
    ! For single mode interface: deepest at x = Lx/2  =>  central spike
    ! y_int(x) = Ly/2 + amp * cos(2*pi*x/Lx)
    do j = 0, ny+1
        y = (dble(j) - 0.5D0) * dy
        
        do i = 0, nx+1
            x = (dble(i) - 0.5D0) * dx
            
            if (y > 0.5D0 * Ly + amp0 * cos(2.0D0 * pi * x / Lx)) then
                rho(i,j)  = rho_H
                mu_c(i,j) = mu_H
            else
                rho(i,j)  = rho_L
                mu_c(i,j) = mu_L
            endif
            
        enddo
    enddo
    rho_new = rho

! Solver information
    write(*,'(A)')         '======================================'
    write(*,'(A)')         '  Rayleigh-Taylor Instability solver  '
    write(*,'(A)')         '======================================'
    write(*,'(A,I5," x ",I5)') '  Grid       : ', nx, ny
    write(*,'(A,F10.4)')   '  Atwood No  : ', (rho_H - rho_L) / (rho_H + rho_L)
    write(*,'(A,F10.4)')   '  dt         : ', dt
    write(*,'(A,I8)')      '  total step : ', max_step

! Start the clock
    t = 0.0d0
    
! Applying boundary conditions
    call apply_bc

! Initialization of output
    call write_vtk(0)

! Start time loop
    do time_step = 1, max_step
        t = t + dt
        
        ! Reappling BCs
        call apply_bc

        ! For advection density 
        call advect_rho

        ! Momentum predictor for u* and v*
        call momentum_predictor

        ! Solve pressure Poisson
        call pressure_solve

        ! Corrector:  u = u* - dt/rho_face * grad(p) 
        call velocity_correct

        ! Swap rho 
        !$omp parallel do collapse(2)
        do j = 0, ny+1
            do i = 0, nx+1
                rho(i, j)  = rho_new(i, j)
                mu_c(i, j) = (rho(i, j) - rho_L) / (rho_H - rho_L) * mu_H + &
                             (rho_H - rho(i, j)) / (rho_H - rho_L) * mu_L
            enddo
        enddo
        !$omp end parallel do

        ! Outputting solution
        if (mod(time_step, output_freq)==0 .or. time_step==1) then
            call max_div(divmax, umax, vmax)
            ! Monitoring convergence and max velocity
            write(*,'(A,I6,A,F8.4,A,ES10.3,A,F8.3,A,F8.3)') &
                    ' step:',time_step,'  t:',t,'  divMax:',divmax,'  uMax:',umax,'  vMax:',vmax
            call write_vtk(time_step)
        endif
        
    enddo ! End of transient loop

    write(*,*) 'Hooray!'

contains


! This subroutine aims at applying the BCs on the ghost cells. 
    subroutine apply_bc
        integer :: i, j
        ! Note that left / right walls : u = 0, v reflected, rho zero-grad
        !$omp parallel do 
        do j = 0, ny+1
            u(1, j) = 0.0D0
            u(nx+1, j) = 0.0D0
            u(0, j) = -u(2, j)           
            u(nx+2, j) = -u(nx, j)
            v(0, j) = -v(1, j)
            v(nx+1, j) = -v(nx, j)
            rho(0, j) = rho(1, j)
            rho(nx+1, j) = rho(nx, j)
            mu_c(0, j) = mu_c(1, j)
            mu_c(nx+1, j) = mu_c(nx, j)
        enddo
        !$omp end parallel do
        ! Note that bottom / top walls: v = 0, u reflected
        !$omp parallel do
        do i = 0, nx+1
            v(i, 1) = 0.0D0
            v(i, ny+1) = 0.0D0
            v(i, 0) = -v(i, 2)
            v(i, ny+2) = -v(i, ny)
            u(i, 0) = -u(i, 1)
            u(i, ny+1) = -u(i, ny)
            rho(i, 0) = rho(i, 1)
            rho(i, ny+1) = rho(i, ny)
            mu_c(i, 0) = mu_c(i, 1)
            mu_c(i, ny+1) = mu_c(i, ny)
        enddo
        !$omp end parallel do
    end subroutine apply_bc

    
! This function aims at calculating the face value of a variable for the convection term, using the QUICK scheme. 
! The function takes in the face velocity (uf) and the values of the variable at the upstream cell (phU), 
! the central cell (phC), and the downstream cell (phD). 
! The function returns the face value (phf) based on the sign of the face velocity, ensuring that the upwind direction is correctly accounted for in the interpolation.
    pure function quick(uf, phU, phC, phD) result(phf)
        real(8), intent(in) :: uf, phU, phC, phD
        real(8) :: phf
        
        if (uf >= 0.0D0) then
            phf = 0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU
        else
            phf = 0.375D0 * phC + 0.75D0 * phD - 0.125D0 * phU  
        endif
        
    end function quick

    
! This subroutine aims at advecting the density field using the QUICK scheme.
    subroutine advect_rho
        integer :: i,j
        real(8) :: uE, uW, vN, vS
        real(8) :: fE, fW, fN, fS
        real(8) :: phC, phD, phU

        !$omp parallel do collapse(2) private(i, j, uE, uW, vN, vS, fE, fW, fN, fS, phC, phD, phU)
        do j = 1, ny
            do i = 1, nx
                uE = u(i+1, j)
                uW = u(i, j)
                vN = v(i, j+1)
                vS = v(i, j)

                ! East face  
                if (uE >= 0.0D0) then
                    phC = rho(i, j)  
                    phD = rho(i+1, j)
                    
                    if (i >= 2) then
                        phU = rho(i-1, j)
                    else
                        phU = rho(i, j) 
                    endif
                    
                else
                    phC = rho(i+1, j)
                    phD = rho(i, j)
                    
                    if (i <= nx-1) then
                        phU = rho(i+2, j)
                    else
                        phU = rho(i+1, j)
                    endif
                    
                endif
                
                fE = uE * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! West face 
                if (uW >= 0.0D0) then
                    phC = rho(i-1, j)
                    phD = rho(i, j)
                    
                    if (i >= 3) then
                        phU = rho(i-2, j)
                    else
                        phU = rho(i-1, j)
                    endif
                    
                else
                    phC = rho(i, j)  
                    phD = rho(i-1, j)
                    
                    if (i <= nx-1) then
                        phU = rho(i+1, j)
                    else
                        phU = rho(i, j)
                    endif
                    
                endif
                
                fW = uW * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! North face 
                if (vN >= 0.0D0) then
                    phC = rho(i, j)
                    phD = rho(i, j+1)
                    
                    if (j >= 2) then
                        phU = rho(i, j-1)
                    else
                        phU = rho(i, j)
                    endif
                    
                else
                    phC = rho(i, j+1)
                    phD = rho(i, j)
                    
                    if (j <= ny-1) then
                        phU = rho(i, j+2)
                    else
                        phU = rho(i, j+1)
                    endif
                    
                endif
                fN = vN * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! South face
                if (vS >= 0.0D0) then
                    phC = rho(i, j-1)
                    phD = rho(i, j)
                    
                    if (j >= 3) then
                        phU = rho(i, j-2)
                    else
                        phU = rho(i, j-1)
                    endif
                    
                else
                    phC = rho(i, j)  
                    phD = rho(i, j-1)
                    
                    if (j <= ny) then
                        phU = rho(i, j+1)
                    else
                        phU = rho(i, j)
                    endif
                    
                endif
                fS = vS * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                rho_new(i, j) = rho(i, j) - dt * ((fE - fW) / dx + (fN - fS) / dy)
                rho_new(i,j) = max(rho_L, min(rho_H, rho_new(i, j)))
            enddo
        enddo
        !$omp end parallel do

        ! For ghost cells (zero-gradient)
        !$omp parallel do
        do j = 0, ny+1
            rho_new(0, j)= rho_new(1, j)
            rho_new(nx+1, j) = rho_new(nx, j)
        enddo
        !$omp end parallel do
        !$omp parallel do
        do i = 0, nx+1
            rho_new(i, 0) = rho_new(i, 1)
            rho_new(i, ny+1) = rho_new(i, ny)
        enddo
        !$omp end parallel do
    end subroutine advect_rho

    
    ! This subroutine aims at predicting the intermediate velocity fields u* and v*. (without pressure correction) => a.k.a SIMPLER algo.
    subroutine momentum_predictor
        integer :: i, j
        real(8) :: ufE, ufW, ufN, ufS, vfE, vfW, vfN, vfS
        real(8) :: rhoE, rhoW, rhoN, rhoS, rhof
        real(8) :: muE, muW, muN, muS
        real(8) :: phC, phD, phU
        real(8) :: cE, cW, cN, cS, conv
        real(8) :: diff
        real(8) :: ucell, vcell
        real(8) :: dudx_e, dudx_w, dudy_n, dudy_s
        real(8) :: dvdx_e, dvdx_w, dvdy_n, dvdy_s

        ! u-momentum predictor
        !$omp parallel do collapse(2) private(i, j, ufE, ufW, ufN, ufS, rhof, muE, muW, muN, muS, phC, phD, phU, &
        !$omp cE, cW, cN, cS, conv, diff, dudx_e, dudx_w, dudy_n, dudy_s)
        do j = 1, ny
            do i = 2, nx
                rhof = 0.5D0 * (rho(i-1, j) + rho(i, j))
                muE = mu_c(i, j)
                muW = mu_c(i-1, j)
                muN = 0.25D0 * (mu_c(i-1, j) + mu_c(i, j) + mu_c(i-1, j+1) + mu_c(i, j+1))
                muS = 0.25D0 * (mu_c(i-1, j) + mu_c(i, j) + mu_c(i-1, j-1) + mu_c(i, j-1))

                ! Advecting velocities at the four edges 
                ufE = 0.5D0 * (u(i, j) + u(i+1, j))
                ufW = 0.5D0 * (u(i, j) + u(i-1, j))
                ufN = 0.5D0 * (v(i-1, j+1)+v(i, j+1))
                ufS = 0.5D0 * (v(i-1, j)  +v(i, j))

                ! QUICK convective flux of u 
                ! East
                if (ufE >= 0.0D0) then
                    phC = u(i, j); 
                    phD = u(i+1, j)
                    
                    if (i >= 2) then
                        phU = u(i-1, j)
                    else
                        phU = u(i, j)
                    endif
                    
                else
                    phC = u(i+1, j) 
                    phD = u(i, j)
                    
                    if (i <= nx-1) then
                        phU = u(i+2, j)
                    else
                        phU = u(i+1, j)
                    endif
                    
                endif
                cE = ufE * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! West
                if (ufW >= 0.0D0) then
                    phC = u(i-1, j) 
                    phD = u(i, j)
                    
                    if (i >= 3) then
                        phU = u(i-2, j)
                    else
                        phU = u(i-1, j)
                    endif
                    
                else
                    phC = u(i, j)
                    phD = u(i-1, j)
                    
                    if (i <= nx) then
                        phU = u(i+1, j)
                    else
                        phU = u(i, j)
                    endif
                    
                endif
                cW = ufW * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! North
                if (ufN >= 0.0D0) then
                    phC = u(i, j) 
                    phD = u(i, j+1)
                    
                    if (j >= 2) then
                        phU = u(i, j-1)
                    else
                        phU = u(i, j)
                    endif
                    
                else
                    phC = u(i, j+1)
                    phD = u(i, j)
                    
                    if (j <= ny-1) then
                        phU = u(i,j+2)
                    else
                        phU = u(i,j+1)
                    endif
                    
                endif
                cN = ufN * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! South
                if (ufS >= 0.0D0) then
                    phC = u(i, j-1) 
                    phD = u(i, j)
                    
                    if (j >= 3) then
                        phU = u(i, j-2)
                    else
                        phU = u(i, j-1)
                    endif
                    
                else
                    phC = u(i, j)
                    phD = u(i, j-1)
                    
                    if (j <= ny) then
                        phU = u(i, j+1)
                    else
                        phU = u(i, j)
                    endif
                    
                endif
                cS = ufS * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                conv = (cE - cW) / dx + (cN - cS) / dy

                ! For diffusion
                dudx_e = (u(i+1, j) - u(i, j)) / dx
                dudx_w = (u(i, j) - u(i-1, j)) / dx
                dudy_n = (u(i, j+1) - u(i, j)) / dy
                dudy_s = (u(i, j) - u(i, j-1)) / dy
                diff = (muE * dudx_e - muW * dudx_w) / dx + (muN * dudy_n - muS * dudy_s) / dy

                us(i, j) = u(i, j) + dt * (-conv + diff / rhof)         ! Note that gravity only appears in v-momentum
            enddo
        enddo
        !$omp end parallel do

        ! v-momentum predictor
        !$omp parallel do collapse(2) private(i, j, vfE, vfW, vfN, vfS, rhof, muE, muW, muN, muS, phC, phD, phU, &
        !$omp cE, cW, cN, cS, conv, diff, dvdx_e, dvdx_w, dvdy_n, dvdy_s)
        do j = 2, ny
            do i = 1, nx
                rhof = 0.5D0 * (rho(i, j-1) + rho(i, j))
                
                muN  = mu_c(i, j)
                muS  = mu_c(i, j-1)
                muE  = 0.25D0 * (mu_c(i, j-1) + mu_c(i, j) + mu_c(i+1, j-1) + mu_c(i+1, j))
                muW  = 0.25D0 * (mu_c(i, j-1) + mu_c(i, j) + mu_c(i-1, j-1) + mu_c(i-1, j))
                
                vfE = 0.5D0 * (u(i+1, j-1) + u(i+1, j))
                vfW = 0.5D0 * (u(i, j-1) + u(i, j))
                vfN = 0.5D0 * (v(i, j) + v(i, j+1))
                vfS = 0.5D0 * (v(i, j) + v(i, j-1))
                
                ! East
                if (vfE >= 0.0D0) then
                    phC = v(i, j)
                    phD = v(i+1, j)
                    
                    if (i>=2) then
                        phU = v(i-1, j)
                    else
                        phU = v(i, j)
                    endif
                    
                else
                    phC = v(i+1,j)
                    phD = v(i,j)
                    
                    if (i<=nx-1) then
                        phU = v(i+2, j)
                    else
                        phU = v(i+1, j)
                    endif
                    
                endif
                cE = vfE * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! West
                if (vfW >= 0.0D0) then
                    phC = v(i-1, j)
                    phD = v(i, j)
                    
                    if (i >= 3) then
                        phU = v(i-2, j)
                    else
                        phU = v(i-1, j)
                    endif
                    
                else
                    phC = v(i, j)
                    phD = v(i-1, j)
                    
                    if (i <= nx) then
                        phU = v(i+1, j)
                    else
                        phU = v(i, j)
                    endif
                    
                endif
                cW = vfW * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! North
                if (vfN >= 0.0D0) then
                    phC = v(i, j)
                    phD = v(i, j+1)
                    
                    if (j >= 2) then
                        phU = v(i, j-1)
                    else
                        phU = v(i, j)
                    endif
                    
                else
                    phC = v(i, j+1)
                    phD = v(i, j)
                    
                    if (j <= ny-1) then
                        phU = v(i, j+2)
                    else
                        phU = v(i, j+1)
                    endif
                    
                endif
                cN = vfN * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                ! South
                if (vfS >= 0.0D0) then
                    phC = v(i, j-1)
                    phD = v(i, j)
                    
                    if (j >= 3) then
                        phU = v(i, j-2)
                    else
                        phU = v(i, j-1)
                    endif
                    
                else
                    phC = v(i, j)
                    phD = v(i, j-1)
                    
                    if (j <= ny) then
                        phU = v(i, j+1)
                    else
                        phU = v(i, j)
                    endif
                    
                endif
                cS = vfS * (0.375D0 * phD + 0.75D0 * phC - 0.125D0 * phU)

                conv = (cE - cW) / dx + (cN - cS) / dy

                dvdx_e = (v(i+1, j) - v(i, j)) / dx
                dvdx_w = (v(i, j) - v(i-1, j)) / dx
                dvdy_n = (v(i, j+1) - v(i, j)) / dy
                dvdy_s = (v(i, j) - v(i, j-1)) / dy
                diff = (muE * dvdx_e - muW * dvdx_w) / dx + (muN * dvdy_n - muS * dvdy_s) / dy

                vs(i, j) = v(i, j) + dt * (-conv + diff / rhof + g_y)
            enddo
        enddo
        !$omp end parallel do

        ! Applying BCs for u* and v*
        do j = 0, ny+1
            us(1, j) = 0.0D0
            us(nx+1, j) = 0.0D0
            us(0, j) = -us(2, j) 
            us(nx+2, j) = -us(nx, j)
            vs(0, j) = -vs(1, j)
            vs(nx+1, j) = -vs(nx, j)
        enddo
        do i = 0, nx+1
            vs(i,1) = 0.0D0
            vs(i,ny+1) = 0.0D0
            vs(i,0) = -vs(i,2)
            vs(i,ny+2) = -vs(i,ny)
            us(i,0) = -us(i,1)
            us(i,ny+1) = -us(i,ny)
        enddo
    end subroutine momentum_predictor


! This subroutine aims at solving the pressure Poisson equation using the red-black SOR method.
    subroutine pressure_solve
        integer :: i, j, it
        real(8) :: aE, aW, aN, aS, aP, bP
        real(8) :: rhoE, rhoW, rhoN, rhoS
        real(8) :: pnew, perr, pmean
        integer :: cnt

        do it = 1, p_max_it
            perr = 0.0D0
            ! Red cells
            !$omp parallel do collapse(2) private(i, j, aE, aW, aN, aS, aP, bP, rhoE, rhoW, rhoN, rhoS, pnew) reduction(max:perr)
            do j = 1, ny
                do i = 1, nx
                    if (mod(i+j, 2) /= 0) cycle     
                    rhoE = 0.5D0 * (rho_new(i, j) + rho_new(i+1, j))
                    rhoW = 0.5D0 * (rho_new(i, j) + rho_new(i-1, j))
                    rhoN = 0.5D0 * (rho_new(i, j) + rho_new(i, j+1))
                    rhoS = 0.5D0 * (rho_new(i, j) + rho_new(i, j-1))
                    
                    aE = 1.0D0 / (rhoE * dx * dx)
                    aW = 1.0D0 / (rhoW * dx * dx)
                    aN = 1.0D0 / (rhoN * dy * dy)
                    aS = 1.0D0 / (rhoS * dy * dy)
                    
                    if (i == 1)  aW = 0.0D0
                    if (i == nx) aE = 0.0D0
                    if (j == 1)  aS = 0.0D0
                    if (j == ny) aN = 0.0D0
                    aP = aE + aW + aN + aS
                    
                    if (aP <= 0.0D0) cycle

                    bP = ((us(i+1, j) - us(i, j)) / dx + (vs(i, j+1) - vs(i, j)) / dy) / dt
                    
                    pnew = (aE * p(i+1, j) + aW * p(i-1, j) + aN * p(i, j+1) + aS * p(i, j-1) - bP) / aP
                    pnew = (1.0D0 - omega_sor) * p(i, j) + omega_sor * pnew
                    
                    if (abs(pnew - p(i, j)) > perr) perr = abs(pnew - p(i, j))
                    p(i, j) = pnew
                    
                enddo
            enddo
            !$omp end parallel do
            
            ! Black cells
            !$omp parallel do collapse(2) private(i, j, aE, aW, aN, aS, aP, bP, rhoE, rhoW, rhoN, rhoS, pnew) reduction(max:perr)
            do j = 1, ny
                do i = 1, nx
                    if (mod(i+j, 2) == 0) cycle     
                    rhoE = 0.5D0 * (rho_new(i, j) + rho_new(i+1, j))
                    rhoW = 0.5D0 * (rho_new(i, j) + rho_new(i-1, j))
                    rhoN = 0.5D0 * (rho_new(i, j) + rho_new(i, j+1))
                    rhoS = 0.5D0 * (rho_new(i, j) + rho_new(i, j-1))

                    aE = 1.0D0 / (rhoE * dx * dx)
                    aW = 1.0D0 / (rhoW * dx * dx)
                    aN = 1.0D0 / (rhoN * dy * dy) 
                    aS = 1.0D0 / (rhoS * dy * dy)
                    
                    if (i == 1)  aW = 0.0D0
                    if (i == nx) aE = 0.0D0
                    if (j == 1)  aS = 0.0D0
                    if (j == ny) aN = 0.0D0
                    aP = aE + aW + aN + aS
                    
                    if (aP <= 0.0D0) cycle

                    bP = ((us(i+1, j) - us(i, j)) / dx + (vs(i, j+1) - vs(i, j)) / dy) / dt
                    
                    pnew = (aE * p(i+1, j) + aW * p(i-1, j) + aN * p(i, j+1) + aS * p(i, j-1) - bP) / aP
                    pnew = (1.0D0 - omega_sor) * p(i, j) + omega_sor * pnew
                    
                    if (abs(pnew - p(i, j)) > perr) perr = abs(pnew - p(i, j))
                    p(i, j) = pnew
                    
                end do
            end do
            !$omp end parallel do

            ! Applying Neumann BC on p
            !$omp parallel do
            do j = 0, ny+1
                p(0, j) = p(1, j)
                p(nx+1, j) = p(nx, j)
            enddo
            !$omp end parallel do
            
            !$omp parallel do
            do i = 0, nx+1
                p(i, 0) = p(i, 1)
                p(i, ny+1) = p(i, ny)
            enddo
            !$omp end parallel do
            
            if (perr < ptol) exit
        enddo

        ! Subtract mean
        pmean = 0.0D0
        do j = 1, ny
            do i = 1, nx
                pmean = pmean + p(i, j)
            enddo
        enddo
        pmean = pmean / dble(nx * ny)
        do j = 0, ny+1
            do i = 0, nx+1
                p(i,j) = p(i, j) - pmean
            enddo
        enddo
    end subroutine pressure_solve


! This subroutine aims at correcting the intermediate velocity field.
    subroutine velocity_correct
        integer :: i, j
        real(8) :: rhof

        !$omp parallel do collapse(2) private(i, j, rhof)
        do j = 1, ny
            do i = 2, nx
                rhof = 0.5D0 * (rho_new(i-1, j) + rho_new(i, j))
                u(i, j) = us(i, j) - dt / rhof * (p(i, j) - p(i-1, j)) / dx
            enddo
        enddo
        !$omp end parallel do

        !$omp parallel do collapse(2) private(i, j, rhof)
        do j = 2, ny
            do i = 1, nx
                rhof = 0.5D0 * (rho_new(i, j-1) + rho_new(i, j))
                v(i, j) = vs(i, j) - dt / rhof * (p(i, j) - p(i, j-1)) / dy
            enddo
        enddo
        !$omp end parallel do

        ! Reapplying  BCs for corrected velocities
        !$omp parallel do
        do j = 0, ny+1
            u(1, j) = 0.0D0 
            u(nx+1, j) = 0.0D0
            u(0, j) = -u(2, j)
            u(nx+2, j) = -u(nx, j)
            v(0, j) = -v(1, j) 
            v(nx+1, j) = -v(nx, j)
        enddo
        !$omp end parallel do
        
        !$omp parallel do
        do i = 0, nx+1
            v(i, 1) = 0.0D0
            v(i, ny+1) = 0.0D0
            v(i, 0) = -v(i, 2)
            v(i, ny+2) = -v(i, ny)
            u(i, 0) = -u(i, 1)
            u(i, ny+1) = -u(i, ny)
        enddo
        !$omp end parallel do
    end subroutine velocity_correct


! This subroutine aims at calculating the maximum divergence and maximum velocity in the domain.
    subroutine max_div(dmax, umax, vmax)
        real(8), intent(out) :: dmax, umax, vmax
        integer :: i,j
        real(8) :: d
        
        dmax = 0.0D0 
        umax = 0.0D0 
        vmax = 0.0D0
        
        do j = 1, ny
            do i = 1, nx
                d = (u(i+1, j) - u(i, j)) / dx + (v(i, j+1) - v(i, j)) / dy
                if (abs(d)>dmax) dmax = abs(d)
            enddo
        enddo
        do j = 1, ny
            do i = 1, nx+1
                if (abs(u(i, j)) > umax) umax = abs(u(i, j))
            enddo
        enddo
        do j = 1, ny+1
            do i = 1, nx
                if (abs(v(i, j)) > vmax) vmax = abs(v(i, j))
            enddo
        enddo
    end subroutine max_div

    
! This subroutine aims at outputting VTK files.
    subroutine write_vtk(istep)
        integer, intent(in) :: istep
        integer :: i, j
        real(8) :: uc, vc
        character(len=80) :: fn
        write(fn, '(A,I6.6,A)') 'RT_', istep, '.vtk'
        open(unit=21, file=trim(fn), status='replace')
        write(21,'(A)') '# vtk DataFile Version 3.0'
        write(21,'(A)') 'Rayleigh-Taylor density field'
        write(21,'(A)') 'ASCII'
        write(21,'(A)') 'DATASET STRUCTURED_POINTS'
        write(21,'(A,I0,1X,I0,1X,I0)') 'DIMENSIONS ', nx, ny, 1
        write(21,'(A,F10.6,1X,F10.6,1X,F10.6)') 'ORIGIN ', 0.5D0 * dx, 0.5D0 * dy, 0.0D0
        write(21,'(A,F10.6,1X,F10.6,1X,F10.6)') 'SPACING ', dx, dy, 1.0D0
        write(21,'(A,I0)') 'POINT_DATA ', nx*ny
        write(21,'(A)') 'SCALARS Density float 1'
        write(21,'(A)') 'LOOKUP_TABLE default'
        do j = 1, ny
            do i = 1, nx
                write(21,*) rho(i, j)
            enddo
        enddo
        write(21,'(A)') 'SCALARS Pressure float 1'
        write(21,'(A)') 'LOOKUP_TABLE default'
        do j = 1, ny
            do i = 1, nx
                write(21,*) p(i, j)
            enddo
        enddo
        write(21,'(A)') 'VECTORS Velocity float'
        do j = 1, ny
            do i = 1, nx
                uc = 0.5D0 * (u(i, j) + u(i+1, j))
                vc = 0.5D0 * (v(i, j) + v(i, j+1))
                write(21,*) uc, vc, 0.0D0
            enddo
        enddo
        close(21)
    end subroutine write_vtk

end program RTI_Solver
