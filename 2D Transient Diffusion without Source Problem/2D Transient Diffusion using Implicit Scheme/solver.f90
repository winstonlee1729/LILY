program Transient_Diffusion
    use omp_lib
    implicit none
    
    ! Declaration of variables
    integer*4 nx, ny, i, j
    integer*4 time_step, max_step, iter, maxiter
    integer*4 output_freq
    parameter(nx = 300, ny = 300)
    parameter(maxiter = 5000, max_step = 5000)
    parameter(output_freq = 25)
    
    real*8 maxres, err, tol, rhs, lhs
    parameter(tol = 1.0D-7)
    
    real*8, allocatable :: phi(:,:), phiold(:,:)
    real*8 ga_w, ga_e, ga_n, ga_s, area_w, area_e, area_n, area_s
    real*8 dx, dy, dt, Lx, Ly, ap, ae, aw, an, as, ap_o
    real*8 rho, Cp, volume, thickness
    real*8 w_node, e_node, n_node, s_node, p_past
    
    parameter(ga_e = 5.0D0, ga_w = 5.0D0, ga_n = 5.0D0, ga_s = 5.0D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, thickness = 1.0D0)
    parameter(dt = 0.05D0)
    parameter(rho = 1000.0D0, Cp = 1000.0D0)
    
    ! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    volume = dx * dy * thickness
    
    area_w = dy * thickness
    area_e = area_w
    area_n = dx * thickness
    area_s = area_n
    
    ap_o = rho * Cp * (volume / dt)
    aw = (ga_w * area_w) / dx
    ae = (ga_e * area_e) / dx
    an = (ga_n * area_n) / dy
    as = (ga_s * area_s) / dy
    ap = aw + ae + an + as + ap_o
    
    ! Allocate memory
    allocate(phi(nx, ny), phiold(nx, ny))
    
    ! Initialize flow field
    phi = 0.0D0
    phiold = 0.0D0
    
    ! Applying boundary condition
    !$omp parallel do
    do i = nx/2-3, nx/2+3
        phi(i, ny/2) = 10000.0D0
        phiold(i, ny/2) = 10000.0D0
    enddo
    !$omp end parallel do
    
    !$omp parallel do 
    do i = nx/2 -2, nx/2+2
        phi(i, ny/2+1) = 10000.0D0
        phiold(i, ny/2+1) = 10000.0D0
    enddo
    !$omp end parallel do
    
    !$omp parallel do 
    do i = nx/2 -2, nx/2+2
        phi(i, ny/2-1) = 10000.0D0
        phiold(i, ny/2-1) = 10000.0D0
    enddo
    !$omp end parallel do
    
    !$omp parallel do 
    do i = nx/2 -1, nx/2+1
        phi(i, ny/2+2) = 10000.0D0
        phiold(i, ny/2+2) = 10000.0D0
    enddo
    !$omp end parallel do
    
    !$omp parallel do 
    do i = nx/2 -1, nx/2+1
        phi(i, ny/2-2) = 10000.0D0
        phiold(i, ny/2-2) = 10000.0D0
    enddo
    !$omp end parallel do
    
    phi(nx/2, ny/2 + 3) = 10000.0D0
    phiold(nx/2, ny/2 + 3) = 10000.0D0
    
    phi(nx/2, ny/2 - 3) = 10000.0D0
    phiold(nx/2, ny/2 - 3) = 10000.0D0
    
    
    ! Open the file to write data
    open(unit=10, file='phi_out.txt', status='replace')
    
    !   Start time iteration
    do time_step = 1, max_step
        
        ! Save flow field
        !$omp parallel do collapse(2)
        do j = 1, ny
            do i = 1, nx
                phiold(i, j) = phi(i, j)
            enddo
        enddo
        !$omp end parallel do
        
        ! Inner iteration
        do iter = 1, maxiter
            maxres = 0.0D0
            
            !$omp parallel do private(i, e_node, w_node, n_node, s_node, p_past, rhs, lhs, err) reduction(max:maxres)
            do j = 2, ny-1
                do i = 2, nx-1
                    e_node = ae * phi(i+1, j)
                    w_node = aw * phi(i-1, j)
                    n_node = an * phi(i, j+1)
                    s_node = as * phi(i, j-1)
                    p_past = ap_o * phiold(i, j)
                    
                    
                    rhs = e_node + w_node + n_node + s_node + p_past
                    lhs = ap * phi(i, j)
                    
                    err = abs(rhs - lhs)
                    if (err >= maxres) maxres = err
                    
                    phi(i, j) = rhs / ap
                enddo
            enddo
            !$omp end parallel do
            
            if (maxres <= tol) exit
        enddo
        
        write(*,*) 'Time step:', time_step, '|Iters:', iter, '|Res:', maxres
        
        ! Outputting solution
        if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
            do j = 1, ny
                write(10, '(1000E15.6)') (phi(i, j), i = 1, nx)
            enddo
        endif
        
    enddo
    
    ! Close the file
    close(10)
    
    ! Deallocate memory
    deallocate(phi, phiold)
        
    stop 
    end program Transient_Diffusion
    
    