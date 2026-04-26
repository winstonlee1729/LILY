program CD_QUICK
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny, nz
    integer*4 iter, maxiter
    integer*4 i, j, k
    parameter(nx = 100, ny = 100, nz = 100)
    parameter(maxiter = 100000)
    
    ! Calculate program time
    real*8 start_time, end_time
    
    ! Flow field
    real*8, allocatable :: phi(:,:,:), phiold(:,:,:)
    real*8 phi_new
    
    ! Parameters for convergence
    real*8 maxres, err
    real*8, parameter :: tol = 1.0D-5, rf = 0.6D0
    
    ! Flow parameters
    real*8 u, v, w, rho
    real*8 Lx, Ly, Lz, dx, dy, dz
    real*8 area_w, area_e, area_n, area_s, area_t, area_b
    real*8 ga_w, ga_e, ga_n, ga_s, ga_t, ga_b
    parameter(Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0)
    parameter(u = 0.0D0, v = 0.1D0, w = 0.0D0)
    parameter(rho = 1.0D0)
    parameter(ga_w = 1.0D0, ga_e = 1.0D0, ga_n = 1.0D0, ga_s = 1.0D0, ga_t = 1.0D0, ga_b = 1.0D0)
    
    ! Equation varaibles
    real*8 De, Dw, Dn, Ds, Dt, Db
    real*8 Fe, Fw, Fn, Fs, Ft, Fb, delta_F
    
    ! For outputting solution
    real*8 xc, yc, zc
    
    ! Equation coefficients
    real*8 alf_w, alf_e, alf_s, alf_n, alf_b, alf_t
    real*8 aw_hb, ae_hb, as_hb, an_hb, at_hb, ab_hb, ap_hb
    real*8 aw_qk, ae_qk, as_qk, an_qk, at_qk, ab_qk, ap_qk
    real*8 aww_qk, aee_qk, ass_qk, ann_qk, att_qk, abb_qk
    real*8 x_term, y_term, z_term
    
! Allocate memory
    allocate(phi(nx, ny, nz), phiold(nx, ny, nz))
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    dz = Lz / dble(nz)
    
    area_e = dy * dz
    area_w = dy * dz
    area_n = dx * dz
    area_s = dx * dz
    area_b = dx * dy
    area_t = dx * dy
    
    Fe = rho * u * area_e
    Fw = rho * u * area_w
    Fn = rho * v * area_n
    Fs = rho * v * area_s
    Ft = rho * w * area_t
    Fb = rho * w * area_b
    delta_F = Fe - Fw + Fn - Fs + Ft - Fb
    
    De = (ga_e / dx) * area_e
    Dw = (ga_w / dx) * area_w
    Dn = (ga_n / dy) * area_n
    Ds = (ga_s / dy) * area_s
    Dt = (ga_t / dz) * area_t
    Db = (ga_b / dz) * area_b
    
    ! For alpha value
    if (Fw > 0) then; alf_w = 1.0D0; else; alf_w = 0.0D0; endif;
    if (Fe > 0) then; alf_e = 1.0D0; else; alf_e = 0.0D0; endif;
    if (Fs > 0) then; alf_s = 1.0D0; else; alf_s = 0.0D0; endif; 
    if (Fn > 0) then; alf_n = 1.0D0; else; alf_n = 0.0D0; endif; 
    if (Ft > 0) then; alf_t = 1.0D0; else; alf_t = 0.0D0; endif; 
    if (Fb > 0) then; alf_b = 1.0D0; else; alf_b = 0.0D0; endif; 
    
    ! Hybrid scheme coefficients
    aw_hb = max(Fw, Dw + (Fw / 2.0D0), 0.0D0)
    ae_hb = max(-Fe, De - (Fe / 2.0D0), 0.0D0)
    as_hb = max(Fs, Ds + (Fs / 2.0D0), 0.0D0)
    an_hb = max(-Fn, Dn - (Fn / 2.0D0), 0.0D0)
    ab_hb = max(Fb, Db + (Fb / 2.0D0), 0.0D0)
    at_hb = max(-Ft, Dt - (Ft / 2.0D0), 0.0D0)
    ap_hb = aw_hb + ae_hb + as_hb + an_hb + at_hb + ab_hb + delta_F
    
    ! QUICK scheme coefficients
    aw_qk = Dw + (6.0D0/8.0D0)*alf_w*Fw + (1.0D0/8.0D0)*alf_e*Fe + (3.0D0/8.0D0)*(1.0D0-alf_w)*Fw
    aww_qk = (-1.0D0/8.0D0)*alf_w*Fw
    
    ae_qk = De - (3.0D0/8.0D0)*alf_e*Fe - (6.0D0/8.0D0)*(1.0D0-alf_e)*Fe - (1.0D0/8.0D0)*(1.0D0-alf_w)*Fw
    aee_qk = (1.0D0/8.0D0)*(1.0D0-alf_e)*Fe

    as_qk = Ds + (6.0D0/8.0D0)*alf_s*Fs + (1.0D0/8.0D0)*alf_n*Fn + (3.0D0/8.0D0)*(1.0D0-alf_s)*Fs
    ass_qk = (-1.0D0/8.0D0)*alf_s*Fs
    
    an_qk = Dn - (3.0D0/8.0D0)*alf_n*Fn - (6.0D0/8.0D0)*(1.0D0-alf_n)*Fn - (1.0D0/8.0D0)*(1.0D0-alf_s)*Fs
    ann_qk = (1.0D0/8.0D0)*(1.0D0-alf_n)*Fn
    
    
    ab_qk = Db + (6.0D0/8.0D0)*alf_b*Fb + (1.0D0/8.0D0)*alf_t*Ft + (3.0D0/8.0D0)*(1.0D0-alf_b)*Fb
    abb_qk = (-1.0D0/8.0D0)*alf_b*Fb
    
    at_qk = Dt - (3.0D0/8.0D0)*alf_t*Ft - (6.0D0/8.0D0)*(1.0D0-alf_t)*Ft - (1.0D0/8.0D0)*(1.0D0-alf_b)*Fb
    att_qk = (1.0D0/8.0D0)*(1.0D0-alf_t)*Ft
    
    ap_qk = aw_qk + ae_qk + as_qk + an_qk + at_qk + ab_qk + aww_qk + aee_qk + ass_qk + ann_qk + att_qk + abb_qk + delta_F
    
! Initialization
    phi = 0.0D0
    phiold = 0.0D0
    ! Applying boundary condition
    !$omp parallel do collapse(2)
    do i = 1, nx
        do k = 1, nz
            phi(i, 1, k) = 100.0D0       
            phiold(i, 1, k) = 100.0D0
        enddo
    enddo
    !$omp end parallel do
    
! Start the clock
    start_time = omp_get_wtime()
    
! Jacobi Iteration
    do iter = 1, maxiter
        maxres = 0.0D0
        
        ! Save old flow field
        !$omp parallel do collapse(3)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    phiold(i, j, k) = phi(i, j, k)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        ! Pass 1 : Red cells
        !$omp parallel do collapse(3) private(i, j, k, err, x_term, y_term, z_term, phi_new) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (mod(i+j+k, 2) == 0) then
                        if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1 .or. k==2 .or. k==nz-1) then
                            x_term = aw_hb * phi(i-1, j, k) + ae_hb * phi(i+1, j, k)
                            y_term = as_hb * phi(i, j-1, k) + an_hb * phi(i, j+1, k)
                            z_term = ab_hb * phi(i, j, k-1) + at_hb * phi(i, j, k+1)
                            phi_new = (x_term + y_term + z_term) / ap_hb
                        else
                            x_term = aw_qk * phi(i-1, j, k) + ae_qk * phi(i+1, j, k) + aww_qk * phi(i-2, j, k) + aee_qk * phi(i+2, j, k)
                            y_term = as_qk * phi(i, j-1, k) + an_qk * phi(i, j+1, k) + ass_qk * phi(i, j-2, k) + ann_qk * phi(i, j+2, k)
                            z_term = ab_qk * phi(i, j, k-1) + at_qk * phi(i, j, k+1) + abb_qk * phi(i, j, k-2) + att_qk * phi(i, j, k+2)
                            phi_new = (x_term + y_term + z_term) / ap_qk
                        endif
                        
                        phi(i, j, k) = (1.0D0 - rf) * phiold(i, j, k) + rf * phi_new
                        
                        err = abs(phi(i, j, k) - phiold(i, j, k))
                        if(err >= maxres) maxres = err
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        ! Pass 2 : Black cells
        !$omp parallel do collapse(3) private(i, j, k, err, x_term, y_term, z_term, phi_new) reduction(max:maxres)
        do i = 2, nx-1
            do j = 2, ny-1
                do k = 2, nz-1
                    if (mod(i+j+k, 2) /= 0) then
                        if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1 .or. k==2 .or. k==nz-1) then
                            x_term = aw_hb * phi(i-1, j, k) + ae_hb * phi(i+1, j, k)
                            y_term = as_hb * phi(i, j-1, k) + an_hb * phi(i, j+1, k)
                            z_term = ab_hb * phi(i, j, k-1) + at_hb * phi(i, j, k+1)
                            phi_new = (x_term + y_term + z_term) / ap_hb
                        else
                            x_term = aw_qk * phi(i-1, j, k) + ae_qk * phi(i+1, j, k) + aww_qk * phi(i-2, j, k) + aee_qk * phi(i+2, j, k)
                            y_term = as_qk * phi(i, j-1, k) + an_qk * phi(i, j+1, k) + ass_qk * phi(i, j-2, k) + ann_qk * phi(i, j+2, k)
                            z_term = ab_qk * phi(i, j, k-1) + at_qk * phi(i, j, k+1) + abb_qk * phi(i, j, k-2) + att_qk * phi(i, j, k+2)
                            phi_new = (x_term + y_term + z_term) / ap_qk
                        endif
                        
                        phi(i, j, k) = (1.0D0 - rf) * phiold(i, j, k) + rf * phi_new
                        
                        err = abs(phi(i, j, k) - phiold(i, j, k))
                        if(err >= maxres) maxres = err
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        write(*,*) 'Iteration: ', iter, '| Max residual:', maxres

        
        if (maxres <= tol) exit
        
    enddo
    
! Stop the clock
    end_time = omp_get_wtime()
    write(*,*) '------------------------------------'
    write(*,*) 'Total iterations: ', iter
    write(*,*) 'Total time (s): ', end_time - start_time
    write(*,*) '------------------------------------'
    
! Open the output file
    open(unit=10, file='flow_field.txt', status='replace')
    
! Outputting solution
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                xc = (dble(i) - 0.5D0) * dx
                yc = (dble(j) - 0.5D0) * dy
                zc = (dble(k) - 0.5D0) * dz
                write(10, '(4E16.6)') xc, yc, zc, phi(i, j, k)
            enddo
        enddo
    enddo
    
! Close the output file
    close(10)
    
! Deallocate memory
    deallocate(phi, phiold)
    
    stop
end program CD_QUICK