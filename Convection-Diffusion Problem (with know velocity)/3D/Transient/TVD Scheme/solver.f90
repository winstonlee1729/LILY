program Transient_CD_TVD
    use omp_lib
    implicit none
    
! Declaration of variables
    ! Grid size and iteration number
    integer*4 nx, ny, nz
    integer*4 i, j, k
    integer*4 time_step, max_step
    integer*4 iter, maxiter, iter_limit
    integer*4 output_freq
    parameter(nx = 100, ny = 100, nz = 100)
    parameter(maxiter = 5000, max_step = 1000)
    parameter(output_freq = 25)
    
    ! Program convergence parameter
    real*8 tol, maxres, err, rhs, lhs

    ! Flow field
    real*8, allocatable :: phi(:,:,:), phiold(:,:,:), temp(:,:,:)
    
    ! Flow parameters
    real*8 rho, u, v, w
    real*8 area_w, area_e, area_n, area_s, area_t, area_b
    real*8 ga_w, ga_e, ga_n, ga_s, ga_t, ga_b
    real*8 dx, dy, dz, Lx, Ly, Lz, volume, delta_t
    real*8 Fe, Fw, Fn, Fs, Ft, Fb, delta_F
    real*8 De, Dw, Dn, Ds, Db, Dt
    parameter(tol = 1.0D-7)
    parameter(rho = 1.0D0, u = 0.0D0, v = 0.1D0, w = 0.0D0)
    parameter(ga_w = 0.01D0, ga_e = 0.01D0, ga_n = 0.01D0, ga_s = 0.01D0, ga_b = 0.01D0, ga_t = 0.01D0)
    parameter(Lx = 1.0D0, Ly = 1.0D0, Lz = 1.0D0)
    parameter(delta_t = 0.05D0)
    
    ! Output string for VTK files
    character(len=50) :: vtk_filename
    
    ! Temporal discretization scheme
    real*8 theta
    ! 1.0 = Fully Implicit | 0.5 = Crank-Nicolson | 0.0 = Explcit
    parameter(theta = 1.0D0)
    
    ! Flow direction
    real*8 alf_w, alf_e, alf_n, alf_s, alf_t, alf_b
    
    ! Hybrid scheme coefficients
    real*8 aw_hb, ae_hb, as_hb, an_hb, at_hb, ab_hb
    real*8 ap_hb_s, ap_hb, ap_o, ap_prime_hb
    
    ! TVD scheme coefficients
    real*8 aw_t, ae_t, an_t, as_t, at_t, ab_t
    real*8 ap_t_s, ap_t, ap_prime_t
    real*8 Su_DC
    
    ! Nodal values for the solver
    real*8 w_node, e_node, s_node, n_node, t_node, b_node, p_past
    
    ! TVD ratios and limiters
    real*8 rw_p, re_p, rn_p, rs_p, rt_p, rb_p
    real*8 rw_m, re_m, rn_m, rs_m, rt_m, rb_m
    real*8 psi_w, psi_e, psi_n, psi_s, psi_t, psi_b, tmp
    real*8 term1, term2, term3, term4, term5, term6, C
    parameter(C = 1.0D-15)
    
! Defining variables
    dx = Lx / dble(nx)
    dy = Ly / dble(ny)
    dz = Lz / dble(nz)
    
    area_w = dy * dz
    area_e = dy * dz
    area_n = dx * dz
    area_s = dx * dz
    area_t = dx * dy
    area_b = dx * dy
    
    volume = dx * dy * dz
    
    Fw = rho * u * area_w
    Fe = rho * u * area_e
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
    
    ap_o = rho * volume / delta_t
    
    ! Flow direction 
    if (Fw > 0.0D0) then; alf_w = 1.0D0; else; alf_w = 0.0D0; endif
    if (Fe > 0.0D0) then; alf_e = 1.0D0; else; alf_e = 0.0D0; endif
    if (Fs > 0.0D0) then; alf_s = 1.0D0; else; alf_s = 0.0D0; endif
    if (Fn > 0.0D0) then; alf_n = 1.0D0; else; alf_n = 0.0D0; endif
    if (Fb > 0.0D0) then; alf_b = 1.0D0; else; alf_b = 0.0D0; endif
    if (Ft > 0.0D0) then; alf_t = 1.0D0; else; alf_t = 0.0D0; endif
    
    ! Hybrid scheme coefficients
    aw_hb = max(Fw, Dw + (Fw / 2.0D0), 0.0D0)
    ae_hb = max(-Fe, De - (Fe / 2.0D0), 0.0D0)
    as_hb = max(Fs, Ds + (Fs / 2.0D0), 0.0D0)
    an_hb = max(-Fn, Dn - (Fn / 2.0D0), 0.0D0)
    ab_hb = max(Fb, Db + (Fb / 2.0D0), 0.0D0)
    at_hb = max(-Ft, Dt - (Ft / 2.0D0), 0.0D0)

    ap_hb_s = aw_hb + ae_hb + as_hb + an_hb + ab_hb + at_hb + delta_F
    ap_hb = ap_o + ap_hb_s * theta
    ap_prime_hb = ap_o - ap_hb_s * (1.0D0 - theta)
    
    ! TVD scheme coefficients
    aw_t = Dw + max(Fw, 0.0D0)
    ae_t = De + max(-Fe, 0.0D0)
    as_t = Ds + max(Fs, 0.0D0)
    an_t = Dn + max(-Fn, 0.0D0)
    ab_t = Db + max(Fb, 0.0D0)
    at_t = Dt + max(-Ft, 0.0D0)
    
    ap_t_s = aw_t + ae_t + as_t + an_t + ab_t + at_t + delta_F
    ap_t = ap_o + ap_t_s * theta
    ap_prime_t = ap_o - ap_t_s * (1.0D0 - theta)
    
! Allocate memory 
    allocate(phi(nx, ny, nz), phiold(nx, ny, nz), temp(nx, ny, nz))
    
! Initialize flow field
    phi = 0.0D0
    phiold = 0.0D0
    temp = 0.0D0
    ! Applying boundary condition
    !$omp parallel do collapse(2)
    do i = 1, nx
        do k = 1, nz
            phi(i, 1, k) = 100.0D0
            phiold(i, 1, k) = 100.0D0
            temp(i, 1, k) = 100.0D0
        enddo
    enddo
    !$omp end parallel do
    
! Check for inner iteration number
    if (theta == 0.0D0) then
        iter_limit = 1
    else
        iter_limit = maxiter
    endif
    
! Start time step
    do time_step = 1, max_step
        
        ! Save old field
        !$omp parallel do collapse(3)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    phiold(i, j, k) = phi(i, j, k)
                enddo
            enddo
        enddo
        !$omp end parallel do
        
        ! Inner iteration loop
        do iter = 1, maxiter
            maxres = 0.0D0
            
            !$omp parallel do collapse(3) private(i, j, k, err, re_p, rw_p, rn_p, rs_p, rt_p, rb_p, &
            !$omp re_m, rw_m, rn_m, rs_m, rt_m, rb_m, tmp, psi_w, psi_e, psi_s, psi_n, psi_t, psi_b, &
            !$omp term1, term2, term3, term4, term5, term6, rhs, lhs, Su_DC, &
            !$omp w_node, e_node, s_node, n_node, b_node, t_node, p_past) reduction(max:maxres)
            do i = 2, nx-1
                do j = 2, ny-1
                    do k = 2, nz-1
                        
                        if (i==2 .or. i==nx-1 .or. j==2 .or. j==ny-1 .or. k==2 .or. k==nz-1) then
                            ! Hybrid scheme
                            w_node = aw_hb * (theta * phi(i-1, j, k) + (1.0D0 - theta) * phiold(i-1, j, k))
                            e_node = ae_hb * (theta * phi(i+1, j, k) + (1.0D0 - theta) * phiold(i+1, j, k))
                            s_node = as_hb * (theta * phi(i, j-1, k) + (1.0D0 - theta) * phiold(i, j-1, k))
                            n_node = an_hb * (theta * phi(i, j+1, k) + (1.0D0 - theta) * phiold(i, j+1, k))
                            b_node = ab_hb * (theta * phi(i, j, k-1) + (1.0D0 - theta) * phiold(i, j, k-1))
                            t_node = at_hb * (theta * phi(i, j, k+1) + (1.0D0 - theta) * phiold(i, j, k+1))
                            p_past = ap_prime_hb * phiold(i, j, k)
                            
                            rhs = w_node + e_node + s_node + n_node + b_node + t_node + p_past
                            lhs = ap_hb * phi(i, j, k)
                            temp(i, j, k) = rhs / ap_hb
                        else
                            ! TVD scheme
                            re_p = (phi(i, j, k) - phi(i-1, j, k)) / (phi(i+1, j, k) - phi(i, j, k) + C)
                            rw_p = (phi(i-1, j, k) - phi(i-2, j, k)) / (phi(i, j, k) - phi(i-1, j, k) + C)
                            rn_p = (phi(i, j, k) - phi(i, j-1, k)) / (phi(i, j+1, k) - phi(i, j, k) + C)
                            rs_p = (phi(i, j-1, k) - phi(i, j-2, k)) / (phi(i, j, k) - phi(i, j-1, k) + C)
                            rt_p = (phi(i, j, k) - phi(i, j, k-1)) / (phi(i, j, k+1) - phi(i, j, k) + C)
                            rb_p = (phi(i, j, k-1) - phi(i, j, k-2)) / (phi(i, j, k) - phi(i, j, k-1) + C)
                            
                            re_m = (phi(i+2, j, k) - phi(i+1, j, k)) / (phi(i+1, j, k) - phi(i, j, k) + C)
                            rw_m = (phi(i+1, j, k) - phi(i, j, k)) / (phi(i, j, k) - phi(i-1, j, k) + C)
                            rn_m = (phi(i, j+2, k) - phi(i, j+1, k)) / (phi(i, j+1, k) - phi(i, j, k) + C) 
                            rs_m = (phi(i, j+1, k) - phi(i, j, k)) / (phi(i, j, k) - phi(i, j-1, k) + C)
                            rt_m = (phi(i, j, k+2) - phi(i, j, k+1)) / (phi(i, j, k+1) - phi(i, j, k) + C)
                            rb_m = (phi(i, j, k+1) - phi(i, j, k)) / (phi(i, j, k) - phi(i, j, k-1) + C)
                            
                            ! Use UMIST model for limiter function
                            if (Fw > 0.0D0) then
                                tmp = min(2.0D0*rw_p, (1.0D0+3.0D0*rw_p)/4.0D0, (3.0D0+rw_p)/4.0D0, 2.0D0)
                                psi_w = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rw_m, (1.0D0+3.0D0*rw_m)/4.0D0, (3.0D0+rw_m)/4.0D0, 2.0D0)
                                psi_w = max(0.0D0, tmp)
                            endif
                            
                            if (Fe > 0.0D0) then
                                tmp = min(2.0D0*re_p, (1.0D0+3.0D0*re_p)/4.0D0, (3.0D0+re_p)/4.0D0, 2.0D0)
                                psi_e = max(0.0D0, tmp) 
                            else
                                tmp = min(2.0D0*re_m, (1.0D0+3.0D0*re_m)/4.0D0, (3.0D0+re_m)/4.0D0, 2.0D0)
                                psi_e = max(0.0D0, tmp) 
                            endif
                          
                            if (Fs > 0.0D0) then
                                tmp = min(2.0D0*rs_p, (1.0D0+3.0D0*rs_p)/4.0D0, (3.0D0+rs_p)/4.0D0, 2.0D0)
                                psi_s = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rs_m, (1.0D0+3.0D0*rs_m)/4.0D0, (3.0D0+rs_m)/4.0D0, 2.0D0)
                                psi_s = max(0.0D0, tmp)
                            endif

                            if (Fn > 0.0D0) then
                                tmp = min(2.0D0*rn_p, (1.0D0+3.0D0*rn_p)/4.0D0, (3.0D0+rn_p)/4.0D0, 2.0D0)
                                psi_n = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rn_m, (1.0D0+3.0D0*rn_m)/4.0D0, (3.0D0+rn_m)/4.0D0, 2.0D0)
                                psi_n = max(0.0D0, tmp)
                            endif
                            
                            if (Fb > 0.0D0) then
                                tmp = min(2.0D0*rb_p, (1.0D0+3.0D0*rb_p)/4.0D0, (3.0D0+rb_p)/4.0D0, 2.0D0)
                                psi_b = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rb_m, (1.0D0+3.0D0*rb_m)/4.0D0, (3.0D0+rb_m)/4.0D0, 2.0D0)
                                psi_b = max(0.0D0, tmp)
                            endif

                            if (Ft > 0.0D0) then
                                tmp = min(2.0D0*rt_p, (1.0D0+3.0D0*rt_p)/4.0D0, (3.0D0+rt_p)/4.0D0, 2.0D0)
                                psi_t = max(0.0D0, tmp)
                            else
                                tmp = min(2.0D0*rt_m, (1.0D0+3.0D0*rt_m)/4.0D0, (3.0D0+rt_m)/4.0D0, 2.0D0)
                                psi_t = max(0.0D0, tmp)
                            endif
                            
                            tmp = (1.0D0 - alf_e) * psi_e - alf_e * psi_e
                            term1 = 0.5D0 * Fe * (phi(i+1, j, k) - phi(i, j, k)) * tmp
                          
                            tmp = alf_w * psi_w - (1.0D0 - alf_w) * psi_w
                            term2 = 0.5D0 * Fw * (phi(i, j, k) - phi(i-1, j, k)) * tmp
                      
                            tmp = (1.0D0 - alf_n) * psi_n - alf_n * psi_n
                            term3 = 0.5D0 * Fn * (phi(i, j+1, k) - phi(i, j, k)) * tmp
                          
                            tmp = alf_s * psi_s - (1.0D0 - alf_s) * psi_s
                            term4 = 0.5D0 * Fs * (phi(i, j, k) - phi(i, j-1, k)) * tmp
                            
                            tmp = (1.0D0 - alf_t) * psi_t - alf_t * psi_t
                            term5 = 0.5D0 * Ft * (phi(i, j, k+1) - phi(i, j, k)) * tmp
                          
                            tmp = alf_b * psi_b - (1.0D0 - alf_b) * psi_b
                            term6 = 0.5D0 * Fb * (phi(i, j, k-1) - phi(i, j, k)) * tmp
                            
                            Su_DC = term1 + term2 + term3 + term4 + term5 + term6
                            
                            w_node = aw_t * (theta * phi(i-1, j, k) + (1.0D0 - theta) * phiold(i-1, j, k))
                            e_node = ae_t * (theta * phi(i+1, j, k) + (1.0D0 - theta) * phiold(i+1, j, k))
                            s_node = as_t * (theta * phi(i, j-1, k) + (1.0D0 - theta) * phiold(i, j-1, k))
                            n_node = an_t * (theta * phi(i, j+1, k) + (1.0D0 - theta) * phiold(i, j+1, k))
                            b_node = ab_t * (theta * phi(i, j, k-1) + (1.0D0 - theta) * phiold(i, j, k-1))
                            t_node = at_t * (theta * phi(i, j, k+1) + (1.0D0 - theta) * phiold(i, j, k+1))
                            p_past = ap_prime_t * phiold(i, j, k)
                            
                            rhs = w_node + e_node + s_node + n_node + b_node + t_node + p_past + Su_DC
                            lhs = ap_t * phi(i, j, k)
                            temp(i, j, k) = rhs / ap_t
                        endif
                        
                        err = abs(rhs - lhs)
                        if (err > maxres) maxres = err
                    enddo
                enddo
            enddo
            !$omp end parallel do

            ! Update the main phi field with the newly calculated values
            !$omp parallel do collapse(3)
            do i = 2, nx-1
                do j = 2, ny-1
                    do k = 2, nz-1
                        phi(i, j, k) = temp(i, j, k)
                    enddo
                enddo
            enddo
            !$omp end parallel do
            
            if (theta > 0.0D0 .and. maxres <= tol) exit
        enddo
        
        if (theta == 0.0D0) then
            write(*,*) 'Time step:', time_step
        else
            write(*,*) 'Time step:', time_step, '|Iters:', iter, '|Res:', maxres
        endif
        
        ! Outputting solution 
        if (mod(time_step, output_freq) == 0) then
            write(vtk_filename, '(A, I4.4, A)') 'flow_', time_step, '.vtk'
            call write_vtk(vtk_filename, phi, nx, ny, nz, dx, dy, dz)
        endif
        
    enddo
    
! Deallocate memory
    deallocate(phi, phiold, temp)
    
    stop
end program Transient_CD_TVD


! This subroutine aims at outputting matrices as Legacy VTK format
subroutine write_vtk(filename, phi, nx, ny, nz, dx, dy, dz)
    implicit none
    character(len=*), intent(in) :: filename
    integer*4, intent(in) :: nx, ny, nz
    real*8, intent(in) :: dx, dy, dz
    real*8, intent(in) :: phi(nx, ny, nz)
    integer*4 :: i, j, k
    
    ! Open file
    open(unit=20, file=filename, status='replace')
    
    write(20, '(A)') '# vtk DataFile Version 3.0'
    write(20, '(A)') 'CFD Transient Output'
    write(20, '(A)') 'ASCII'
    write(20, '(A)') 'DATASET STRUCTURED_POINTS'
    
    ! Grid properties
    write(20, '(A, 3I8)') 'DIMENSIONS ', nx, ny, nz
    write(20, '(A, 3F10.6)') 'ORIGIN ', 0.0, 0.0, 0.0
    write(20, '(A, 3F10.6)') 'SPACING ', dx, dy, dz
    
    ! Data assignment 
    write(20, '(A, I12)') 'POINT_DATA ', nx * ny * nz
    write(20, '(A)') 'SCALARS phi float 1'
    write(20, '(A)') 'LOOKUP_TABLE default'
    
    ! Write phi values 
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                write(20, '(E16.6)') phi(i, j, k)
            enddo
        enddo
    enddo
    
    close(20)
end subroutine write_vtk