      program CD_QUICK_omp
      implicit none
      
      ! Declaration of variables
      integer*4 n, iter, maxiter, i, j
      parameter(maxiter = 1000000)
      
      real*8, allocatable :: phi(:, :), phiold(:, :)
      real*8 tol, maxres, err, rf
      real*8 u, v, rho, gamma, dx, dy, A
      real*8 phi_new 
      
      real*8 De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, delta_F
      real*8 alf_w, alf_e, alf_s, alf_n
      real*8 aw_up, ae_up, as_up, an_up, ap_up
      real*8 aw_qk, ae_qk, as_qk, an_qk, ap_qk
      real*8 aww_qk, aee_qk, ass_qk, ann_qk
      real*8 start_time, end_time, omp_get_wtime
      
      parameter(tol = 1.0D-3, dx = 0.1, dy = 0.1, A = 1)
      parameter(u = 0.001, v = 0.01, rho = 1, gamma = 0.1, rf = 0.6)
      
      ! Grid allocation
      n = 100
      allocate(phi(n,n), phiold(n,n))
      
      Fw = rho * u * A
      Fe = rho * u * A
      Fs = rho * v * A
      Fn = rho * v * A
      
      Dw = (gamma / dx) * A
      De = (gamma / dx) * A
      Ds = (gamma / dy) * A
      Dn = (gamma / dy) * A
      
      delta_F = (Fe - Fw) + (Fn - Fs)
      
      ! Specifying flow direction 
      ! For alpha w
      if (Fw > 0) then
          alf_w = 1.0D0
      else
          alf_w = 0.0D0
      endif
      ! For alpha e
      if (Fe > 0) then
          alf_e = 1.0D0
      else
          alf_e = 0.0D0
      endif
      ! For alpha s
      if (Fs > 0) then
          alf_s = 1.0D0
      else
          alf_s = 0.0D0
      endif
      ! For alpha n
      if (Fn > 0) then
          alf_n = 1.0D0
      else
          alf_n = 0.0D0
      endif
      
      ! Upwind Scheme Coeff. (Safety Ring)
      aw_up = max(Fw, (Dw+(Fw/2.0D0)), 0.0D0)
      ae_up = max(-Fe, (De-(Fe/2.0D0)), 0.0D0)
      as_up = max(Fs, (Ds+(Fs/2.0D0)), 0.0D0)
      an_up = max(-Fn, (Dn-(Fn/2.0D0)), 0.0D0)
      ap_up = aw_up + ae_up + as_up + an_up + delta_F
      
      ! QUICK Coeff. 
      aw_qk = Dw + (6.0D0/8.0D0)*alf_w*Fw + (1.0D0/8.0D0)*alf_e*Fe
      aw_qk = aw_qk + (3.0D0/8.0D0)*(1.0D0-alf_w)*Fw
      aww_qk = (-1.0D0/8.0D0)*alf_w*Fw
    
      ae_qk = De - (3.0D0/8.0D0)*alf_e*Fe
      ae_qk = ae_qk - (6.0D0/8.0D0)*(1.0D0-alf_e)*Fe 
      ae_qk = ae_qk - (1.0D0/8.0D0)*(1.0D0-alf_w)*Fw
      aee_qk = (1.0D0/8.0D0)*(1.0D0-alf_e)*Fe

      as_qk = Ds + (6.0D0/8.0D0)*alf_s*Fs + (1.0D0/8.0D0)*alf_n*Fn
      as_qk = as_qk + (3.0D0/8.0D0)*(1.0D0-alf_s)*Fs
      ass_qk = (-1.0D0/8.0D0)*alf_s*Fs
    
      an_qk = Dn - (3.0D0/8.0D0)*alf_n*Fn 
      an_qk = an_qk - (6.0D0/8.0D0)*(1.0D0-alf_n)*Fn 
      an_qk = an_qk - (1.0D0/8.0D0)*(1.0D0-alf_s)*Fs
      ann_qk = (1.0D0/8.0D0)*(1.0D0-alf_n)*Fn

      ap_qk = aw_qk + ae_qk + as_qk + an_qk + aww_qk + aee_qk
      ap_qk = ap_qk + ass_qk + ann_qk + delta_F
      
      ! Initialization of field
      phi = 0.0D0
      phiold = 0.0D0
      
      ! Applying Boundary condition
      do i = 1, n
          phi(n, i) = 100.0D0
          phiold(n, i) = 100.0D0
      enddo
      

      ! Start the clock
      start_time = omp_get_wtime()
      
      do iter = 1, maxiter
          maxres = 0.0D0
          
C         Pass 1: Red cells
C$OMP PARALLEL DO PRIVATE(i, j, err, phi_new) REDUCTION(MAX:maxres)
          do i = 2, n-1
              do j = 2, n-1
                  if (mod(i+j, 2) == 0) then
                      
                      if (i==2 .or. i==n-1 .or. j==2 .or. j==n-1) then
                          phi_new = (aw_up*phi(i, j-1) + 
     +                    ae_up*phi(i, j+1) + as_up*phi(i+1, j) + 
     +                    an_up*phi(i-1, j)) / ap_up
                      else
                          phi_new = (aw_qk*phi(i, j-1) + 
     +                    ae_qk*phi(i, j+1) + as_qk*phi(i+1, j) + 
     +                    an_qk*phi(i-1, j) + aww_qk*phi(i, j-2) + 
     +                    aee_qk*phi(i, j+2) + ass_qk*phi(i+2, j) + 
     +                    ann_qk*phi(i-2, j)) / ap_qk
                      endif
                      
                      phi(i,j) = (1.0D0-rf)*phiold(i,j) + rf*phi_new
                      
                      err = abs(phi(i, j) - phiold(i, j))
                      if (err > maxres) maxres = err
                      
                  endif
              enddo
          enddo
C$OMP END PARALLEL DO

C         Pass 2: Black cells
C$OMP PARALLEL DO PRIVATE(i, j, err, phi_new) REDUCTION(MAX:maxres)
          do i = 2, n-1
              do j = 2, n-1
                  if (mod(i+j, 2) /= 0) then
                      
                      if (i==2 .or. i==n-1 .or. j==2 .or. j==n-1) then
                          phi_new = (aw_up*phi(i, j-1) + 
     +                    ae_up*phi(i, j+1) + as_up*phi(i+1, j) + 
     +                    an_up*phi(i-1, j)) / ap_up
                      else
                          phi_new = (aw_qk*phi(i, j-1) + 
     +                    ae_qk*phi(i, j+1) + as_qk*phi(i+1, j) + 
     +                    an_qk*phi(i-1, j) + aww_qk*phi(i, j-2) + 
     +                    aee_qk*phi(i, j+2) + ass_qk*phi(i+2, j) + 
     +                    ann_qk*phi(i-2, j)) / ap_qk
                      endif
                      
                      phi(i,j) = (1.0D0-rf)*phiold(i,j) + rf*phi_new
                      
                      err = abs(phi(i, j) - phiold(i, j))
                      if (err > maxres) maxres = err
                      
                  endif
              enddo
          enddo
C$OMP END PARALLEL DO
          
          phiold = phi
          
          write(*,*) 'Iterations: ', iter, '|Max residual: ', maxres
          
          if (maxres <= tol) exit
          
      enddo
      
      ! Stop the clock
      end_time = omp_get_wtime()
      
      write(*,*) '------------------------------------'
      write(*,*) 'Total iterations: ', iter
      write(*,*) 'Total time (s): ', end_time - start_time
      write(*,*) '------------------------------------'
      

      ! Outputting
      open(unit=10, file='flow_field_QUICK.txt', status='replace')
      do i = 1, n
          write(10, '(1000(F12.5, 1X))') (phi(i, j), j = 1, n)
      enddo
      close(10)
      
      deallocate(phi, phiold)
      
      stop
      end