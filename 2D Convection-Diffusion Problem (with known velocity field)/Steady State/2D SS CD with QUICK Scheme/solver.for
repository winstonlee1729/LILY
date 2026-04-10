      program CD_QUICK
      implicit none
      
      ! Declaration of variables
      integer*4 n, iter, maxiter, i, j
      parameter(n = 100, maxiter = 1000000)
      
      ! Variables for timing the run
      integer*8 count1, count2, count_rate
      real*8 time_elapsed
      
      real*8 phi(n, n), phiold(n, n), tmp(n, n)
      real*8 tol, maxres, err, rf
      real*8 u, v, rho, gamma, dx, dy, A
      
      real*8 De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, delta_F
      real*8 alf_w, alf_e, alf_s, alf_n
      real*8 aw_up, ae_up, as_up, an_up, ap_up
      real*8 aw_qk, ae_qk, as_qk, an_qk, ap_qk
      real*8 aww_qk, aee_qk, ass_qk, ann_qk
      parameter(tol = 1.0D - 3, dx = 0.1, dy = 0.1, A = 1)
      parameter(u = 0.001, v = 0.01, rho = 1, gamma = 0.1, rf = 0.6)
      
      ! Initialization of F and D
      Fw = rho * u * A
      Fe = rho * u * A
      Fs = rho * v * A
      Fn = rho * v * A
      
      Dw = (gamma / dx) * A
      De = (gamma / dx) * A
      Ds = (gamma / dy) * A
      Dn = (gamma / dy) * A
      
      delta_F = (Fe - Fw) + (Fn - Fs)
      
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
      
      ! Upwind Scheme coeff. for safety ring
      aw_up = max(Fw, (Dw+(Fw/2)), 0.0D0)
      ae_up = max(-Fe, (De-(Fe/2)), 0.0D0)
      as_up = max(Fs, (Ds+(Fs/2)), 0.0D0)
      an_up = max(-Fn, (Dn-(Fn/2)), 0.0D0)
      ap_up = aw_up + ae_up + as_up + an_up + delta_F
      
      ! QUICK coeff.
      aw_qk  = Dw + (6.0D0/8.0D0)*alf_w*Fw + (1.0D0/8.0D0)*alf_e*Fe + 
     + (3.0D0/8.0D0)*(1.0D0-alf_w)*Fw
      aww_qk = (-1.0D0/8.0D0)*alf_w*Fw
    
      ae_qk  = De - (3.0D0/8.0D0)*alf_e*Fe - (6.0D0/8.0D0)* + 
     + (1.0D0-alf_e)*Fe - (1.0D0/8.0D0)*(1.0D0-alf_w)*Fw
      aee_qk = (1.0D0/8.0D0)*(1.0D0-alf_e)*Fe

      as_qk  = Ds + (6.0D0/8.0D0)*alf_s*Fs + (1.0D0/8.0D0)*alf_n*Fn + 
     + (3.0D0/8.0D0)*(1.0D0-alf_s)*Fs
      ass_qk = (-1.0D0/8.0D0)*alf_s*Fs
    
      an_qk  = Dn - (3.0D0/8.0D0)*alf_n*Fn - (6.0D0/8.0D0)*
     + (1.0D0-alf_n)*Fn - (1.0D0/8.0D0)*(1.0D0-alf_s)*Fs
      ann_qk = (1.0D0/8.0D0)*(1.0D0-alf_n)*Fn

      ap_qk = aw_qk + ae_qk + as_qk + an_qk + aww_qk + aee_qk + 
     + ass_qk + ann_qk + delta_F
      
      ! Initialization of field
      phi = 0.0D0
      phiold = 0.0D0
      
      ! Applying boundary condtion
      do i = 1, n
          phi(n, i) = 100.0D0
          phiold(n, i) = 100.0D0
      enddo
      
      ! Start the clock
      call system_clock(count1, count_rate)

      ! Gauss-Seidel Iteration
      do iter = 1, maxiter
          maxres = 0.0D0
          
          do i = 2, n-1
              do j = 2, n-1
                  
                  if (i==2 .or. i==n-1 .or. j==2 .or. j==n-1) then
                      tmp(i, j) = (aw_up * phi(i, j-1) + ae_up *
     +                phi(i, j+1) + as_up * phi(i+1, j) + an_up * +
     +                 phi(i-1, j)) / ap_up
                      
                  else
                      tmp(i, j) = (aw_qk*phi(i, j-1) + 
     +                             ae_qk*phi(i, j+1) + 
     +                             as_qk*phi(i+1, j) + 
     +                             an_qk*phi(i-1, j) + 
     +                             aww_qk*phi(i, j-2) + 
     +                             aee_qk*phi(i, j+2) + 
     +                             ass_qk*phi(i+2, j) + 
     +                             ann_qk*phi(i-2, j)) / ap_qk
                  endif
                  
                  phi(i,j) = (1.0D0 - rf) * phiold(i,j) + rf * tmp(i,j)
                  
                  err = abs(phi(i, j) - phiold(i, j))
                  if (err > maxres) maxres = err
                  
              enddo
          enddo
          
          phiold = phi
          
          ! Specifying iterations and residual
          if (mod(iter, 100) == 0) then
              write(*,*) 'Iteration ', iter, '| Max Residual: ', maxres
          endif
          
          if (maxres <= tol) exit
          
      enddo
      
      ! Stop the clock and calculate time
      call system_clock(count2, count_rate)
      time_elapsed = real(count2 - count1) / real(count_rate)

      write(*,*) '------------------------------------'
      write(*,*) 'Total Iterations: ', iter
      write(*,*) 'Total Time (Seconds): ', time_elapsed
      write(*,*) '------------------------------------'

      ! Outputting solution
      open(unit=10, file='flow_field_QUICK.txt', status='replace')
      do i = 1, n
          write(10, '(1000(F12.5, 1X))') (phi(i, j), j = 1, n)
      enddo
      close(10)
      
      stop
      end
      
                      
      