      program TVD_solver
      implicit none
      
      ! Declaration of varibles
      integer*4 n, iter, maxiter, i, j, k
      parameter(n = 200, maxiter = 100000)
      integer*4 count1, count2, count_rate
      
      real*8 time_elapsed
      real*8 phi(n, n), phiold(n, n)
      real*8 tol, maxres, err
      real*8 u, v, rho, gamma, dx, dy, A
      real*8 Dw, De, Ds, Dn, Fw, Fe, Fs, Fn, delta_F
      real*8 alf_w, alf_e, alf_s, alf_n
      real*8 aw_up, ae_up, as_up, an_up, ap_up
      real*8 aw_t, ae_t, as_t, an_t, ap_t, Su_DC
      real*8 rw_p, re_p, rs_p, rn_p
      real*8 rw_m, re_m, rs_m, rn_m
      real*8 psi_w, psi_e, psi_s, psi_n, tmp
      real*8 term1, term2, term3, term4, C
      real*8 phi_new, rf
      
      parameter(tol = 1.0D-3, dx = 0.1, dy = 0.1, A = 1.0D0)
      parameter(u = 0.001, v = 0.01, rho = 1.0, gamma = 0.1)
      parameter(C = 1.0D-15, rf = 0.6D0)
      
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
      if (Fw > 0) then; alf_w = 1.0D0; else; alf_w = 0.0D0; endif
      if (Fe > 0) then; alf_e = 1.0D0; else; alf_e = 0.0D0; endif
      if (Fs > 0) then; alf_s = 1.0D0; else; alf_s = 0.0D0; endif
      if (Fn > 0) then; alf_n = 1.0D0; else; alf_n = 0.0D0; endif
      
      ! Upwind Scheme coeff. for safety ring
      aw_up = max(Fw, (Dw+(Fw/2.0D0)), 0.0D0)
      ae_up = max(-Fe, (De-(Fe/2.0D0)), 0.0D0)
      as_up = max(Fs, (Ds+(Fs/2.0D0)), 0.0D0)
      an_up = max(-Fn, (Dn-(Fn/2.0D0)), 0.0D0)
      ap_up = aw_up + ae_up + as_up + an_up + delta_F
      
      ! TVD coefficients 
      aw_t = Dw + max(Fw, 0.0D0)
      ae_t = De + max(-Fe, 0.0D0)
      as_t = Ds + max(Fs, 0.0D0)
      an_t = Dn + max(-Fn, 0.0D0)
      ap_t = aw_t + ae_t + as_t + an_t + delta_F
      
      ! Initialization of field
      phi = 0.0D0
      phiold = 0.0D0
      
      ! Applying boundary condition
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
                      ! Safety Ring uses standard Upwind
                      phi_new = (aw_up * phi(i, j-1) + ae_up *
     +                phi(i, j+1) + as_up * phi(i+1, j) + an_up * +
     +                phi(i-1, j)) / ap_up
                      
                  else
                      ! Interior nodes use TVD
                      re_p=(phi(i,j)-phi(i,j-1))/(phi(i,j+1)-phi(i,j)+C)
                      rw_p=(phi(i,j-1)-phi(i,j-2))/(phi(i,j)-phi(i,j-1)+
     +                     C)
                      rs_p=(phi(i+1,j)-phi(i+2,j))/(phi(i,j)-phi(i+1,j)+
     +                     C)
                      rn_p=(phi(i,j)-phi(i+1,j))/(phi(i-1,j)-phi(i,j)+C)
                      
                      re_m=(phi(i,j+2)-phi(i,j+1))/(phi(i,j+1)-phi(i,j)+
     +                     C)
                      rw_m=(phi(i,j+1)-phi(i,j))/(phi(i,j)-phi(i,j-1)+C)
                      rs_m=(phi(i-1,j)-phi(i,j))/(phi(i,j)-phi(i+1,j)+C)
                      rn_m=(phi(i-2,j)-phi(i-1,j))/(phi(i-1,j)-phi(i,j)+
     +                     C)
                      
                      ! Flux limiter functions
                      if (Fw > 0.0D0) then
                          tmp = min(2.0D0*rw_p, (1.0D0+3.0D0*rw_p)/4.0D0
     +                          , (3.0D0+rw_p)/4.0D0, 2.0D0)
                          psi_w = max(0.0D0, tmp)
                      else
                          tmp = min(2.0D0*rw_m, (1.0D0+3.0D0*rw_m)/4.0D0
     +                     , (3.0D0+rw_m)/4.0D0, 2.0D0)
                          psi_w = max(0.0D0, tmp)
                      endif
                      
                      if (Fe > 0.0D0) then
                          tmp = min(2.0D0*re_p, (1.0D0+3.0D0*re_p)/4.0D0
     +                     , (3.0D0+re_p)/4.0D0, 2.0D0)
                          psi_e = max(0.0D0, tmp) 
                      else
                          tmp = min(2.0D0*re_m, (1.0D0+3.0D0*re_m)/4.0D0
     +                     , (3.0D0+re_m)/4.0D0, 2.0D0)
                          psi_e = max(0.0D0, tmp) 
                      endif
                      
                      if (Fs > 0.0D0) then
                          tmp = min(2.0D0*rs_p, (1.0D0+3.0D0*rs_p)/4.0D0
     +                     , (3.0D0+rs_p)/4.0D0, 2.0D0)
                          psi_s = max(0.0D0, tmp)
                      else
                          tmp = min(2.0D0*rs_m, (1.0D0+3.0D0*rs_m)/4.0D0
     +                     , (3.0D0+rs_m)/4.0D0, 2.0D0)
                          psi_s = max(0.0D0, tmp)
                      endif

                      if (Fn > 0.0D0) then
                          tmp = min(2.0D0*rn_p, (1.0D0+3.0D0*rn_p)/4.0D0
     +                     , (3.0D0+rn_p)/4.0D0, 2.0D0)
                          psi_n = max(0.0D0, tmp)
                      else
                          tmp = min(2.0D0*rn_m, (1.0D0+3.0D0*rn_m)/4.0D0
     +                     , (3.0D0+rn_m)/4.0D0, 2.0D0)
                          psi_n = max(0.0D0, tmp)
                      endif
                      
                      ! Deferred correction 
                      tmp = (1.0D0 - alf_e) * psi_e - alf_e * psi_e
                      term1 = 0.5D0 * Fe * (phi(i,j+1)-phi(i,j)) * tmp
                      
                      tmp = alf_w * psi_w - (1.0D0 - alf_w) * psi_w
                      term2 = 0.5D0 * Fw * (phi(i,j)-phi(i,j-1)) * tmp
                  
                      tmp = (1.0D0 - alf_n) * psi_n - alf_n * psi_n
                      term3 = 0.5D0 * Fn * (phi(i-1,j)-phi(i,j)) * tmp
                      
                      tmp = alf_s * psi_s - (1.0D0 - alf_s) * psi_s
                      term4 = 0.5D0 * Fs * (phi(i,j)-phi(i+1,j)) * tmp
                      
                      Su_DC = term1 + term2 + term3 + term4
                      
                      phi_new = (aw_t * phi(i, j-1) + ae_t * + 
     +                          phi(i, j+1) + as_t * phi(i+1, j) + 
     +                          an_t * phi(i-1, j) + Su_DC) / ap_t
                  endif
                  
                  ! Under-Relaxation 
                  phi(i,j) = (1.0D0 - rf) * phiold(i,j) + rf * phi_new
                  
                  ! Check Error
                  err = abs(phi(i,j) - phiold(i,j))
                  if (err > maxres) maxres = err
                  
              enddo
          enddo
          
          phiold = phi
          
          if (mod(iter, 100) == 0) then
              write(*,*) 'Iteration: ', iter, '| Max Residual: ', maxres
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
      open(unit=10, file='flow_field_TVD.txt', status='replace')
      do i = 1, n
          write(10, '(1000(F12.5, 1X))') (phi(i, j), j = 1, n)
      enddo
      close(10)
      
      stop
      end   
          
              