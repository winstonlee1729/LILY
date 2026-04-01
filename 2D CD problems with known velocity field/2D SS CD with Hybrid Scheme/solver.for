      program CD_Hybrid
      implicit none
      
      ! Declaration of variables
      integer*4 n, iter, maxiter, i, j
      parameter(n = 100, maxiter = 10000)
      integer*4 count1, count2, count_rate
      
      real*8 phi(n, n), phiold(n, n), A
      real*8 tol, maxres, err, time_elapsed
      real*8 u, v, rho, gamma, dx, dy
      real*8 De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, delta_F
      real*8 ap, aw, ae, an, as, C1, C2, C3, C4
      real*8 term1, term2
      parameter(tol = 1.0D-5, dx = 0.1, dy = 0.1, A = 1)
      parameter(u = 0, v = 0.01, rho = 1, gamma = 0.1)
      
      ! Intialization of variables

      Fe = rho * u * A
      Fw = rho * u * A
      Fn = rho * v * A
      Fs = rho * v * A
      
      De = (gamma / dx) * A
      Dw = (gamma / dx) * A
      Dn = (gamma / dy) * A
      Ds = (gamma / dy) * A
      
      delta_F = Fe - Fw + Fn - Fs
      
      aw = max(Fw, (Dw+(Fw/2)), 0.0)
      ae = max(-Fe, (De-(Fe/2)), 0.0)
      as = max(Fs, (Ds+(Fs/2)), 0.0)
      an = max(-Fn, (Dn-(Fn/2)), 0.0)
      ap = aw + ae + as + an + delta_F
      
      C1 = aw / ap
      C2 = ae / ap
      C3 = as / ap
      C4 = an / ap
      
      ! Initialization of flow field
      do i = 1, n
          do j = 1, n
              phi(i, j) = 0.0D0
          enddo
      enddo
      
      ! Applying boundary condition
      do i = 1, n
          phi(n, i) = 100.0D0
      enddo
      
      ! Start the clock
      call SYSTEM_CLOCK(count1, count_rate)
      
      ! Gauss-Seidel Iteration
      do iter = 1, maxiter
          maxres = 0.0D0
          
          ! Save old flow field to calculate the residual
          do i = 1, n
              do j = 1, n
                  phiold(i, j) = phi(i, j)
              enddo
          enddo
          
          do i = 2, n-1
              do j = 2, n-1
                  term1 = C1 * phi(i, j-1) + C2 * phi(i, j+1)
                  term2 = C3 * phi(i+1, j) + C4 * phi(i-1, j)
                  phi(i,j) = term1 + term2
                  
                  err = abs(phi(i,j) - phiold(i,j))
                  if (err >= maxres) maxres = err
              enddo
          enddo
          
          write(*,*) 'Iteration: ', iter, '|Max Residual: ', maxres
          
          if (maxres <= tol) exit
          
      enddo
      
      ! Stop the clock
      call SYSTEM_CLOCK(count2, count_rate)
      time_elapsed = real(count2 - count1) / real(count_rate)

      write(*,*) '------------------------------------'
      write(*,*) 'Total Iterations: ', iter
      write(*,*) 'Total Time (Seconds): ', time_elapsed
      write(*,*) '------------------------------------'
      
      ! Open the output file
      open(unit=10, file='flow_field.txt', status='unknown')
      
      ! Write data from row 1 to n using a strict format
      ! (1000(F12.5, 1X)) tells Fortran: "Print up to 1000 numbers on a SINGLE line. 
      ! Format each as a float with 5 decimal places, followed by 1 space."
      do i = 1, n
          write(10, '(1000(F12.5, 1X))') (phi(i, j), j = 1, n)
      enddo
      
      ! Close the file
      close(10)
      
      
      stop
      end    