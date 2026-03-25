      program convection_diffusion
      implicit none
      
      ! Declaration of variables
      integer*4 n, maxiter, iter, i, j
      parameter(n = 800, maxiter = 10000)
      
      real*8 phi(n, n), phiold(n, n)
      real*8 tol, maxres, err
      real*8 u, v, rho, gamma, delta, D, F, C1, C2, term1, term2
      parameter(tol = 1.0D -3, u = 1, v = 1, rho = 1, gamma = 0.1)
      parameter(delta = 0.1)
      
      D = gamma / delta
      F = rho * u
      C1 = D / (2*F + 4*D)
      C2 = (D + F) / (2*F + 4*D)
      
      ! Initialize flow field
      do i = 1, n
          do j = 1, n
              phi(i, j) = 0.0D0
          enddo
      enddo
      
      ! Applying boundary condition
      do i = 1, n
          phi(n, i) = 100
      enddo
      
      ! Gauss-Seidel Iteration
      do iter = 1, maxiter
          maxres = 0.0D0
          
          ! Save the old flow field to calculate the residual
          do i = 1, n
              do j = 1, n
                  phiold(i, j) = phi(i, j)
              enddo
          enddo
          
          do i = 2, n-1
              do j = 2, n-1
                  term1 = C1 * ( phi(i, j+1) + phi(i-1, j) )
                  term2 = C2 * ( phi(i, j-1) + phi(i+1, j) )
                  phi(i,j) = term1 + term2
                  
                  err = abs(phi(i, j) - phiold(i,j))
                  if (err >= maxres) maxres = err
              enddo
          enddo
          
          write(*,*) 'Iteration ', iter, 'Max Residual = ', maxres
          
          if (maxres <= tol) exit
          
      enddo
      
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