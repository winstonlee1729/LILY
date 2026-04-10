      program heat_diffusion
      implicit none
      
      ! Declaration of variables
      integer*4 n, maxiter, iter, i, j
      parameter(n = 100, maxiter = 10000)

      real*8 T(n, n), Told(n, n)
      real*8 tol, maxres, err
      parameter(tol = 1.0D -3)
      
      ! This section aims at initializing temperature field
      do i = 1, n
          do j = 1, n
              T(i, j) = 0.0D0
          enddo
      enddo
      
      !Applying boundary condition
      do i = 1, n
          T(1, i) = 100
      enddo
      
      ! Gauss-Seidel Iteration
      do iter = 1, maxiter
          maxres = 0.0D0
          
          ! Save the old temperature field to calculate the residual
          do i = 1, n
              do j = 1, n
                  Told(i, j) = T(i, j)
              enddo
          enddo
          
          do i = 2, n-1
              do j = 2, n-1
                  T(i, j) =0.25D0*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1))
                  
                  err = abs ( T(i, j) - Told(i, j))
                  if (err >= maxres) maxres = err
                  
              enddo
          enddo
          
          write(*,*) 'Iteration: ', iter, '|Max Residual = ', maxres
          
          if (maxres <= tol) exit
          
      enddo

      ! Open the output file
      open(unit=10, file='temperature_field.txt', status='unknown')
      
      ! Write data from row 1 to n using a strict format
      ! (1000(F12.5, 1X)) tells Fortran: "Print up to 1000 numbers on a SINGLE line. 
      ! Format each as a float with 5 decimal places, followed by 1 space."
      do i = 1, n
          write(10, '(1000(F12.5, 1X))') (T(i, j), j = 1, n)
      enddo
      
      ! Close the file
      close(10)
      
      
      stop
      end