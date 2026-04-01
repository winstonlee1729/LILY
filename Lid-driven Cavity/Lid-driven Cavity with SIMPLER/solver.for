      program SIMPLER_Lid_Driven_Cavity
      implicit none
      
      ! Declaration of variables
      integer*4 nx, ny, maxiter, i, j, k, iter
      integer*4 count1, count2, count_rate
      parameter (nx = 100, ny = 100)
      parameter (maxiter = 500000)
      
      real*8 tol, time_elapsed
      parameter (tol = 1.0D-6)
      
      real*8 rho, mu, ulid, Lx, Ly
      real*8 alph_u, alph_v, dx, dy
      parameter (rho = 1.0D0, mu = 0.01D0, ulid = 5.0D0)
      parameter (Lx = 1.0D0, Ly = 1.0D0)
      parameter (alph_u = 0.7D0, alph_v = 0.7D0)
      
      real*8 u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1)
      real*8 v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1)
      real*8 uold(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1)
      real*8 uhat(0:nx+1, 0:ny+1), vhat(0:nx+1, 0:ny+1)
      real*8 p(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1)
      real*8 bprime(0:nx+1, 0:ny+1), b_p(0:nx+1, 0:ny+1)
      real*8 du(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1)
      
      real*8 ae, aw, an, as, ap, delta_F
      real*8 Fe, Fw, Fn, Fs
      real*8 De, Dw, Dn, Ds, bmax
      real*8 xc, yc, uc, vc, pc

      dx = Lx / dble(nx)
      dy = Ly / dble(ny)

      ! Initialization
      u = 0.0D0
      v = 0.0D0
      p = 0.0D0
      ustar = 0.0D0
      vstar = 0.0D0
      
      ! Applying boundary condition
      do i = 0, nx+1
         u(i, ny+1) = ulid
      enddo
      
      ! Start the clock
      call system_clock(count1, count_rate)
      
      ! SIMPLER algorithm starts here
      do iter = 1, maxiter
      
         ! Save velocities from previous iteration to old velocity field
         do i = 0, nx+1
            do j = 0, ny+1
               uold(i,j) = u(i,j)
               vold(i,j) = v(i,j)
            enddo
         enddo
         
         De = mu * (dy / dx)
         Dw = De
         Dn = mu * (dx / dy)
         Ds = Dn
         
         ! Calculate pseudo-velocities (u_hat, v_hat)
         do i = 2, nx
            do j = 1, ny
               Fe = 0.5D0 * rho * dy * (u(i,j) + u(i+1,j))
               Fw = 0.5D0 * rho * dy * (u(i-1,j) + u(i,j))
               Fn = 0.5D0 * rho * dx * (v(i-1,j+1) + v(i,j+1))
               Fs = 0.5D0 * rho * dx * (v(i-1,j) + v(i,j))
               delta_F = Fe - Fw + Fn - Fs
               
               aw = Dw + max(Fw, 0.0D0)
               ae = De + max(-Fe, 0.0D0)
               as = Ds + max(Fs, 0.0D0)
               an = Dn + max(-Fn, 0.0D0)
               ap = aw + ae + as + an + delta_F
               
               du(i,j) = dy / ap
               uhat(i,j) = (aw*u(i-1,j) + ae*u(i+1,j) + 
     +                      as*u(i,j-1) + an*u(i,j+1)) / ap
            enddo
         enddo

         do i = 1, nx
            do j = 2, ny
               Fe = 0.5D0 * rho * dy * (u(i+1,j-1) + u(i+1,j))
               Fw = 0.5D0 * rho * dy * (u(i,j-1) + u(i,j))
               Fn = 0.5D0 * rho * dx * (v(i,j) + v(i,j+1))
               Fs = 0.5D0 * rho * dx * (v(i,j-1) + v(i,j))
               delta_F = Fe - Fw + Fn - Fs
               
               aw = Dw + max(Fw, 0.0D0)
               ae = De + max(-Fe, 0.0D0)
               as = Ds + max(Fs, 0.0D0)
               an = Dn + max(-Fn, 0.0D0)
               ap = aw + ae + as + an + delta_F
               
               dv(i,j) = dx / ap
               vhat(i,j) = (aw*v(i-1,j) + ae*v(i+1,j) + 
     +                      as*v(i,j-1) + an*v(i,j+1)) / ap
            enddo
         enddo

         ! Solve discretized pressure equation from discretized continuity equation
         do i = 1, nx
            do j = 1, ny
               b_p(i,j) = (rho*dy*uhat(i,j)) + (rho*dx*vhat(i,j))
     +                  - (rho*dy*uhat(i+1,j)) - (rho*dx*vhat(i,j+1))
            enddo
         enddo

         do k = 1, 100
            do i = 1, nx
               do j = 1, ny
                  aw = rho * dy * du(i,j)
                  ae = rho * dy * du(i+1,j)
                  as = rho * dx * dv(i,j)
                  an = rho * dx * dv(i,j+1)
                  
                  if (i == 1) aw = 0.0D0
                  if (i == nx) ae = 0.0D0
                  if (j == 1) as = 0.0D0
                  if (j == ny) an = 0.0D0
                  ap = aw + ae + as + an
                  
                  p(i,j) = (ae*p(i+1,j) + aw*p(i-1,j) + 
     +                      an*p(i,j+1) + as*p(i,j-1) + b_p(i,j)) / ap
               enddo
            enddo
            
            p(1,1) = 0.0D0 ! Pin pressure
            
         enddo

         ! Solve momentum equation for u* and v* using the newly obtained pressure field
         do i = 2, nx
            do j = 1, ny
               ustar(i,j) = uhat(i,j) + du(i,j)*(p(i-1,j) - p(i,j))
            enddo
         enddo
         
         do i = 1, nx
            do j = 2, ny
               vstar(i,j) = vhat(i,j) + dv(i,j)*(p(i,j-1) - p(i,j))
            enddo
         enddo

         ! Solve pressure correction equation 
         bmax = 0.0D0
         do i = 1, nx
            do j = 1, ny
               pprime(i,j) = 0.0D0
               
               bprime(i,j) = (rho*dy*ustar(i,j)) + (rho*dx*vstar(i,j))
     +                - (rho*dy*ustar(i+1,j)) - (rho*dx*vstar(i,j+1))
               
               if (abs(bprime(i,j)) > bmax) then
                   bmax = abs(bprime(i,j))
               endif
               
            enddo
         enddo

         do k = 1, 100
            do i = 1, nx
               do j = 1, ny
                  aw = rho * dy * du(i,j)
                  ae = rho * dy * du(i+1,j)
                  as = rho * dx * dv(i,j)
                  an = rho * dx * dv(i,j+1)
                  
                  if (i == 1) aw = 0.0D0
                  if (i == nx) ae = 0.0D0
                  if (j == 1) as = 0.0D0
                  if (j == ny) an = 0.0D0
                  ap = aw + ae + as + an
                  
                  pprime(i,j) = (ae*pprime(i+1,j) + aw*pprime(i-1,j) + 
     +            an*pprime(i,j+1) + as*pprime(i,j-1) + bprime(i,j))/ap
               enddo
            enddo
            
            pprime(1,1) = 0.0D0 ! Pin pressure
            
         enddo

         ! Correcting velocities
         do i = 1, nx
            do j = 1, ny
                
               if (i > 1) then
                  uc=ustar(i,j)+du(i,j)*(pprime(i-1,j) - pprime(i,j))
                  u(i,j) = alph_u * uc + (1.0D0 - alph_u) * uold(i,j)
               endif
               
               if (j > 1) then
                  vc=vstar(i,j)+dv(i,j)*(pprime(i,j-1) - pprime(i,j))
                  v(i,j) = alph_v * vc + (1.0D0 - alph_v) * vold(i,j)
               endif
               
            end do
         end do
         
         !if (mod(iter, 100) == 0) then
            write(*,*) 'Iteration: ', iter, '|Max mass error: ', bmax
         !end if
            
         if (bmax <= tol) exit
         
      end do  ! SIMPLER algorithm stops here
      

      ! Stop the clock
      call system_clock(count2, count_rate)
      time_elapsed = real(count2 - count1) / real(count_rate)

      write(*,*) '------------------------------------'
      write(*,*) 'SIMPLER Converged in: ', iter, ' iterations'
      write(*,*) 'Total Time: ', time_elapsed, ' seconds'
      write(*,*) '------------------------------------'

      ! Outputting solution
      open(unit=10, file='cavity_SIMPLER.txt', status='unknown')
      do j = 1, ny
         do i = 1, nx
            xc = (dble(i) - 0.5D0) * dx
            yc = (dble(j) - 0.5D0) * dy
            uc = 0.5D0 * (u(i,j) + u(i+1,j))
            vc = 0.5D0 * (v(i,j) + v(i,j+1)) 
            pc = p(i,j)
            write(10, '(5E16.6)') xc, yc, uc, vc, pc
         enddo
      enddo
      close(10)

      stop
      end    