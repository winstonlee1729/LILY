      program SIMPLER_Lid_Driven_Cavity_omp
      use omp_lib
      implicit none
      
      ! Declaration of variables
      integer*4 nx, ny, maxiter, i, j, k, iter
      
      real*8 tol
      real*8 start_time, end_time
      real*8 rho, mu, ulid, Lx, Ly
      real*8 alph_u, alph_v, dx, dy
      
      real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), v(:,:)
      real*8, allocatable :: vstar(:,:), vold(:,:), uhat(:,:), vhat(:,:)
      real*8, allocatable :: p(:,:), pprime(:,:), bprime(:,:), b_p(:,:)
      real*8, allocatable :: du(:,:), dv(:,:)
      
      real*8 ae, aw, an, as, ap, delta_F
      real*8 Fe, Fw, Fn, Fs
      real*8 De, Dw, Dn, Ds, bmax
      real*8 xc, yc, uc, vc, pc

      ! Assigning values for variables
      nx = 100
      ny = 100
      maxiter = 500000
      tol = 1.0D-6
      rho = 1.0D0
      mu = 0.01D0
      ulid = 5.0D0
      Lx = 1.0D0
      Ly = 1.0D0
      alph_u = 0.7D0
      alph_v = 0.7D0

      ! Allocation
      allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1))
      allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1))
      allocate(uold(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1))
      allocate(uhat(0:nx+1, 0:ny+1), vhat(0:nx+1, 0:ny+1))
      allocate(p(0:nx+1, 0:ny+1), pprime(0:nx+1, 0:ny+1))
      allocate(bprime(0:nx+1, 0:ny+1), b_p(0:nx+1, 0:ny+1))
      allocate(du(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))

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
      start_time = omp_get_wtime()
      
      ! SIMPLER algorithm starts here
      do iter = 1, maxiter
      
         ! Save velocities 
!$OMP PARALLEL DO PRIVATE(j)
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
!$OMP PARALLEL DO PRIVATE(j, Fe, Fw, Fn, Fs, delta_F, 
!$OMP+ aw, ae, as, an, ap)
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

!$OMP PARALLEL DO PRIVATE(j, Fe, Fw, Fn, Fs, delta_F, 
!$OMP+ aw, ae, as, an, ap)
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

         ! Solve discretized pressure equation source
!$OMP PARALLEL DO PRIVATE(j)
         do i = 1, nx
            do j = 1, ny
               b_p(i,j) = (rho*dy*uhat(i,j)) + (rho*dx*vhat(i,j))
     +                  - (rho*dy*uhat(i+1,j)) - (rho*dx*vhat(i,j+1))
            enddo
         enddo

         ! Solve pressure equation 
         do k = 1, 100
!$OMP PARALLEL DO PRIVATE(j, aw, ae, as, an, ap)
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

         ! Solve momentum predictor
!$OMP PARALLEL DO PRIVATE(j)
         do i = 2, nx
            do j = 1, ny
               ustar(i,j) = uhat(i,j) + du(i,j)*(p(i-1,j) - p(i,j))
            enddo
         enddo
         
!$OMP PARALLEL DO PRIVATE(j)
         do i = 1, nx
            do j = 2, ny
               vstar(i,j) = vhat(i,j) + dv(i,j)*(p(i,j-1) - p(i,j))
            enddo
         enddo

         ! Solve pressure correction equation source
         bmax = 0.0D0
!$OMP PARALLEL DO PRIVATE(j) REDUCTION(MAX:bmax)
         do i = 1, nx
            do j = 1, ny
               pprime(i,j) = 0.0D0
               bprime(i,j) = (rho*dy*ustar(i,j)) + (rho*dx*vstar(i,j))
     +                  - (rho*dy*ustar(i+1,j)) - (rho*dx*vstar(i,j+1))
               if (abs(bprime(i,j)) > bmax) bmax = abs(bprime(i,j))
            enddo
         enddo

         ! Solve Pressure Correction (Parallel)
         do k = 1, 100
!$OMP PARALLEL DO PRIVATE(j, aw, ae, as, an, ap)
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

         ! Correcting velocities (Parallel)
!$OMP PARALLEL DO PRIVATE(j, uc, vc)
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
               
            enddo
         enddo
         
         write(*,*) 'Iteration: ', iter, '|Max mass error: ', bmax
         if (bmax <= tol) exit
         
      enddo  ! SIMPLER algorithm stops here
      
      ! Stop the clock
      end_time = omp_get_wtime()
      
      write(*,*) '------------------------------------'
      write(*,*) 'Total iterations: ', iter
      write(*,*) 'Total time (s): ', end_time - start_time
      write(*,*) '------------------------------------'

      ! Outputting solution
      open(unit=10, file='cavity_SIMPLER_omp.txt', status='unknown')
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

      ! Deallocation
      deallocate(u, ustar, uold, v, vstar, vold, uhat, vhat)
      deallocate(p, pprime, bprime, b_p, du, dv)

      stop
      end