      program SIMPLEC_Lid_Driven_Cavity_omp
      use omp_lib
      implicit none
      
      ! Declaration of variables
      integer*4 nx, ny, maxiter, i, j, k, iter
      
      real*8 tol
      real*8 start_time, end_time
      real*8 rho, mu, ulid, Lx, Ly
      real*8 alph_p, alph_u, alph_v, dx, dy
      
      real*8, allocatable :: u(:,:), ustar(:,:), uold(:,:), du(:,:)
      real*8, allocatable :: v(:,:), vstar(:,:), vold(:,:), dv(:,:)
      real*8, allocatable :: p(:,:), pstar(:,:), pprime(:,:),bprime(:,:)
      
      real*8 ae, aw, an, as, ap, delta_f
      real*8 Fe, Fw, Fn, Fs
      real*8 De, Dw, Dn, Ds, bmax, uc, vc
      real*8 xc, yc, pc

      ! Assigning values of variables
      nx = 100
      ny = 100
      maxiter = 500000
      tol = 1.0D-6
      
      rho = 1.0D0
      mu = 0.01D0
      ulid = 1.0D0
      Lx = 1.0D0
      Ly = 1.0D0
      
      ! For SIMPLEC, alph_p is assumed to be 1
      alph_p = 1.0D0
      alph_u = 0.5D0
      alph_v = 0.5D0

      ! Memory allocation
      allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1))
      allocate(uold(0:nx+1, 0:ny+1), du(0:nx+1, 0:ny+1))
      allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1))
      allocate(vold(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
      allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1))
      allocate(pprime(0:nx+1, 0:ny+1), bprime(0:nx+1, 0:ny+1))

      dx = Lx / dble(nx)
      dy = Ly / dble(ny)

      ! Initialization
      u = 0.0D0
      v = 0.0D0
      p = 0.0D0
      ustar = 0.0D0
      vstar = 0.0D0
      pstar = 0.0D0
      
      ! Applying boundary condition
      do i = 0, nx+1
         u(i, ny+1) = ulid
      enddo
      
      ! Start the clock
      start_time = omp_get_wtime()
      
      ! SIMPLEC algorithm
      do iter = 1, maxiter
         
         ! Save old velocities for explicit under-relaxation
!$omp parallel do private(j)
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
         
         ! u-momentum predictor
!$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, 
!$omp+ aw, ae, as, an, ap)
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
               
               du(i,j) = dy / (ap/alph_u - (aw + ae + as + an))
               
               ustar(i,j) = (aw * uold(i-1,j) + ae * uold(i+1,j) + 
     +                       as * uold(i,j-1) + an * uold(i,j+1) + 
     +                       (pstar(i-1,j) - pstar(i,j)) * dy) / ap
            enddo
         enddo
         
         ! v-momentum predictor
!$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, 
!$omp+ aw, ae, as, an, ap)
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
               
               dv(i,j) = dx / (ap/alph_v - (aw + ae + as + an))
               
               vstar(i,j) = (aw * vold(i-1,j) + ae * vold(i+1,j) + 
     +                       as * vold(i,j-1) + an * vold(i,j+1) + 
     +                       (pstar(i,j-1) - pstar(i,j)) * dx) / ap
            enddo
         enddo
         
         ! Reapplying boundary conditions
         do i = 0, nx+1
            ustar(i, 0)    = 0.0D0    
            ustar(i, ny+1) = ulid     
            vstar(i, 1)    = 0.0D0    
            vstar(i, ny+1) = 0.0D0    
         enddo
         
         do j = 0, ny+1
            ustar(1, j)    = 0.0D0    
            ustar(nx+1, j) = 0.0D0    
            vstar(0, j)    = 0.0D0
            vstar(nx+1, j) = 0.0D0
         enddo
         
         ! Calculate mass source and reduction for bmax
         bmax = 0.0D0
!$omp parallel do private(j) reduction(max:bmax)
         do i = 1, nx
            do j = 1, ny
               pprime(i,j) = 0.0D0
               bprime(i,j) = (rho * dy * ustar(i,j)) +
     +                       (rho * dx * vstar(i,j)) -
     +                       (rho * dy * ustar(i+1,j)) -
     +                       (rho * dx * vstar(i,j+1))
               
               if (abs(bprime(i,j)) > bmax) then
                   bmax = abs(bprime(i,j))
               endif
               
            enddo
         enddo
         
         ! Gauss-Seidel Iteration for pressure correction
         do k = 1, 100
!$omp parallel do private(j, aw, ae, as, an, ap)
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
                  ap = aw + ae + an + as
                  
                  pprime(i,j) = (ae * pprime(i+1,j) + 
     +                           aw * pprime(i-1,j) + 
     +                           an * pprime(i,j+1) + 
     +                           as * pprime(i,j-1) + 
     +                           bprime(i,j)) / ap
               enddo
            enddo
            
            pprime(1,1) = 0.0D0   ! Pin pressure
            
         enddo
         
         ! Field updates
!$omp parallel do private(j, uc, vc)
         do i = 1, nx
            do j = 1, ny
               p(i,j) = pstar(i,j) + alph_p * pprime(i,j)
               pstar(i,j) = p(i,j)
               
               if (i > 1) then
                  uc = ustar(i,j)+du(i,j)*(pprime(i-1,j)-pprime(i,j))
                  u(i,j) = alph_u * uc + (1.0D0 - alph_u) * uold(i,j)
               endif
               
               if (j > 1) then
                  vc = vstar(i,j)+dv(i,j)*(pprime(i,j-1)-pprime(i,j))
                  v(i,j) = alph_v * vc + (1.0D0 - alph_v) * vold(i,j)
               endif
            enddo
         enddo
         
         !if (mod(iter, 100) == 0) then
            write(*,*) 'Iteration: ', iter, '|Max mass error: ', bmax
         !endif
         
         if (bmax <= tol) exit
         
      enddo  ! SIMPLEC algorithm stops here
      
      ! Stop the clock
      end_time = omp_get_wtime()
      
      write(*,*) '------------------------------------'
      write(*,*) 'Total iterations: ', iter
      write(*,*) 'Total time (s): ', end_time - start_time
      write(*,*) '------------------------------------'

      ! Outputting solution
      open(unit=10, file='cavity_SIMPLEC_omp.txt', status='unknown')
      do j = 1, ny
         do i = 1, nx
            xc = (dble(i) - 0.5D0) * dx
            yc = (dble(j) - 0.5D0) * dy
            uc = 0.5D0 * (u(i,j) + u(i+1,j))
            vc = 0.5D0 * (v(i,j) + v(i,j+1)) 
            pc = p(i,j)
            write(10, '(5E16.6)') xc, yc, uc, vc, pc
         end do
      end do
      close(10)

      ! Memory deallocation
      deallocate(u, ustar, uold, du)
      deallocate(v, vstar, vold, dv)
      deallocate(p, pstar, pprime, bprime)

      stop
      end