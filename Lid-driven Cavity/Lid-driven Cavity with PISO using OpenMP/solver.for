      program PISO_Lid_Driven_Cavity_omp
      use omp_lib
      implicit none
      
      ! Declaration of variables
      integer*4 nx, ny, maxiter, i, j, k, iter
      parameter(nx = 100, ny = 100, maxiter = 500000)
      
      real*8 tol, start_time, end_time
      parameter(tol = 1.0D - 6)
      
      real*8 rho, mu, ulid, dx, dy, Lx, Ly
      parameter(rho = 1.0D0, mu = 0.01D0, ulid = 1.0D0)
      parameter(Lx = 1.0D0, Ly = 1.0D0)
      
      real*8 alph_p, alph_u, alph_v
      parameter(alph_p = 1.0D0, alph_u = 0.7D0, alph_v = 0.7D0)
      
      real*8, allocatable ::u(:,:),ustar(:,:),ustarstar(:,:),uold(:,:)
      real*8, allocatable ::v(:,:),vstar(:,:),vstarstar(:,:),vold(:,:)
      real*8, allocatable ::p(:,:),pstar(:,:),pstarstar(:,:)
      
      real*8, allocatable :: pprime(:,:), pprime2(:,:)
      real*8, allocatable :: bprime(:,:), bprime2(:,:)
      
      real*8, allocatable :: uprime(:,:), vprime(:,:)
      real*8, allocatable :: ucorr(:,:), vcorr(:,:)
      real*8, allocatable :: du(:,:), dv(:,:)
      
      real*8 Fe, Fw, Fn, Fs, delta_F
      real*8 De, Dw, Dn, Ds
      real*8 ae, aw, an, as, ap, bmax
      real*8 xc, yc, uc, vc, pc
      

      ! Memory allocation
      allocate(u(0:nx+1, 0:ny+1), ustar(0:nx+1, 0:ny+1), 
     +  ustarstar(0:nx+1, 0:ny+1))
      allocate(v(0:nx+1, 0:ny+1), vstar(0:nx+1, 0:ny+1), 
     + vstarstar(0:nx+1, 0:ny+1))
      allocate(uold(0:nx+1, 0:ny+1), vold(0:nx+1, 0:ny+1))
      
      allocate(p(0:nx+1, 0:ny+1), pstar(0:nx+1, 0:ny+1), 
     + pstarstar(0:nx+1, 0:ny+1))
      allocate(pprime(0:nx+1, 0:ny+1), pprime2(0:nx+1, 0:ny+1))
      allocate(bprime(0:nx+1, 0:ny+1), bprime2(0:nx+1, 0:ny+1))
      
      allocate(uprime(0:nx+1, 0:ny+1), vprime(0:nx+1, 0:ny+1))
      allocate(ucorr(0:nx+1, 0:ny+1), vcorr(0:nx+1, 0:ny+1))
      allocate(du(0:nx+1, 0:ny+1), dv(0:nx+1, 0:ny+1))
      
      dx = Lx / dble(nx)
      dy = Ly / dble(ny)
      
      ! Initialization
      u = 0.0D0
      v = 0.0D0
      p = 0.0D0
      ustar = 0.0D0
      vstar = 0.0D0
      pstar = 0.0D0
      
      ! Applying boundary layer (upper lid velocity)
      do i = 0, nx+1
          u(i, ny+1) = ulid
      enddo
      
      ! Start the clock
      start_time = omp_get_wtime()
      
      ! PISO algorithm starts here
      do iter = 1, maxiter
          
          ! Save old velocities
!$omp parallel do private(j)
          do i = 0, nx+1
              do j = 0, ny+1
                  uold(i, j) = u(i, j)
                  vold(i, j) = v(i, j)
              enddo
          enddo
          
          ! Restart matrix of correction
!$omp parallel do private(j)
          do i = 0, nx+1
              do j = 0, ny+1
                  uprime(i, j) = 0.0D0
                  vprime(i, j) = 0.0D0
                  ucorr(i ,j) = 0.0D0
                  vcorr(i, j) = 0.0D0
              enddo
          enddo
          
          De = mu * (dy / dx)
          Dw = De
          Dn = mu * (dx / dy)
          Ds = Dn
          
c         Predictor step
          ! u-momentum predictor
!$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, 
!$omp+ aw, ae, an, as, ap)
          do i = 2, nx
              do j = 1, ny
                  Fe = 0.5D0 * rho * dy * (u(i, j) + u(i+1, j))
                  Fw = 0.5D0 * rho * dy * (u(i, j) + u(i-1, j))
                  Fn = 0.5D0 * rho * dx * (v(i-1, j+1) + v(i, j+1))
                  Fs = 0.5D0 * rho * dx * (v(i-1, j) + v(i, j))
                  delta_F = Fe - Fw + Fn - Fs
                  
                  aw = Dw + max(Fw, 0.0D0)
                  ae = De + max(-Fe, 0.0D0)
                  an = Dn + max(-Fn, 0.0D0)
                  as = Ds + max(Fs, 0.0D0)
                  
                  ap = (aw + ae + an + as + delta_F) / alph_u
                  du(i, j) = dy / ap
                  ustar(i,j) = (aw * uold(i-1,j) + ae * uold(i+1,j) + 
     +                        as * uold(i,j-1) + an * uold(i,j+1) + 
     +                        (pstar(i-1,j) - pstar(i,j)) * dy +
     +                        (1.0D0 - alph_u) * ap * uold(i,j)) / ap
              enddo
          enddo
          
          ! v-momentum predictor
!$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, 
!$omp+ aw, ae, as, an, ap)
          do i = 1, nx
              do j = 2, ny
                  Fe = 0.5D0 * rho * dy * (u(i+1, j-1) + u(i+1, j))
                  Fw = 0.5D0 * rho * dy * (u(i, j-1) + u(i, j))
                  Fn = 0.5D0 * rho * dx * (v(i, j+1) + v(i, j))
                  Fs = 0.5D0 * rho * dx * (v(i, j-1) + v(i, j))
                  delta_F = Fe - Fw + Fn - Fs
                  
                  aw = Dw + max(Fw, 0.0D0)
                  ae = De + max(-Fe, 0.0D0)
                  as = Ds + max(Fs, 0.0D0)
                  an = Dn + max(-Fn, 0.0D0)
                  
                  ap = (aw + ae + as + an + delta_F) / alph_v
                  dv(i, j) = dx / ap
                  vstar(i,j) = (aw * vold(i-1,j) + ae * vold(i+1,j) + 
     +                        as * vold(i,j-1) + an * vold(i,j+1) + 
     +                        (pstar(i,j-1) - pstar(i,j)) * dx +
     +                        (1.0D0 - alph_v) * ap * vold(i,j)) / ap
                  
              enddo
          enddo
          
          ! Applying boundary conditions to u* and v*
          do i = 0, nx+1
              ustar(i, 0)    = 0.0D0          ! Bottom BC
              vstar(i, 1)    = 0.0D0          ! Bottom BC
              ustar(i, ny+1) = ulid           ! Top BC
              vstar(i, ny+1) = 0.0D0          ! Top BC
          enddo
          
          do j = 0 , ny+1
              ustar(1, j)    = 0.0D0          ! Left BC
              vstar(0, j)    = 0.0D0          ! Left BC
              ustar(nx+1, j) = 0.0D0          ! Right BC
              vstar(nx+1, j) = 0.0D0          ! Right BC
          enddo
          
c         First corrector step
          ! Calculate bprime
          bmax = 0.0D0
!$omp parallel do private(j) reduction(max:bmax)
          do i = 1, nx
              do j = 1, ny
                  pprime(i,j) = 0.0D0
                  bprime(i,j) = (rho * dy * ustar(i,j)) +
     +                        (rho * dx * vstar(i,j)) -
     +                        (rho * dy * ustar(i+1,j)) -
     +                        (rho * dx * vstar(i,j+1))
                
                  if (abs(bprime(i,j)) > bmax) then
                      bmax = abs(bprime(i,j))
                  endif
              enddo
          enddo
          
          ! First pressure correction equation 
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
     +                               aw * pprime(i-1,j) + 
     +                               an * pprime(i,j+1) + 
     +                               as * pprime(i,j-1) + 
     +                               bprime(i,j)) / ap
                  enddo
              enddo
              pprime(1, 1) = 0.0D0    ! Pin pressure
          enddo
          
          ! First field update
!$omp parallel do private(j)
          do i = 1, nx
              do j = 1, ny
                  pstarstar(i, j) = pstar(i, j) + pprime(i, j)
                  
                  if (i > 1) then
                      uprime(i,j) =du(i,j)*(pprime(i-1,j) - pprime(i,j))
                      ustarstar(i,j) = ustar(i,j) + uprime(i,j)
                  endif
                  
                  if (j > 1) then
                      vprime(i,j) =dv(i,j)*(pprime(i,j-1) - pprime(i,j))
                      vstarstar(i,j) = vstar(i,j) + vprime(i,j)   
                  endif
              enddo
          enddo
          
          ! Applying boundary conditions for u** and v**
          do i = 0, nx+1
              ustarstar(i, 0)    = 0.0D0          ! Bottom BC
              vstarstar(i, 1)    = 0.0D0          ! Bottom BC
              ustarstar(i, ny+1) = ulid           ! Top BC
              vstarstar(i, ny+1) = 0.0D0          ! Top BC
          enddo
          
          do j = 0 , ny+1
              ustarstar(1, j)    = 0.0D0          ! Left BC
              vstarstar(0, j)    = 0.0D0          ! Left BC
              ustarstar(nx+1, j) = 0.0D0          ! Right BC
              vstarstar(nx+1, j) = 0.0D0          ! Right BC
          enddo
          
c         Second predictor step
!$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, 
!$omp+ aw, ae, as, an, ap)
          do i = 2, nx
              do j = 1, ny
                  Fe = 0.5D0 * rho * dy * (uold(i,j) + uold(i+1,j))
                  Fw = 0.5D0 * rho * dy * (uold(i-1,j) + uold(i,j))
                  Fn = 0.5D0 * rho * dx * (vold(i-1,j+1) + vold(i,j+1))
                  Fs = 0.5D0 * rho * dx * (vold(i-1,j) + vold(i,j))
                  delta_F = Fe - Fw + Fn - Fs
                  
                  aw = Dw + max(Fw, 0.0D0)
                  ae = De + max(-Fe, 0.0D0)
                  as = Ds + max(Fs, 0.0D0)   
                  an = Dn + max(-Fn, 0.0D0)  
                  ap = (aw + ae + as + an + delta_F) / alph_u
               
                  ucorr(i,j) = (aw*uprime(i-1,j) + ae*uprime(i+1,j) + 
     +                       as*uprime(i,j-1) + an*uprime(i,j+1)) / ap
              enddo
          enddo
          
!$omp parallel do private(j, Fe, Fw, Fn, Fs, delta_F, 
!$omp+ aw, ae, as, an, ap)
          do i = 1, nx
              do j = 2, ny
                  Fe = 0.5D0 * rho * dy * (uold(i+1,j-1) + uold(i+1,j))
                  Fw = 0.5D0 * rho * dy * (uold(i,j-1) + uold(i,j))
                  Fn = 0.5D0 * rho * dx * (vold(i,j) + vold(i,j+1))
                  Fs = 0.5D0 * rho * dx * (vold(i,j-1) + vold(i,j))
                  delta_F = Fe - Fw + Fn - Fs
                  
                  aw = Dw + max(Fw, 0.0D0)
                  ae = De + max(-Fe, 0.0D0)
                  as = Ds + max(Fs, 0.0D0)
                  an = Dn + max(-Fn, 0.0D0)
                  ap = (aw + ae + as + an + delta_F) / alph_v
               
                  vcorr(i,j) = (aw*vprime(i-1,j) + ae*vprime(i+1,j) + 
     +                       as*vprime(i,j-1) + an*vprime(i,j+1)) / ap
              enddo
          enddo
          
          ! Second mass source
!$omp parallel do private(j)
          do i = 1, nx
              do j = 1, ny
                  pprime2(i, j) = 0.0D0
                  bprime2(i,j) = (rho * dy * ucorr(i,j)) +
     +                           (rho * dx * vcorr(i,j)) -
     +                           (rho * dy * ucorr(i+1,j)) -
     +                           (rho * dx * vcorr(i,j+1))
              enddo
          enddo
          
          ! Second pressure correction equation
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
                  
                      pprime2(i,j) = (ae * pprime2(i+1,j) + 
     +                                aw * pprime2(i-1,j) + 
     +                                an * pprime2(i,j+1) + 
     +                                as * pprime2(i,j-1) + 
     +                                bprime2(i,j)) / ap   
                  enddo
              enddo
              pprime2(1, 1) = 0.0D0
          enddo
          
          ! Final update u***, v***
!$omp parallel do private(j)
          do i = 1, nx
              do j = 1, ny
                  p(i,j) = pstarstar(i,j) + pprime2(i,j)
                  pstar(i,j) = p(i,j)
               
                  if (i > 1) then
                      u(i,j) = ustarstar(i,j) + ucorr(i,j) + 
     +                     du(i,j)*(pprime2(i-1,j) - pprime2(i,j))
                  endif
               
                  if (j > 1) then
                      v(i,j) = vstarstar(i,j) + vcorr(i,j) + 
     +                     dv(i,j)*(pprime2(i,j-1) - pprime2(i,j))
                  endif            
              enddo
          enddo
          
          write(*,*) 'Iteration: ', iter, '|Max mass error: ', bmax
          
          if (bmax <= tol) exit
          
      enddo       ! PISO algorithm stops here
      
      ! Stop the clock
      end_time = omp_get_wtime()
      
      write(*,*) '------------------------------------'
      write(*,*) 'Total iterations: ', iter
      write(*,*) 'Total time (s): ', end_time - start_time
      write(*,*) '------------------------------------'

      ! Outputting solution
      open(unit=10, file='cavity_PISO_omp.txt', status='unknown')
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
      
      
      ! Memory deallocation
      deallocate(u, ustar, ustarstar, uold)
      deallocate(v, vstar, vstarstar, vold)
      deallocate(p, pstar, pstarstar)
      deallocate(pprime, pprime2, bprime, bprime2)
      deallocate(uprime, vprime, ucorr, vcorr)
      deallocate(du, dv)
      
      stop
      end