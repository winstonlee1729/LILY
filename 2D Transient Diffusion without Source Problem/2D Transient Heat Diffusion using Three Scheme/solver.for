      program Transient_Heat_Diffusion
      use omp_lib
      implicit none

      ! Declaration of variables
      integer*4 i, j, iter, nx, ny 
      integer*4 time_step, max_step, maxiter, iter_limit
      integer*4 output_freq
      parameter(nx = 300, ny = 300) 
      parameter(maxiter = 5000, max_step = 2000)
      
      real*8 maxres, err, tol, rhs, lhs
      parameter(tol = 1.0D-6)
      
      real*8, allocatable :: T(:,:), Told(:,:)
      real*8 ke, kw, kn, ks, Area_e, Area_w, Area_s, Area_n
      real*8 dx, dy, dt, Lx, Ly, ap, ae, aw, an, as, ap_o, ap_prime
      real*8 rho, Cp, theta, volume, thickness
      real*8 w_node, e_node, n_node, s_node, p_past
      
      parameter(ke = 5.0D0, kw = 5.0D0, ks = 5.0D0, kn = 5.0D0)
      parameter(Lx = 1.0D0, Ly = 1.0D0, thickness = 1.0D0)
      parameter(dt = 0.05D0)
      parameter(rho = 1000.0D0, Cp = 1000.0D0)
      
      ! 1.0 = Fully Implicit | 0.5 = Crank-Nicolson | 0.0 = Explicit
      parameter(theta = 0.5D0) 
      
      ! Output frequency: Save data to file every 50 time steps
      parameter(output_freq = 50)
      
      ! Defining variables
      dx = Lx / dble(nx)
      dy = Ly / dble(ny)
      volume = dx * dy * thickness
      
      Area_e = dy * thickness
      Area_w = Area_e
      Area_n = dx * thickness
      Area_s = Area_n
      
      ap_o = rho * Cp * (volume / dt)
      ae = (ke * Area_e) / dx
      aw = (kw * Area_w) / dx
      an = (kn * Area_n) / dy
      as = (ks * Area_s) / dy
      
      ap = ap_o + (aw + ae + an + as) * theta
      ap_prime = ap_o - (ae + aw + an + as) * (1.0D0 - theta)
      
      ! Allocate memory
      allocate(T(nx, ny), Told(nx, ny))
      
      ! Initialize flow field
      T = 0.0D0
      Told = 0.0D0
      
      ! Applying boundary condition 
!$omp parallel do
      do i = 1, nx
          T(i, ny) = 1000.0D0
          Told(i, ny) = 1000.0D0
      enddo
!$omp end parallel do
      
      ! Check for iteration number using different schemes
      if (theta == 0.0D0) then
          iter_limit = 1
      else
          iter_limit = maxiter
      endif
      
      ! Open file to write data
      open(unit=10, file='T_out.txt', status='replace')

      ! Start time iteration
      do time_step = 1, max_step
          
          ! Save temperature field
!$omp parallel do collapse(2)
          do j = 1, ny
              do i = 1, nx
                  Told(i, j) = T(i, j)
              enddo
          enddo
!$omp end parallel do
          
          ! Inner iteration
          do iter = 1, iter_limit
              maxres = 0.0D0

!$omp parallel do private(i, e_node, w_node, n_node, 
!$omp+ s_node, p_past, rhs, lhs, err) reduction(max:maxres)
              do j = 2, ny-1
                  do i = 2, nx-1
                      e_node = ae * (theta * T(i+1, j) + (1.0D0 - theta)
     +                        * Told(i+1, j))
                      w_node = aw * (theta * T(i-1, j) + (1.0D0 - theta)
     +                        * Told(i-1, j))
                      n_node = an * (theta * T(i, j+1) + (1.0D0 - theta)
     +                        * Told(i, j+1))
                      s_node = as * (theta * T(i, j-1) + (1.0D0 - theta)
     +                        * Told(i, j-1))
                      p_past = ap_prime * Told(i, j)
                             
                      rhs = e_node + w_node + n_node + s_node + p_past
                      lhs = ap * T(i, j)
                          
                      err = abs(rhs - lhs)
                      if (err > maxres) maxres = err
                      
                      ! Update temperature field
                      T(i, j) = rhs / ap
                  enddo
              enddo
!$omp end parallel do

              if (theta > 0.0D0 .and. maxres < tol) exit
          enddo
          
          if (theta == 0.0D0) then
              write(*,*) 'Time step:', time_step
          else
              write(*,*) 'Time step:', time_step, '|Iters:', iter
     +                    ,'|Res:', maxres
          endif 
          
          ! Outputting solution
          if (mod(time_step, output_freq) == 0 .or. time_step == 1) then
              do j = 1, ny
                  write(10, '(1000E15.6)') (T(i, j), i=1, nx)
              enddo
          endif
    
      enddo
      
      ! Close the file
      close(10)
      
      ! Deallocate memory
      deallocate(T, Told)
     
      stop
      end
                  
                  
              
          
          
      
      