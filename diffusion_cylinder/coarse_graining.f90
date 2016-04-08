program coarse_grained_diffusion

  use random
  use mybessel
  use module_ana_dat_diffusion


  implicit none
  real(8), parameter :: TOL = 0.000001d0
  !, pi=acos(-1.0d0)
  real(8) :: x(3), x0(3),x_old(3),tmp_x(3), rho, z, D, dt, u_x, u_y, u_z, ur, sqdt, f, tau, avjmp, tau_arr(npoints), t
  integer :: it, nit,i
  integer, allocatable :: seed(:)
  real(8) :: rho0, a , b, theta, t_save, u(2)
  character(50) :: avjmp_label, outdir, rep, fname
  real :: points(npoints,3)
  !for ana_dat
  !logical :: to_measure
  
  ! real(4) :: avg_L(npoints-1),avg_L2(npoints-1),&
  !      avg_Lz(npoints-1),avg_Lz2(npoints-1),&
  !      avg_dl(npoints-1),&
  !      avg_dz(npoints-1),&
  !      avg_Rg2(npoints-1),avg_Rg4(npoints-1),&
  !      avg_Ree(npoints-1),avg_Ree2(npoints-1),&
  !      L(npoints-1),&
  !      Lz(npoints-1),&
  !      nth_Ree(npoints-1),&
  !      nth_Rg2(npoints-1),&
  !      dz(npoints-1),&
  !      dl(npoints-1)



  call seed_from_urandom(seed)

  

  call getarg(1,avjmp_label)
  read(avjmp_label,*), avjmp

  call getarg(2,outdir)
  call system("mkdir -p "//trim(outdir))

  call getarg(3,rep)
  fname=trim(outdir)//"/rep_"//trim(rep)//".dat"
  
  !read(avjmp_label,*), avjmp

  
  
  m = 2./R/avjmp; a = 0.000001*m; b = 2.41/R

  if (m < b) b = m - 1.e-10

  call bisect(f_lambda,a,b,m,TOL,lambda)


  !print*, m
  !  print*, f_rho(1.d0)


  

  !the simulation for the diffusion process has been performed with Ddt=10^-5
  D = 0.001
  dt = 0.001
  tau = avjmp**2/4./D
  

  u_z = 2.*lambda
  sqdt = sqrt(2.*D*dt)
  
  !L = 0.
  !rep = 0


  ! do i=1,100000
  !    call random_number(ur)
  !    write(100,*), -tau*log(1-ur)
  ! end do
  ! stop
  !  do rep=1,nreps!10000!nreps

  !to_measure = .true.

  a = 0.d0
  b = R
  call random_number(w) !passed to f_rho

  call bisect(f_rho,a,b,w,TOL,rho0)
  call random_number(w)
  theta = 2.*pi*(w-0.5)

  x0 = rho0*[cos(theta),sin(theta),0.d0]
  x = x0
  points(1,:) = x0
  
  it = 1
  call random_number(ur)
  t_save = -tau*log(1-ur)
  t=0
  rho=rho0

  do while (it<=npoints)

     x_old = x

     u_x = -u_z*besj1( lambda*rho )/besj0( lambda*rho )*x(1)/rho
     u_y = -u_z*besj1( lambda*rho )/besj0( lambda*rho )*x(2)/rho

     x(1) = x(1) + u_x*D*dt + sqdt*random_normal()
     x(2) = x(2) + u_y*D*dt + sqdt*random_normal()
     x(3) = x(3) + u_z*D*dt + sqdt*random_normal()

     rho = sqrt( x(1)**2 + x(2)**2 )

     !reflect if outside
     if (rho>R) then
        
!        delta = norm2(x(1:2)) - R
!        u = x(1:2)/norm2(x(1:2))
!        x = [(R - delta)*u(1),(R - delta)*u(2), x(3)]

        tmp_x = reflect( x_old, x )
        x = tmp_x
        rho = sqrt( x(1)**2 + x(2)**2 )
     end if
     
     t = t + dt !elapsed time  
     

     !it = it + 1
     !points(it,:) = x
     !write(40,*), points(it,:)
     
     !sample
     if( t > t_save ) then

        !print*, points(1,:)
        t = 0 !reset elapsed time
        it = it + 1
!        print*, it
        call random_number(ur)
        t_save = -tau*log(1-ur)
        points(it,:) = x
        !           print*, points(2,:)
 !       write(20,*), points(it,:)
        write(100,*), sqrt( points(it,1)**2 + points(it,2)**2   )
     end if
     
  end do
  
  points(1,:) = x0
  !  open(unit=10, form='unformatted', file='prova.bin')
!  open(10, form="unformatted",file=trim(outdir)//"/rep_"//trim(rep)//".dat",status="unknown")
!  
!  open(10, form="unformatted", file=trim(fname),status="unknown")
!  write(10) points
!  close(10)


  !if (to_measure) then
  !        print*, 'measuring...'

  !rep = rep + 1

  !     call measure(&
  ! points,& !in
  ! L,& !out
  ! Lz,&
  ! nth_Ree,&
  ! nth_Rg2,&
  ! dz,&
  ! dl)

  !        print*, 'done measuring'
  !end_to_end = nth_Ree(points)
  !     avg_Ree = avg_Ree + nth_Ree !norm2(  points(1,:) - points(npoints,:)  )/nreps
  !     avg_Ree2 = avg_Ree2 + nth_Ree**2
  !center of mass position
  !Rcm = sum( points,dim=1 )/npoints 

  !average end to end distance 

  !gyration radius
  !gyradius = nth_Rg2(points)
  !     avg_Rg2 = avg_Rg2 + nth_Rg2 !/ nreps
  !     avg_Rg4 = avg_Rg4 + nth_Rg2**2 !/ nreps
  !tot length up to n-th point
  !tot_l = L(points)
  !     avg_L = avg_L + L !/ nreps
  !     avg_L2 = avg_L2 + L**2 !/ nreps
  !distance travelled on the z-axis up to n-th point
  !tot_lz =  Lz(points)
  !     avg_Lz = avg_Lz + Lz !/ nreps
  !     avg_Lz2 = avg_Lz2 + Lz**2 !/ nreps
  !     avg_dz = avg_dz + dz !/ nreps
  !     avg_dl = avg_dl + dl !/ nreps

  !     close(10)

  !  end do
  !print*, 'ciao'
  !print*, rep, nreps
  ! print*, 'writing ...'

  !  call write_results(npoints,&
  !       avg_L/rep,avg_L2/rep,&
  !       avg_Lz/rep,avg_Lz2/rep,&
  ! avg_Rg2/rep,avg_Rg4/rep,&
  ! avg_Ree/rep,avg_Ree2/rep,&
  ! avg_dl/rep,&
  ! avg_dz/rep,&
  ! avjmp_label)


  ! print*, 'done writing'



  ! contains

  !   function int_2_str(int) result(str)
  !   !input a real. If integer cast it to real before
  !     implicit none

  !     integer,intent(in) :: int
  !     character(len=200) :: str

  !     write(str,'(I0)'), int

  !   end function int_2_str



end program coarse_grained_diffusion
