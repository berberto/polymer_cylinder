module module_ana_dat_diffusion

  integer, parameter :: nreps = 1
  integer, parameter :: npoints = 10000
  real(8), parameter :: z01=2.4048, R=1.d0, pi=acos(-1.0d0)
  real(8) :: w, lambda, m
 

contains

  ! subroutine averages
  
  !   implicit none

!   !integer, parameter :: nreps = 20
!   integer, parameter :: npoints = 1000000
!   integer :: seed, k
!   integer(4) :: values(13),file_results
!   integer(4) :: f,i,blocksize,fsize,c
!   logical :: exist
!   real(4) :: tmp(3*npoints), x(3)
!   real(4) :: points(npoints,3),over_n
!   !real(4) :: t(npoints-1,3),avg_t(npoints-1,3)
!   real(4) :: normt,Ree
!   character(100) :: arg1,main_input_dir
!   real(4) :: avg_L(npoints-1),avg_L2(npoints-1),&
!        avg_Lz(npoints-1),avg_Lz2(npoints-1),&
!        avg_dl(npoints-1),&
!        avg_dz(npoints-1),&
!        avg_Rg2(npoints-1),avg_Rg4(npoints-1),&
!        avg_Ree(npoints-1),avg_Ree2(npoints-1),&
!        L(npoints-1),&
!        Lz(npoints-1),&
!        nth_Ree(npoints-1),&
!        nth_Rg2(npoints-1),&
!        dz(npoints-1),&
!        dl(npoints-1)






!   !call getarg(1,arg1)
!   !read(arg1,*), avgjmp
!   !main_input_dir = '/media/USB-HDD/data_cylinder/avjmp_'//trim(arg1)//'/'
!   !main_input_dir = '/scratch/madorisi/results/polymers/avjmp_'//trim(arg1)//'/'
!   main_input_dir = './results/'
!   !main_output_dir = '../output/'
!   !  print*, 'working on ',trim( main_input_dir )


!   avg_Ree = 0. !array
!   avg_Ree2 = 0. !array
!   avg_Rg2 = 0. !array
!   avg_Rg4 = 0. !array
!   avg_L = 0. ! array
!   avg_L2 = 0. !array
!   avg_Lz = 0. !array
!   avg_Lz2 = 0. !array
!   avg_dz = 0.
!   avg_dl = 0.


!   c = 0


!   do f = 1,nreps
!      !print*, f
!      !check file
!      !call stat(trim(main_input_dir)//'rep_'//trim( int_2_str(f) )//'.dat', values) 
!      !fsize=values(8)

!      !if (fsize==0) then !in case of zero byte files
!      !c = c + 1
!      !cycle
!      !end if

!      ! continue

!      ! open(10,file=trim(main_input_dir)//'traj.dat',&
!      !      access='stream',&
!      !      form='unformatted')

!      ! read(10), seed
!      ! read(10), tmp
!      ! close(10)
!      inquire(file='./results/traj_'//trim(int_2_str(f))//'.dat',exist=exist)
!      print*, './results/traj_'//trim(int_2_str(f))//'.dat', exist
     
!      if (exist) then
!         open(10,file='./results/traj_'//trim(int_2_str(f))//'.dat',status='unknown')

!         !--------------------
!         !
!         !--------------------
!         do k=1,npoints
!            read(10,*), points(k,:)  !reshape(tmp,[npoints,3],order=[2,1])
!         end do
!         close(10)
!         ! do i=1,npoints-1
!         !    tot_contour = tot_contour + norm2(points(i+1,:) - points(i,:))
!         ! end do

!         ! avg_L = L(points)
!         ! do i =1,size(avg_L,1)
!         !    print*, avg_L(i)
!         ! end do
!         ! print*, tot_contour, avg_L(npoints-1)
!         ! stop

!         call measure(&
!              points,& !in
!              L,& !out
!              Lz,&
!              nth_Ree,&
!              nth_Rg2,&
!              dz,&
!              dl)


!         !end_to_end = nth_Ree(points)
!         avg_Ree = avg_Ree + nth_Ree / nreps!norm2(  points(1,:) - points(npoints,:)  )/nreps
!         avg_Ree2 = avg_Ree2 + nth_Ree**2 / nreps
!         !center of mass position
!         !Rcm = sum( points,dim=1 )/npoints 

!         !average end to end distance 

!         !gyration radius
!         !gyradius = nth_Rg2(points)
!         avg_Rg2 = avg_Rg2 + nth_Rg2 / nreps
!         avg_Rg4 = avg_Rg4 + nth_Rg2**2 / nreps
!         !tot length up to n-th point
!         !tot_l = L(points)
!         avg_L = avg_L + L / nreps
!         avg_L2 = avg_L2 + L**2 / nreps
!         !distance travelled on the z-axis up to n-th point
!         !tot_lz =  Lz(points)
!         avg_Lz = avg_Lz + Lz / nreps
!         avg_Lz2 = avg_Lz2 + Lz**2 / nreps
!         avg_dz = avg_dz + dz / nreps
!         avg_dl = avg_dl + dl / nreps

!      end if
!   end do

!   ! do i=1,npoints-1
!   !    avg1_avgi(i) = dot_product( avg_t(1,:), avg_t(i,:) )
!   ! end do
!   print*, 'ciao'
!   call write_results(npoints,&
!        avg_L,avg_L2,&
!        avg_Lz,avg_Lz2,&
!        avg_Rg2,avg_Rg4,&
!        avg_Ree,avg_Ree2,&
!        avg_dl,&
!        avg_dz)


!   print*, trim( main_input_dir ),' done'

!   !write(6,*), c, ' files with zero bytes ', arg1 !6=stdout

!end subroutine averages


  function int_2_str(int) result(str)
    !input a real. If integer cast it to real before
    implicit none
    
    integer,intent(in) :: int
    character(len=200) :: str
    
    write(str,'(I0)'), int
    
  end function int_2_str
  
  
  subroutine measure(points,& !in
       L,& !out
       Lz,&
       nth_Ree,&
       nth_Rg2,&
       dz,&
       dl )
    
    implicit none
    
    real, intent(in) :: points(:,:)
    real, intent(out) :: L(size(points,1)-1),&
         Lz(size(points,1)-1),&
         nth_Ree(size(points,1)-1),&
         nth_Rg2(size(points,1)-1),&
         dz(size(points,1)-1),&
         dl(size(points,1)-1)
    
    integer :: i,k
    real :: sum_r(3), sum_rr, over_n
    
    !init
    L = 0.
    L(1) = norm2(points(2,:) - points(1,:))
    Lz = 0.
    Lz(1) = points(2,3) - points(1,3)
    nth_Ree = 0.
    nth_Ree(1) = L(1)
    dz = 0.
    dz(1) = points(2,3) - points(1,3)
    dl = 0.
    dl(1) = norm2(points(2,:)-points(1,:))
    !nth_Rcm = points(1,:) + points(2,:)
    sum_rr = dot_product(points(1,:),points(1,:)) + dot_product(points(2,:),points(2,:))
    sum_r = points(1,:) + points(2,:)
    over_n = 0.5
    nth_Rg2 = 0.
    nth_Rg2(1) = sum_rr*over_n - dot_product(sum_r,sum_r)*over_n*over_n
    
    do i = 1,npoints-2
       
       L(i+1) = L(i) + norm2( points(i+2,:) - points(i+1,:) )
       Lz(i+1) = Lz(i) + points(i+2,3) - points(i+1,3)
       nth_Ree(i+1) = norm2(  points(i+2,:) - points(1,:)  )
       dz(i+1) = points(i+2,3) - points(i+1,3)
       dl(i+1) = norm2( points(i+2,:) - points(i+1,:) )
       
       sum_r = sum_r + points(i+2,:)
       sum_rr = sum_rr + dot_product(points(i+2,:),points(i+2,:))
       
       over_n = 1./real(i+2)
       nth_Rg2(i+1) = sum_rr*over_n - dot_product(sum_r*over_n,sum_r*over_n)
       !nth_Rg2(i+1) = nth_Rg2(i) + 
       !gyradius
       !do k=1,i+2
       !nth_Rg2(i+1) = nth_Rg2(i+1) + dot_product(points(k,:),points(k,:))
       !end do
       !nth_Rcm = nth_Rcm + sum( points(1:i+1,:),dim=1 ) / real(i+1)
       !nth_Rcm = nth_Rcm + points(i+2,:)
       
       !nth_Rg2(i+1) =  nth_Rg2(i+1) / real(i+1)
       !nth_Rg2(i+1) = nth_Rg2(i+1) - dot_product(nth_Rcm,nth_Rcm)/real(i+2)/real(i+2)
       
    end do
    
    
  end subroutine measure
  
  subroutine write_results(npoints,&
     avg_L,avg_L2,&
     avg_Lz,avg_Lz2,&
     avg_Rg2,avg_Rg4,&
     avg_Ree,avg_Ree2,&
     avg_dl,&
     avg_dz,&
     D_label,&
     avjmp_label)
    
    implicit none
    
    integer, intent(in) :: npoints
    character(len=20) :: label, D_label, avjmp_label
    real(4), intent(in) :: avg_L(npoints-1),avg_L2(npoints-1),&
         avg_Lz(npoints-1),avg_Lz2(npoints-1),&
         avg_dl(npoints-1),&
         avg_dz(npoints-1),&
         avg_Rg2(npoints-1),avg_Rg4(npoints-1),&
         avg_Ree(npoints-1),avg_Ree2(npoints-1)
    
    integer :: i, file_results
    
    !    write(label,'(es12.3)'), avgjmp  !convert to string
    
    !print*, '/scratch/madorisi/results/polymers/results_freejmp_'//trim(label)//'.dat'
    !stop
    
    open(newunit=file_results,&
         file='./results/coarse_grained/results_D'//trim(D_label)//'_avjmp'//trim(avjmp_label)//'.dat',&
         status='unknown')
    write(file_results,*), '1it ', '2avg_Ree ','3var_Ree ',&
         '4avg_Rg2 ', '5var_Rg2 ',&
         '6avg_L ', '7var_L ',&
         '8avg_Lz ','9var_Lz ',&
         '10avg_dl ', &
         '11avg_dz ',&
         '12avgjmp '
    
    do i=1,npoints-1
       !print*, i
       write(file_results,*), i,avg_Ree(i),(avg_Ree2(i)-avg_Ree(i)**2)/real(nreps),&
            avg_Rg2(i),(avg_Rg4(i)-avg_Rg2(i)**2)/real(nreps),&
            avg_L(i),(avg_L2(i)-avg_L(i)**2)/real(nreps),&
            avg_Lz(i),(avg_Lz2(i)-avg_Lz(i)**2)/real(nreps),&
            avg_dl(i), avg_dz(i),&
            sum(avg_dl)/size(avg_dl)
    end do
    close(file_results)
    
    
    ! inquire (file='./results/polymers/results_Ree_Rg.dat', exist=exist)
    ! if (exist) then
    !    open(newunit=file_results, file='/scratch/madorisi/results/polymers/results_Ree_Rg.dat',&
    !         status="old", position="append", action="write")
    !    write(file_results,*), avg_Rg2, avg_Ree, avg_L(npoints-1), trim(label)!avgjmp
    !    close(file_results)
    ! else
    !    !create file and write the first time
    !    open(newunit=file_results, file='./results/polymers/results_Ree_Rg.dat', status="new", action="write")
    !    write(file_results,*), '# avgRg2 ','avgRee ','avgL ','avgjmp'
    !    write(file_results,*), avg_Rg2, avg_Ree, avg_L(npoints-1), trim(label)!avgjmp
    !    close(file_results)
    ! end if
    
    
    
    
  end subroutine write_results
  
  
  

  ! function costheta(t,npoints)
  !   implicit none

  !   real,intent(in) :: t(:,:)
  !   integer, intent(in) :: npoints
  !   real :: costheta(npoints-1)

  !   integer :: j

  !   do i = 1,npoints-1
  !      costheta(i) = dot_product( t(1,:) , t(i,:) )  
  !   end do

  ! end function costheta


  ! function Rg2(points,Rcm)
  !   implicit none

  !   real, intent(in) :: points(:,:), Rcm(3)
  !   real :: Rg2
  !   integer :: i

  !   Rg2 = 0.
  !   do i=1,size(points,1)
  !      Rg2 = Rg2 + dot_product(points(i,:),points(i,:))
  !   end do

  !   Rg2 = Rg2/size(points,1)
  !   Rg2 = Rg2 - dot_product(Rcm,Rcm)

  ! end function Rg2

  ! ! function contour(points,upto)
  ! !   implicit none

  ! !   real, intent(in) :: points(:,:)
  ! !   integer, intent(in) :: upto
  ! !   real :: contour
  ! !   integer :: i

  ! !   contour = 0.


  ! !   do i=1,upto
  ! !      contour = contour + norm2( points(i+1,:) - points(i,:) )
  ! !   end do




  ! ! end function contour
  ! function L(points)
  !   implicit none

  !   real, intent(in) :: points(:,:)
  !   real :: L(size(points,1)-1)

  !   L = 0.
  !   L(1) = norm2(points(2,:) - points(1,:))

  !   do i = 1,npoints-2
  !      L(i+1) = L(i) + norm2( points(i+2,:) - points(i+1,:) )
  !   end do

  ! end function L

  ! function Lz(points)
  !   implicit none
  !   !compute the jump along z 
  !   !Lz(i): jump on z up to the i-th jump

  !   real, intent(in) :: points(:,:)
  !   real :: Lz(size(points,1)-1)

  !   Lz = 0.
  !   Lz(1) = points(2,3) - points(1,3)

  !   do i = 1,npoints-2
  !      Lz(i+1) = Lz(i) + points(i+2,3) - points(i+1,3) 
  !   end do

  ! end function Lz

  ! function nth_Ree(points)
  !   implicit none
  !   !compute the jump along z 
  !   !Lz(i): jump on z up to the i-th jump

  !   real, intent(in) :: points(:,:)
  !   real :: nth_Ree(size(points,1)-1)

  !   nth_Ree = 0.

  !   do i = 1,npoints-1
  !      nth_Ree(i) = norm2(  points(i+1,:) - points(1,:)  )
  !   end do

  ! end function nth_Ree

  ! function nth_Rg2(points)
  !   implicit none

  !   real, intent(in) :: points(:,:)
  !   real :: nth_Rg2(npoints-1),nth_Rcm(3)
  !   integer :: np,k

  !   nth_Rg2 = 0.

  !   do np=2,size(points,1)

  !      nth_Rcm = sum( points(1:np,:),dim=1 ) / real(np)

  !      do k=1,np
  !         nth_Rg2(np-1) = nth_Rg2(np-1) + dot_product(points(k,:),points(k,:))
  !      end do

  !      nth_Rg2(np-1) = nth_Rg2(np-1)/real(np)
  !      nth_Rg2(np-1) = nth_Rg2(np-1) - dot_product(nth_Rcm,nth_Rcm)

  !   end do


  ! end function nth_Rg2


  ! function dz(points)
  !   implicit none
  !   !compute the jump along z 
  !   !Lz(i): jump on z up to the i-th jump

  !   real, intent(in) :: points(:,:)
  !   real :: dz(size(points,1)-1)

  !   dz = 0.

  !   do i = 1,npoints-1
  !      dz(i) = points(i+1,3) - points(i,3) 
  !   end do

  ! end function dz


  ! function dl(points)
  !   implicit none
  !   !compute the jump along z 
  !   !Lz(i): jump on z up to the i-th jump

  !   real, intent(in) :: points(:,:)
  !   real :: dl(size(points,1)-1)

  !   dl = 0.

  !   do i = 1,npoints-1
  !      dl(i) = norm2( points(i+1,:) - points(i,:) ) 
  !   end do

  ! end function dl


  function f_rho(rho)  

    use mybessel
    
    real(8), intent (in) :: rho
    !real(8), parameter :: z01=2.4048,R=1.
    !real(8) :: lambda
    
    !lambda = 1.00076!z01/R
!    print*, 'w',w
    !cumulative of p_st
    f_rho = (rho/R)**2*( besj0( lambda*rho )**2 + besj1( lambda*rho )**2 )/ (besj0( lambda*R )**2 + besj1(lambda*R)**2)  - w

  end function f_rho
  
  function f_lambda(lambda)  
    
    use mybessel
    
    real(8), intent (in) :: lambda
    !real(8), parameter :: z01=2.4048,R=1.
    !real(8) :: lambda
    !print*,'m', m, lambda
!    lambda = 1.00076!z01/R
    
    u = sqrt(m*m - lambda*lambda);
    f_lambda = lambda*besj1(lambda*R) - u*besj0(lambda*R)*besek1(R*u)/besek0(R*u);
    !cumulative of p_st
     !rho/R)**2*( besj0( lambda*rho )**2 + besj1( lambda*rho )**2 )/ (besj0( lambda*R )**2 + besj1(lambda*R)**2)  - w
    
  end function f_lambda
  

  subroutine bisect(f,a,b,w,TOL,x)                                            

    implicit none

    real(8), intent(in) :: TOL,w
    real(8), intent(inout) :: a,b !inout 'cause they are modified here inside
    real(8), intent(out) :: x
    integer :: i
    real(8) :: f0,f1,fx

    interface                                                             
       function f(x)                                                   
         real(8), intent(in) :: x
       end function f
    end interface

    
    f0 = f(a) 
    f1 = f(b)
   
    
    if ( f0==f1 .or. (f0*f1>0) ) then
       print*, 'impossible to use bisection'
       stop
    end if

    i=0
    do while ( abs(a-b) > TOL )
       
       i = i + 1
       x = 0.5*(a + b)
       
       fx = f(x)

       if ( f0*fx < 0 ) then
          b = x
          f1 = fx
       else if (f0*fx > 0) then
          a = x
          f0 = fx
       else
          stop
       end if
       
       if (i==200) then
          print*, 'bisection does not converge'
       end if


    end do

  end subroutine bisect

  function rotate2d(point,theta,origin) result(new_point)
    implicit none

    real(8), dimension(2), intent(in) :: point,origin
    real(8),intent(in) :: theta
    real(8), dimension(2):: new_point
    real :: s,c
    new_point = [0.,0.]

    c = cos(theta)
    s = sin(theta)

    new_point(1) = (point(1) - origin(1))*c - (point(2) - origin(2))*s + origin(1)
    new_point(2) = (point(1) - origin(1))*s + (point(2) - origin(2))*c + origin(2)

  end function rotate2d


  function relative_angle(v1,v2) result(ang)
    implicit none
    !compute the angle of v1 respect to v2 (in 2D)
    ! v1 and v2 must be normalized
    real(8), intent(in) :: v1(2),v2(2)
    real(8) :: ang
    
    ang = acos( dot_product(v1,v2) )
    if (cross(v2,v1)>0) then
       ang = 2.*pi - ang
    end if
    
  end function relative_angle

  function cross(a,b) result(c_z)
    
    implicit none
    
    real(8),intent(in) :: a(2),b(2)
    real :: c_z
    
    c_z = a(1)*b(2)-a(2)*b(1)
    
    
  end function cross
  
  function intersection_point_par(x,x_new) result(int_point)
    implicit none
    
    real(8), intent(in) ::  x(2),x_new(2)
    real(8) :: dx(2), a, c, b, s
    real(8) :: i
    real(8) :: int_point(2)

    dx = x_new-x
    a = dot_product(dx,dx)
    b = dot_product(dx,x)
    c = - 1.d0 + dot_product(x,x) 
    
    s = (-b + sqrt( b**2 - a*c ) ) / a

    int_point = x + (x_new-x)*s

  end function intersection_point_par


  function intersection_point(x,x_new) result(c)

    implicit none
    real(8), intent(in) ::  x(2),x_new(2)
    real(8) :: c(2), xc1, xc2, xc, yc, xmin, xmax
    real :: m,m2, t
    real :: d1(2), d2(2)
    real :: sign

    m = ( x_new(2) - x(2) )/(  x_new(1) - x(1)   )
    
    t = x(2) - m*x(1)
    
    m2 = 1.+m**2

    xc1 = ( -m*t + sqrt(  m2*R**2 - t**2 ) ) / m2
    xc2 = ( -m*t - sqrt(  m2*R**2 - t**2 ) ) / m2
!    y1c = *sqrt( R**2 - x1c**2 )
!    y2c = *sqrt( R**2 - x2c**2 )
    

    xmin = x(1)
    xmax = x_new(1)
    if ( x(1)>x_new(1) ) then
       xmin = x_new(1)
       xmax = x(1)
    end if

    xc = xc1
    if ( xc2 > xmin .and. xc2 < xmax) xc = xc2
    sign = m*(xc-x(1)) + x(2)
    sign =  sign/abs(sign)
    yc = sign*sqrt( R**2 - xc**2 )
    

    ! d1(1) = norm2([x1c,y1c]-x) ; d1(2) = norm2([x2c,y2c]-x)
    ! d2(1) = norm2([x1c,y1c]-x_new) ; d2(2) = norm2([x2c,y2c]-x_new)
    
    ! if ( norm2([x1c,y1c]-x) > norm2([x2c,y2c]-x)) then
    !    c = [ x2c, y2c ]
    ! else
    !    c = [ x1c, y1c ]
    ! end if
    
    c = [xc,yc]
    
  end function intersection_point
  
  function reflect(x,x_new) result(p)
    implicit none

    real(8), intent(in) :: x(3), x_new(3)
    real(8) :: x_(2), x_new_(2)
    real(8) :: p(3)
    real(8) :: c(2),theta
    
    x_ = x(1:2)
    x_new_ = x_new(1:2)
    !c = intersection_point( x_, x_new_ ) 
    c = intersection_point_par( x_, x_new_ )
!    print*, 'c',c
!    write(10,*), c,'3'
    theta = relative_angle( c, [0.d0,1.d0] )
    x_ = rotate2d(x_, theta, [0.d0,0.d0])
    x_new_ = rotate2d(x_new_, theta, [0.d0,0.d0])
    x_new_ = [x_new_(1),2.*R-x_new_(2)]
    x_new_ = rotate2d(x_new_, -theta, [0.d0,0.d0])
    
    p = [x_new_(1), x_new_(2), x_new(3)]
    
  end function reflect
  
  ! subroutine seed_from_urandom(i)
  !   implicit none
  !   integer :: n,k
  !   integer(4),intent(out),allocatable :: i(:)
  !   real(8) :: r
  !   call random_seed(size=n)
  !   allocate(i(n))
  !   open(89, file='/dev/urandom', access='stream', form='UNFORMATTED')
  !   read(89), i
  !   close(89)
    
  !   call random_seed(put=i)
    
  ! end subroutine seed_from_urandom

end module module_ana_dat_diffusion
