program ana_dat

  implicit none

  integer, parameter :: nreps = 10000
  integer, parameter :: npoints = 10000
  integer :: seed
  integer(4) :: values(13),file_results
  integer(4) :: f,i,blocksize,fsize,c
  real(4) :: tmp(3*npoints)
  real(4) :: points(npoints,3)
  real(4) :: t(npoints-1,3),avg_t(npoints-1,3)
  real(4) :: avgjmp, normt,avg_Ree,Ree,Rcm(3),avg_Rg2
  character(100) :: arg1,main_input_dir
  real(8) :: avg_L(npoints-1), avg_costheta(npoints-1), avg_dz(npoints-1), avg1_avgi(npoints-1)

  call getarg(1,arg1)
  read(arg1,*), avgjmp

  main_input_dir = '/scratch/madorisi/results/polymers/avjmp_'//trim(arg1)//'/'
  !main_output_dir = '../output/'
  print*, 'working on ',trim( main_input_dir )
  

  avg_Ree = 0.
  avg_Rg2 = 0.
  avg_L = 0.
  c = 0
  
  
  do f = 0,nreps-1
     !call stat(trim(main_input_dir)//'rep_'//trim( int_2_str(f) )//'.dat', values) 
     !fsize=values(8)
     
     !if (fsize==0) then !in case of zero byte files
        !c = c + 1
        !cycle
     !end if
     
     ! continue
     
     open(10,file=trim(main_input_dir)//'rep_'//trim( int_2_str(f) )//'.dat',&
          access='stream',&
          form='unformatted')

     read(10), seed
     read(10), tmp
     close(10)

     
     !--------------------
     !
     !--------------------
     points = reshape(tmp,[npoints,3],order=[2,1])
     
     !average end to end distance 
     avg_Ree = avg_Ree + norm2(  points(1,:) - points(npoints,:)  )/nreps
     !center of mass position
     Rcm = sum( points,dim=1 )/npoints 
     avg_Rg2 = avg_Rg2 + Rg2(points,Rcm) / nreps
     !contour length up to n-th point
     avg_L = avg_L + L(points) / nreps
     !distance travelled on the z-axis up to n-th point
     avg_Lz = avg_Lz + Lz(points) / nreps
     
  end do

  ! do i=1,npoints-1
  !    avg1_avgi(i) = dot_product( avg_t(1,:), avg_t(i,:) )
  ! end do
     

  open(newunit=file_results,file='/scratch/madorisi/results/polymers/results_'//trim( arg1 )//'.dat',status='replace')
  do i=1,npoints-1
    write(file_results,*), avg_L(i)/(nreps-c), avg_dz(i)/(nreps-c) !( avg_costheta(i) - avg1_avgi(i)/(nreps - c) ) / ( avg_costheta(1) - avg1_avgi(1)/(nreps-c) )
  end do
  close(file_results)

  !write(6,*), c, ' files with zero bytes ', arg1 !6=stdout
  

contains
  
  function int_2_str(int) result(str)
    !input a real. If integer cast it to real before
    implicit none
    
    integer,intent(in) :: int
    character(len=200) :: str
    
    write(str,'(I0)'), int
    
  end function int_2_str
  

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


  function Rg2(points,Rcm)
    implicit none

    real, intent(in) :: points(:,:), Rcm(3)
    real :: Rg2
    integer :: i
    
    Rg2 = 0.
    do i=1,size(points,1)
       Rg2 = Rg2 + dot_product(points(i,:),points(i,:))
    end do
    
    Rg2 = Rg2/size(points,1)
    Rg2 = Rg2 - dot_product(Rcm,Rcm)

  end function Rg2

  ! function contour(points,upto)
  !   implicit none

  !   real, intent(in) :: points(:,:)
  !   integer, intent(in) :: upto
  !   real :: contour
  !   integer :: i
    
  !   contour = 0.


  !   do i=1,upto
  !      contour = contour + norm2( points(i+1,:) - points(i,:) )
  !   end do
    
    


  ! end function contour
  function L(points)
    implicit none
    
    real, intent(in) :: points(:,:)
    real :: L(size(points,1)-1)
    
    L = 0.
    L(1) = norm2(points(2,:) - points(1,:))
    
    do i = 1,npoints-2
       L(i+1) = L(i) + norm2( points(i+2,:) - points(i+1,:) )
    end do
    
  end function L
  
  function Lz(points)
    implicit none
    !compute the jump along z 
    !Lz(i): jump on z up to the i-th jump

    real, intent(in) :: points(:,:)
    real :: Lz(size(points,1)-1)
    
    Lz = 0.
    Lz(1) = norm2(points(2,3) - points(1,3))
    
    do i = 1,npoints-2
       Lz(i+1) = Lz(i) + norm2( points(i+2,:) - points(i+1,:) )
    end do
    
  end function Lz
  
end program ana_dat
