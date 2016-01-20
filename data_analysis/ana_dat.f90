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
  real(4) :: avgjmp, normt
  character(100) :: arg1,main_input_dir
  real(8) :: avg_l(npoints-1), avg_costheta(npoints-1), avg_dz(npoints-1), avg1_avgi(npoints-1)

  call getarg(1,arg1)
  read(arg1,*), avgjmp

  main_input_dir = '/scratch/madorisi/results/polymers/avjmp_'//trim(arg1)//'/'
  !main_output_dir = '../output/'
  print*, 'working on ',trim( main_input_dir )
  avg_l = 0.
  c = 0
  do f = 0,nreps-1
     !print*, f
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

     
     
     points = reshape(tmp,[npoints,3],order=[2,1])
     
     
     
     do i = 1,npoints-1
        
        normt = norm2( points(i+1,:) - points(i,:) )
        !t(i,:) = ( points(i+1,:) - points(i,:) )/normt     !versor between two consecutive points
        avg_l(i) = avg_l(i) + normt/(nreps-1)
        avg_dz(i) = avg_dz(i) + ( points(i+1,3) - points(i,3) )/(nreps-1) !jump on z direction
     
     end do
     
     ! do i=1,npoints-1
     !    avg_t(i,:) = avg_t(i,:) + t(i,:)   
     ! end do
     
     
     !avg_costheta = avg_costheta + costheta(t,npoints) 
     
  end do

  ! do i=1,npoints-1
  !    avg1_avgi(i) = dot_product( avg_t(1,:), avg_t(i,:) )
  ! end do
     

  open(newunit=file_results,file='/scratch/madorisi/results/polymers/results_'//trim( arg1 )//'.dat',status='replace')
  do i=1,npoints-1
    write(file_results,*), avg_l(i)/(nreps-c), avg_dz(i)/(nreps-c) !( avg_costheta(i) - avg1_avgi(i)/(nreps - c) ) / ( avg_costheta(1) - avg1_avgi(1)/(nreps-c) )
  end do
  close(file_results)

  write(6,*), c, ' files with zero bytes ', arg1 !6=stdout
  

contains
  
  function int_2_str(int) result(str)
    !input a real. If integer cast it to real before
    implicit none
    
    integer,intent(in) :: int
    character(len=200) :: str
    
    write(str,'(I0)'), int
    
  end function int_2_str
  

  function costheta(t,npoints)
    implicit none

    real,intent(in) :: t(:,:)
    integer, intent(in) :: npoints
    real :: costheta(npoints-1)
    
    integer :: j
        
    do i = 1,npoints-1
       costheta(i) = dot_product( t(1,:) , t(i,:) )  
    end do
        
  end function costheta


end program ana_dat
