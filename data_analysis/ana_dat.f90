program ana_dat

  implicit none

  integer, parameter :: nreps = 10000
  integer, parameter :: npoints = 10000
  integer :: seed
  integer :: f,i
  real(4) :: tmp(3*npoints)
  real(4) :: points(npoints,3)
  real(4) :: avgjmp
  character(100) :: arg1,main_input_dir
  real(4) :: avg_l(npoints-1)

  call getarg(1,arg1)
  read(arg1,*), avgjmp


  main_input_dir = '/media/usb/data_cylinder/avjmp_'//trim(arg1)//'/'
  !main_output_dir = '../output/'


  avg_l = 0.
  
  do f = 1,nreps
     
     !print*, trim(main_output_dir)//'rep_'//trim( int_2_str(f) )//'.dat'
     open(10,file=trim(main_input_dir)//'rep_'//trim( int_2_str(f) )//'.dat',&
          access='stream',&
          form='unformatted')

     read(10), seed
     read(10), tmp
     close(10)

     points = reshape(tmp,[npoints,3],order=[2,1])

     do i =1,npoints-1
        avg_l(i) = avg_l(i) + norm2( points(i+1,:) - points(i,:) ) 
     end do
  end do

  do i=1,npoints-1
     write(10,*), avg_l(i)/nreps
  end do


contains
  
  function int_2_str(int) result(str)
    !input a real. If integer cast it to real before
    implicit none
    
    integer,intent(in) :: int
    character(len=200) :: str
    
    write(str,'(I0)'), int
    
  end function int_2_str
  
end program ana_dat
