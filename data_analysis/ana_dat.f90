program ana_dat

  implicit none

  integer, parameter :: nreps = 10000
  integer, parameter :: npoints = 10000
  integer :: seed
  integer(4) :: values(13),file_results
  integer(4) :: f,i,blocksize,fsize,c
  logical :: exist
  real(4) :: tmp(3*npoints)
  real(4) :: points(npoints,3)
  real(4) :: t(npoints-1,3),avg_t(npoints-1,3)
  real(4) :: avgjmp, normt,avg_Ree,Ree,Rcm(3),avg_Rg2
  character(100) :: arg1,main_input_dir
  real(4) :: avg_L(npoints-1),& !cumulative total length
             avg_Lz(npoints-1),& !cumulative jump on z
             avg_Ree(npoints-1),& !end-to-end up to the n-th point
             avg_Rg2(npoints-1),& !gyradius up to the n-th point
             avg_dz(npoints-1),&  !z jump
             avg_dl(npoints-1),& !jump in space

  call getarg(1,arg1)
  read(arg1,*), avgjmp
  main_input_dir = '/media/USB-HDD/data_cylinder/avjmp_'//trim(arg1)//'/'
  !main_input_dir = '/scratch/madorisi/results/polymers/avjmp_'//trim(arg1)//'/'
  !main_output_dir = '../output/'
  print*, 'working on ',trim( main_input_dir )
  

  avg_Ree = 0. !array
  avg_Rg2 = 0. !array
  avg_L = 0. ! array
  avg_Lz = 0. !array
  avg_dz = 0.
  avg_dl = 0.
  tot_contour = 0.
  c = 0
  
  
  do f = 0,nreps-1
     print*, f
     !check file
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
     
     ! do i=1,npoints-1
     !    tot_contour = tot_contour + norm2(points(i+1,:) - points(i,:))
     ! end do

     ! avg_L = L(points)
     ! do i =1,size(avg_L,1)
     !    print*, avg_L(i)
     ! end do
     ! print*, tot_contour, avg_L(npoints-1)
     ! stop

     
     !avg_Ree = avg_Ree + Ree(points)!norm2(  points(1,:) - points(npoints,:)  )/nreps
     !center of mass position
     !Rcm = sum( points,dim=1 )/npoints 
     
     !average end to end distance 
     
     !gyration radius
     avg_Rg2 = avg_Rg2 + Rg2(points,Rcm) / nreps
     !contour length up to n-th point
     avg_L = avg_L + L(points) / nreps
     !distance travelled on the z-axis up to n-th point
     avg_Lz = avg_Lz + Lz(points) / nreps
     avg_dz = avg_dz + dz(points) / nreps
     avg_dl = avg_dl + dl(points) / nreps
     
  end do

  ! do i=1,npoints-1
  !    avg1_avgi(i) = dot_product( avg_t(1,:), avg_t(i,:) )
  ! end do
  call write_results(arg1,npoints,avg_L,avg_Lz,avg_dl,avg_dz,avg_Rg2,avg_Ree)

  
  print*, trim( main_input_dir ),' done'

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
    Lz(1) = points(2,3) - points(1,3)
    
    do i = 1,npoints-2
       Lz(i+1) = Lz(i) + points(i+2,3) - points(i+1,3) 
    end do
    
  end function Lz

  function nth_Ree(points)
    implicit none
    !compute the jump along z 
    !Lz(i): jump on z up to the i-th jump

    real, intent(in) :: points(:,:)
    real :: Ree(size(points,1)-1)
    
    Ree = 0.
        
    do i = 1,npoints-1
       Ree(i) = norm2(  points(i+1,:) - points(1,:)  )
    end do
    
  end function nth_Ree

  function nth_Rg2(points)
    implicit none

    real, intent(in) :: points(:,:)
    real :: Rg2(npoints-1),Rcm(3)
    integer :: np,k
    

    do np=2,size(points,1)
       
       Rcm = sum( points(1:np,:),dim=1 ) / real(np)
       Rg2 = 0.
       do k=1,np
          Rg2(np-1) = Rg2(np-1) + dot_product(points(k,:),points(k,:))
       end do

       Rg2(np-1) = Rg2(np-1)/real(np)
       Rg2(np-1) = Rg2(np-1) - dot_product(Rcm,Rcm)
       
    end do
    

  end function Rg2

  
  function dz(points)
    implicit none
    !compute the jump along z 
    !Lz(i): jump on z up to the i-th jump

    real, intent(in) :: points(:,:)
    real :: dz(size(points,1)-1)
    
    dz = 0.
    
    do i = 1,npoints-1
       dz(i) = points(i+1,3) - points(i,3) 
    end do
    
  end function dz
  
  
  function dl(points)
    implicit none
    !compute the jump along z 
    !Lz(i): jump on z up to the i-th jump

    real, intent(in) :: points(:,:)
    real :: dl(size(points,1)-1)
    
    dl = 0.
    
    do i = 1,npoints-1
       dl(i) = norm2( points(i+1,:) - points(i,:) ) 
    end do
    
  end function dl
  

  subroutine write_results(label,npoints,avg_L,avg_Lz,avg_dl,avg_dz,avg_Rg2,avg_Ree)
    
    implicit none

    integer, intent(in) :: npoints 
    character(len=100) :: label
    real(4), intent(in) :: avg_L(npoints-1), avg_Lz(npoints-1), avg_dl(npoints-1), avg_dz(npoints-1)
    real(4), intent(in) :: avg_Ree,avg_Rg2
    real(4) :: avgjmp
    integer :: i, file_results


    !read(arg1,*), avgjmp


    open(newunit=file_results,file='/scratch/madorisi/results/polymers/results_L_Lz_'//trim( label )//'.dat',status='unknown')
    write(file_results,*), 'avg_L ', 'avg_Lz ', 'avg_dl ', 'avg_dz ', 'avgjmp'
    do i=1,npoints-1
       write(file_results,*), avg_L(i), avg_Lz(i), avg_dl(i), avg_dz(i), trim(label)!avgjmp
    end do
    close(file_results)


    inquire (file='./results/polymers/results_Ree_Rg.dat', exist=exist)
    if (exist) then
       open(newunit=file_results, file='/scratch/madorisi/results/polymers/results_Ree_Rg.dat',&
            status="old", position="append", action="write")
       write(file_results,*), avg_Rg2, avg_Ree, avg_L(npoints-1), trim(label)!avgjmp
       close(file_results)
    else
       !create file and write the first time
       open(newunit=file_results, file='./results/polymers/results_Ree_Rg.dat', status="new", action="write")
       write(file_results,*), '# avgRg2 ','avgRee ','avgL ','avgjmp'
       write(file_results,*), avg_Rg2, avg_Ree, avg_L(npoints-1), trim(label)!avgjmp
       close(file_results)
    end if




  end subroutine write_results


end program ana_dat
