!ifort -o test test.f90 -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lm
! 3D complex to complex
  include 'mkl_dfti.f90'

program main

  use MKL_DFTI

  implicit none
  real,parameter::pi=dacos(-1.0d+0)
  complex,parameter::IM=(0,1.0)
  character(len=4),parameter::PDB_ID='1Y5Y'
  real,parameter::a=120.659d+0
  real,parameter::b=120.659d+0
  real,parameter::c=205.947d+0
  real,parameter::alpha=90.00d+0/180.0*pi
  real,parameter::beta=90.00d+0/180.0*pi
  real,parameter::gamm=90.00d+0/180.0*pi
  real,parameter::S1=(b*c*sin(alpha))**2
  real,parameter::S2=(a*c*sin(beta))**2
  real,parameter::S3=(a*b*sin(gamm))**2
  real,parameter::S4=a*b*c**2*(cos(alpha)*cos(beta)-cos(gamm))
  real,parameter::S5=a*b**2*c*(cos(gamm)*cos(alpha)-cos(beta))
  real,parameter::S6=a**2*b*c*(cos(beta)*cos(gamm)-cos(alpha))
  real,parameter::nu=sqrt( 1.0-cos(alpha)**2-cos(beta)**2-cos(gamm)**2   &
                             +2.0*cos(alpha)*cos(beta)*cos(gamm) )
  real,parameter::Vol=a*b*c*nu
  real,parameter::sigma_weight=1.50  !A
  integer,parameter::anum=120
  integer,parameter::bnum=120
  integer,parameter::cnum=208 !mod(cnum,4)=0  for P43212 space group
  integer,parameter::hnum=anum
  integer,parameter::knum=bnum
  integer,parameter::lnum=cnum
  real,parameter::forward_scale_factor=1.0/Vol
  real,parameter::backward_scale_factor=Vol/real(anum*bnum*cnum)
  real,parameter::solvent_content=0.60 !68%
  real,parameter::cutoff_resolution=2.250 !same resolution as 2UXJ
  integer i,j,ios
  integer h,k,l,obsstatus
  integer ia,ib,ic,ia0,ib0,ic0,da,db,dc
  integer ih,ik,il,ih0,ik0,il0
  real sign_of_phase,delta_phase
  integer num_of_grid_points_in_high_density_region
  integer irefl,num_of_uniq_refl
  real S,resolution
  real x,y,z,x0,y0,z0,cache(3)
  real density,density_min,density_ave,density_max
  real density1,density2,cutoff_density
  real obssf,calsf,phase
  character(len=80) filename,string
  character tab

  integer::Status
  type(dfti_descriptor),pointer::my_descriptor_handle
  complex::real_space_3D(0:anum-1,0:bnum-1,0:cnum-1)
  complex::real_space(anum*bnum*cnum)
  equivalence(real_space_3D, real_space)
  complex::reciprocal_space_3D(0:hnum-1,0:knum-1,0:lnum-1)
  complex::reciprocal_space(hnum*knum*lnum)
  equivalence(reciprocal_space_3D,reciprocal_space)
  integer::length_of_each_dimension(3)

  complex::complex_structure_factors(-hnum/2:hnum/2-1,-knum/2:knum/2-1,-lnum/2:lnum/2-1)
  real::density_on_grid_points_in_ASU(0:anum/2-1,0:bnum/2-1,0:cnum/2-1)
  real::list_of_uniq_refl(0:hnum/2,0:knum/2,0:lnum/2,3)

  write(*,*)
  write(*,*)'Cutoff resolution=',cutoff_resolution
  write(*,*)'Assume solvent content=',solvent_content


  length_of_each_dimension(1)=anum !==hnum
  length_of_each_dimension(2)=bnum !==knum
  length_of_each_dimension(3)=cnum !==lnum

  tab=char(9)

  do h=0,hnum/2,1
    do k=0,knum/2,1
      do l=0,lnum/2,1
        list_of_uniq_refl(h,k,l,1)=0
        list_of_uniq_refl(h,k,l,2)=0
        list_of_uniq_refl(h,k,l,3)=0
      end do
    end do
  end do

  !(1)preparing experimental observed structure factors
  open(unit=10,file='syn_F_atoms.txt',status='old')
  do
    read(10,*,iostat=ios)h,k,l,obssf,phase
    if(ios/=0)exit
    S=sqrt( h**2*S1+k**2*S2+l**2*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6 )/Vol
    resolution=1.0/S
    phase=phase/180.0*pi
    if(resolution>=cutoff_resolution)then
      list_of_uniq_refl(h,k,l,1)=obsstatus
      list_of_uniq_refl(h,k,l,2)=obssf*sqrt(2.0*pi*sigma_weight**2)*exp(-2.0*(pi*sigma_weight*S)**2)
      list_of_uniq_refl(h,k,l,3)=phase
      if(h>k)then
        list_of_uniq_refl(k,h,l,1)=obsstatus
        list_of_uniq_refl(k,h,l,2)=obssf*sqrt(2.0*pi*sigma_weight**2)*exp(-2.0*(pi*sigma_weight*S)**2)
        phase=mod(l,2)*pi-phase
        phase=mod(phase,2.0*pi)
        if(phase<0)phase=phase+2.0*pi
        list_of_uniq_refl(k,h,l,3)=phase
      end if
    end if
  end do
  close(10)

  write(*,*)'F(0,0,0)=',list_of_uniq_refl(0,0,0,2)
  write(*,*)

  !(3)preparing complex_structure_factors(-h:h,-k:k,-l:l)

  do h=-hnum/2,hnum/2-1,1
    do k=-knum/2,knum/2-1,1  
      do l=-lnum/2,lnum/2-1,1
        obssf=list_of_uniq_refl(abs(h),abs(k),abs(l),2) !in case of too small grid resulting in big h,k,l
        phase=list_of_uniq_refl(abs(h),abs(k),abs(l),3)
        call expand_reflections(h,k,l,sign_of_phase,delta_phase)
        phase=sign_of_phase*phase+delta_phase*pi
        phase=mod(phase,2.0*pi)
        if(phase<0) phase=phase+2.0*pi
        complex_structure_factors(h,k,l)=obssf*(cos(phase)+IM*sin(phase))
      end do
    end do
  end do

  !(4)preparing DFTI FFT
  !shift half grid in order to make density on the center of each grid for plotting with symmetry
  !here (ih,ik,il) mean address, (h,k,l) are Miller indices
  do ih=0,hnum-1,1
    do ik=0,knum-1,1
      do il=0,lnum-1,1
        call convert_index(ih,ik,il,h,k,l,hnum,knum,lnum)
        reciprocal_space_3D(ih,ik,il)=complex_structure_factors(h,k,l)         &
          *exp(-2.0*pi*IM*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum)) )
      end do
    end do  
  end do

  !(5)preparing initial values
  !Perform a complex to complex transform
  Status = DftiCreateDescriptor( my_descriptor_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, length_of_each_dimension )
  Status = DftiSetValue( my_descriptor_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
  Status = DftiSetValue( my_descriptor_handle, DFTI_FORWARD_SCALE, forward_scale_factor )
  Status = DftiCommitDescriptor( my_descriptor_handle )
  Status = DftiComputeForward( my_descriptor_handle, reciprocal_space, real_space )
  Status = DftiFreeDescriptor( my_descriptor_handle )

  do ia=0,anum/2-1,1
    do ib=0,bnum/2-1,1
      do ic=0,cnum/2-1,1
        density_on_grid_points_in_ASU(ia,ib,ic)=real(real_space_3D(ia,ib,ic))
      end do
    end do
  end do

  density_min=huge(density_min)
  density_ave=0
  density_max=-huge(density_max)
  do ia=0,anum/2-1,1
    do ib=0,bnum/2-1,1
      do ic=0,cnum/2-1,1   
        density=density_on_grid_points_in_ASU(ia,ib,ic)
        if(density<density_min) density_min=density
        if(density>density_max) density_max=density
        density_ave=density_ave+density
      end do
    end do
  end do
  density_ave=density_ave/real(anum/2*bnum/2*cnum/2)
  write(*,*)'density_min_in_whole_unitcell=',density_min
  write(*,*)'density_ave_in_whole_unitcell=',density_ave
  write(*,*)'density_max_in_whole_unitcell=',density_max
  write(*,*)

  num_of_grid_points_in_high_density_region=int(anum/2*bnum/2*cnum/2*(1.0-solvent_content))  !*****ATTENTION****

  density1=density_min
  density2=density_max
  do
    cutoff_density=(density2+density1)/2.0
    i=0
    do ia=0,anum/2-1,1
      do ib=0,bnum/2-1,1
        do ic=0,cnum/2-1,1
          if(density_on_grid_points_in_ASU(ia,ib,ic)>=cutoff_density)then
            i=i+1
          end if
        end do
      end do
    end do
    if(i>num_of_grid_points_in_high_density_region)then
      density1=cutoff_density
      if(abs(density2-density1)<10.0**(-5))exit
    else if(i<num_of_grid_points_in_high_density_region)then
      density2=cutoff_density
      if(abs(density2-density1)<10.0**(-5))exit
    else if(i==num_of_grid_points_in_high_density_region)then
      exit
    end if
  end do
  num_of_grid_points_in_high_density_region=i

  density_min=huge(density_min)
  density_ave=0
  density_max=-huge(density_max)
  do ia=0,anum/2-1,1
    do ib=0,bnum/2-1,1
      do ic=0,cnum/2-1,1
        density=density_on_grid_points_in_ASU(ia,ib,ic)
        if(density>=cutoff_density)then
          if(density<density_min) density_min=density
          if(density>density_max) density_max=density
          density_ave=density_ave+density
        end if
      end do
    end do
  end do
  density_ave=density_ave/real(num_of_grid_points_in_high_density_region)
  write(*,*)'density_min_in_molecular_region=',density_min
  write(*,*)'density_ave_in_molecular_region=',density_ave
  write(*,*)'density_max_in_molecular_region=',density_max
  write(*,*)

  density_min=huge(density_min)
  density_ave=0
  density_max=-huge(density_max)
  i=0
  do ia=0,anum/2-1,1
    do ib=0,bnum/2-1,1
      do ic=0,cnum/2-1,1
        density=density_on_grid_points_in_ASU(ia,ib,ic)
        if(density<cutoff_density)then
          i=i+1
          if(density<density_min) density_min=density
          if(density>density_max) density_max=density
          density_ave=density_ave+density
        end if
      end do
    end do
  end do
  density_ave=density_ave/real(i)
  write(*,*)'density_min_in_solvent_region=',density_min
  write(*,*)'density_ave_in_solvent_region=',density_ave
  write(*,*)'density_max_in_solvent_region=',density_max
  write(*,*)

  open(unit=12,file='boundary.pdb',status='replace')
  write(12,'(A6,3F9.3,3F7.2,A11,4X,A2,9X)')  &
    'CRYST1',a,b,c,alpha/pi*180.0,beta/pi*180.0,gamm/pi*180.0,' P 43 21 2 ','40'
  i=0
  do ia=0,anum/2-1,1
    do ib=0,bnum/2-1,1
      do ic=0,cnum/2-1,1
        if(density_on_grid_points_in_ASU(ia,ib,ic)>=cutoff_density)then
          cache(1)=real(ia+0.5)/real(anum)
          cache(2)=real(ib+0.5)/real(bnum)
          cache(3)=real(ic+0.5)/real(cnum)
          i=i+1
          x=cache(1)
          y=cache(2)
          z=cache(3)
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z 
          x=-cache(1)+1.0
          y=-cache(2)+1.0
          z= cache(3)+0.5
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          i=i+1
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
          x=-cache(2)+0.5
          y= cache(1)+0.5
          z= cache(3)+0.75
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          i=i+1
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
          x= cache(2)+0.5
          y=-cache(1)+0.5
          z= cache(3)+0.25
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          i=i+1
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
          x=-cache(1)+0.5
          y= cache(2)+0.5
          z=-cache(3)+0.75
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          i=i+1
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
          x= cache(1)+0.5
          y=-cache(2)+0.5
          z=-cache(3)+0.25
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          i=i+1
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
          x= cache(2)
          y= cache(1)
          z=-cache(3)+1.0
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          i=i+1
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
          x=-cache(2)+1.0
          y=-cache(1)+1.0
          z=-cache(3)+0.5
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          call convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)
          i=i+1
          write(12,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
        end if
      end do  
    end do
  end do
  close(12)


  stop
end


!----------------------------------------------------------------------------80

subroutine initialize_random_seed()

  implicit none

  integer :: i
  integer :: n
  integer :: clock
  integer, dimension(:), allocatable :: seed
          
  call random_seed(size = n)

  allocate(seed(n))
   
  call system_clock(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(PUT = seed)
         
  deallocate(seed)

  return

end subroutine initialize_random_seed


!-----------------------------------------------------------------------------80

subroutine adjust_fractional_coordinates_into_unitcell(x,y,z)

  implicit none
  integer i
  real x,y,z
  real fract(3)

  if(x>=0 .and. x<=1.0 .and. y>=0 .and. y<=1.0 .and. z>=0 .and. z<=1.0)return

  fract(1)=x
  fract(2)=y
  fract(3)=z

  !make tranlations
  do i=1,3,1
    if(fract(i)<-2.0)then
      fract(i)=fract(i)+3.0
    else if(fract(i)<-1.0)then
      fract(i)=fract(i)+2.0
    else if(fract(i)<0)then
      fract(i)=fract(i)+1.0
    else if(fract(i)>2.0)then
      fract(i)=fract(i)-2.0
    else if(fract(i)>1.0)then
      fract(i)=fract(i)-1.0
    end if 
  end do

  x=fract(1)
  y=fract(2)
  z=fract(3)

  return
end

!-----------------------------------------------------------------------------80

subroutine convert_orthogonal_coordinates_to_fractional_coordinates(a,b,c,alpha,beta,gamm,x,y,z)

  implicit none
  real,parameter::pi=dacos(-1.0d+0)
  real nu
  real a,b,c,alpha,beta,gamm,x,y,z
  real cart(3),fract(3)
  real matrix(3,3)

  nu=sqrt( 1.0-cos(alpha)**2.0-cos(beta)**2.0-cos(gamm)**2.0 &
    +2.0*cos(alpha)*cos(beta)*cos(gamm) )

  matrix(1,1)=1.0/a
  matrix(1,2)=-cos(gamm)/(a*sin(gamm))
  matrix(1,3)=(cos(alpha)*cos(gamm)-cos(beta))/(a*sin(gamm)*nu)
  matrix(2,1)=0
  matrix(2,2)=1.0/(b*sin(gamm))
  matrix(2,3)=-(cos(alpha)-cos(beta)*cos(gamm))/(b*sin(gamm)*nu)
  matrix(3,1)=0
  matrix(3,2)=0
  matrix(3,3)=sin(gamm)/(c*nu)

  cart(1)=x
  cart(2)=y
  cart(3)=z
  
  !convert fractional coordinates to orthogonal coordinates
  fract=MATMUL(matrix,cart)

  x=fract(1)
  y=fract(2)
  z=fract(3)

  return
end

!----------------------------------------------------------------------------80

subroutine convert_fractional_coordinates_to_orthogonal_coordinates(a,b,c,alpha,beta,gamm,x,y,z)

  implicit none
  real,parameter::pi=dacos(-1.0d+0)
  real nu
  real a,b,c,alpha,beta,gamm,x,y,z
  real cart(3),fract(3)
  real matrix(3,3)

  nu=sqrt( 1.0-cos(alpha)**2-cos(beta)**2-cos(gamm)**2 &
    +2.0*cos(alpha)*cos(beta)*cos(gamm) )

  matrix(1,1)=a
  matrix(1,2)=b*cos(gamm)
  matrix(1,3)=c*cos(beta)
  matrix(2,1)=0
  matrix(2,2)=b*sin(gamm)
  matrix(2,3)=c*( cos(alpha)-cos(beta)*cos(gamm) )/sin(gamm)
  matrix(3,1)=0
  matrix(3,2)=0
  matrix(3,3)=c*nu/sin(gamm)

  fract(1)=x
  fract(2)=y
  fract(3)=z
  
  !convert fractional coordinates to orthogonal coordinates
  cart=MATMUL(matrix,fract)

  x=cart(1)
  y=cart(2)
  z=cart(3)

  return
end

!-----------------------------------------------------------------------------80

subroutine expand_reflections(h,k,l,sign_of_phase,delta_phase)

  implicit none
  integer h,k,l
  real sign_of_phase,delta_phase

!  input h,k,l
!  output sign_of_phase, delta_phase 
!  phase:  sign_of_phase*phase+delta_phase*pi

  if(h>=0 .and. k>=0 .and. l>=0)then
    sign_of_phase=1
    delta_phase=0
    return
  end if
 
  if(h>=0 .and. k>=0 .and. l<=0)then
    sign_of_phase=-1
    delta_phase=mod(abs(real(l)),2.0)
    return
  end if
 
  if(h>=0 .and. k<=0 .and. l>=0)then
    sign_of_phase=-1
    delta_phase=mod(abs(h)+abs(k)+1.5*abs(l),2.0)
    return
  end if
 
  if(h<=0 .and. k>=0 .and. l>=0)then
    sign_of_phase=-1
    delta_phase=mod(abs(h)+abs(k)+0.5*abs(l),2.0)
    return
  end if

  if(h>=0 .and. k<=0 .and. l<=0)then
    sign_of_phase=1
    delta_phase=mod(abs(h)+abs(k)-0.5*abs(l),2.0)
    return
  end if

  if(h<=0 .and. k>=0 .and. l<=0)then
    sign_of_phase=1
    delta_phase=mod(abs(h)+abs(k)-1.5*abs(l),2.0)
    return
  end if

  if(h<=0 .and. k<=0 .and. l>=0)then
    sign_of_phase=1
    delta_phase=mod(abs(real(l)),2.0)
    return
  end if

  if(h<=0 .and. k<=0 .and. l<=0)then
    sign_of_phase=-1
    delta_phase=0
    return
  end if

end

!-----------------------------------------------------------------------------80

subroutine convert_index(ih,ik,il,h,k,l,hnum,knum,lnum)

  implicit none
  integer ih,ik,il,h,k,l,hnum,knum,lnum

  !input 0<=ih<=hnum-1, 0<=ik<=knum-1, 0<=il<=lnum-1
  !output h,k,l

  if(ih<=(hnum/2-1)) h=ih
  if(ih>(hnum/2-1)) h=ih-hnum

  if(ik<=(knum/2-1)) k=ik
  if(ik>(knum/2-1)) k=ik-knum

  if(il<=(lnum/2-1)) l=il
  if(il>(lnum/2-1)) l=il-lnum

  return
end

!-----------------------------------------------------------------------------80
