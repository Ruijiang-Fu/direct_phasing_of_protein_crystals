!ifort -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lm
! 3D complex to complex
  include 'mkl_dfti.f90'

program main

  use MKL_DFTI

  implicit none
  real,parameter::pi=dacos(-1.0d+0)
  complex,parameter::IM=(0,1.0)
  character(len=4),parameter::PDB_ID='1UII'
  real,parameter::a=37.160d+0
  real,parameter::b=94.460d+0
  real,parameter::c=102.560d+0
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
  real,parameter::R_pseudo=0.0
  integer,parameter::anum=38
  integer,parameter::bnum=96
  integer,parameter::cnum=104
  integer,parameter::hnum=anum
  integer,parameter::knum=bnum
  integer,parameter::lnum=cnum
  real,parameter::forward_scale_factor=1.0/Vol
  real,parameter::backward_scale_factor=Vol/real(anum*bnum*cnum)
  real,parameter::cutoff_resolution=2.0
  integer,parameter::num_of_origin_choices=16
  integer i,j,ios
  integer h,k,l,obsstatus
  integer ia,ib,ic,ia0,ib0,ic0,da,db,dc
  integer ih,ik,il,ih0,ik0,il0
  real sign_of_phase,delta_phase
  integer num_of_grid_points_in_high_density_region
  integer irefl,num_of_uniq_refl
  integer iorigin
  real S,resolution
  real x,y,z,x0,y0,z0,cache(3)
  real density,density_min,density_ave,density_max
  real solvent_density_ave
  real density1,density2,cutoff_density
  real obssf,calsf,phase
  real imag_part,real_part
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
  real::density_on_grid_points_in_unitcell(0:anum-1,0:bnum-1,0:cnum-1)
  real::density_on_grid_points_in_unitcell_tmp(0:anum-1,0:bnum-1,0:cnum-1)
  real::list_of_uniq_refl(0:hnum/2,0:knum/2,0:lnum/2,3)


  length_of_each_dimension(1)=anum !==hnum
  length_of_each_dimension(2)=bnum !==knum
  length_of_each_dimension(3)=cnum !==lnum

  tab=char(9)

  do h=0,hnum/2,1
    do k=0,knum/2,1
      do l=0,lnum/2,1
        list_of_uniq_refl(h,k,l,1)=0  !obsstatus
        list_of_uniq_refl(h,k,l,2)=0  !obssf
        list_of_uniq_refl(h,k,l,3)=0  !phase
      end do
    end do
  end do

  !(1)preparing experimental observed structure factors
  write(filename,'(A4,A21)')PDB_ID,'_fmodel.hkl,.HKL,.sca'
  open(unit=10,file=filename,status='old')
  do
     read(10,'(A80)',iostat=ios)string
    if(ios/=0)exit
    if(string(1:6)==' 1 1 1')then
      read(string(10:50),'(3I5,3X,F11.3,F12.1)')h,k,l,obssf,phase
      obsstatus=1 !default
      phase=phase/180.0*pi
      phase=mod(phase,2.0*pi)
      list_of_uniq_refl(h,k,l,1)=obsstatus
      list_of_uniq_refl(h,k,l,2)=obssf !*exp(-(pi*R_pseudo*S)**2)
      list_of_uniq_refl(h,k,l,3)=phase
    end if
  end do
  close(10)

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

  do ia=0,anum-1,1
    do ib=0,bnum-1,1
      do ic=0,cnum-1,1
        density_on_grid_points_in_unitcell(ia,ib,ic)=real(real_space_3D(ia,ib,ic))
      end do
    end do
  end do


  do iorigin=1,num_of_origin_choices,1

    do ia=0,anum-1,1
      do ib=0,bnum-1,1
        do ic=0,cnum-1,1
          density_on_grid_points_in_unitcell_tmp(ia,ib,ic)=0
        end do
      end do
    end do

    do ia=0,anum-1,1
      do ib=0,bnum-1,1
        do ic=0,cnum-1,1
          x=real(ia+0.5)/real(anum)
          y=real(ib+0.5)/real(bnum)
          z=real(ic+0.5)/real(cnum)
          select case(iorigin)
            case(1)
              x=x
              y=y
              z=z
            case(2)
              x=x+0.5
              y=y
              z=z
            case(3)
              x=x 
              y=y+0.5
              z=z 
            case(4)
              x=x
              y=y
              z=z+0.5
            case(5)
              x=x+0.5
              y=y+0.5
              z=z
            case(6)
              x=x+0.5
              y=y
              z=z+0.5
            case(7)
              x=x
              y=y+0.5 
              z=z+0.5
            case(8)
              x=x+0.5
              y=y+0.5
              z=z+0.5
            case(9)
              x=-x+1.0
              y=-y+1.0
              z=-z+1.0
            case(10)
              x=-x+0.5
              y=-y+1.0
              z=-z+1.0
            case(11)
              x=-x+1.0
              y=-y+0.5
              z=-z+1.0
            case(12)
              x=-x+1.0
              y=-y+1.0
              z=-z+0.5
            case(13)
              x=-x+0.5
              y=-y+0.5
              z=-z+1.0
            case(14)
              x=-x+0.5
              y=-y+1.0
              z=-z+0.5
            case(15)
              x=-x+1.0 
              y=-y+0.5
              z=-z+0.5
            case(16)
              x=-x+0.5
              y=-y+0.5
              z=-z+0.5
          end select
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinates_into_unitcell(x,y,z)
          end if
          ia0=int(x*anum)
          ib0=int(y*bnum)
          ic0=int(z*cnum)
          density_on_grid_points_in_unitcell_tmp(ia0,ib0,ic0)=density_on_grid_points_in_unitcell(ia,ib,ic)
        end do
      end do
    end do

    do ia=0,anum-1,1
      do ib=0,bnum-1,1
        do ic=0,cnum-1,1
          real_space_3D(ia,ib,ic)=density_on_grid_points_in_unitcell_tmp(ia,ib,ic)
        end do
      end do
    end do

    !Perform a complex to complex transform
    Status = DftiCreateDescriptor( my_descriptor_handle, DFTI_SINGLE, DFTI_COMPLEX, 3, length_of_each_dimension )
    Status = DftiSetValue( my_descriptor_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE )
    Status = DftiSetValue( my_descriptor_handle, DFTI_BACKWARD_SCALE, backward_scale_factor )
    Status = DftiCommitDescriptor( my_descriptor_handle )
    Status = DftiComputeBackward( my_descriptor_handle, real_space, reciprocal_space )
    Status = DftiFreeDescriptor( my_descriptor_handle )

    do ih=0,hnum-1,1
      do ik=0,knum-1,1
        do il=0,lnum-1,1
          if(ih<=(hnum/2-1))then
            h=ih
          else if(ih>(hnum/2-1))then
            h=ih-hnum
          end if
          if(ik<=(knum/2-1))then
            k=ik
          else if(ik>(knum/2-1))then
            k=ik-knum
          end if
          if(il<=(lnum/2-1))then
            l=il
          else if(il>(lnum/2-1))then
            l=il-lnum
          end if
          reciprocal_space_3D(ih,ik,il)=reciprocal_space_3D(ih,ik,il)  &
            *exp( 2.0*pi*IM*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum)) ) 
        end do
      end do
    end do

    do ih=0,hnum-1,1
      do ik=0,knum-1,1
        do il=0,lnum-1,1
          if(ih<=(hnum/2-1))then
            h=ih
          else if(ih>(hnum/2-1))then
            h=ih-hnum
          end if
          if(ik<=(knum/2-1))then
            k=ik
          else if(ik>(knum/2-1))then
            k=ik-knum
          end if
          if(il<=(lnum/2-1))then
            l=il
          else if(il>(lnum/2-1))then
            l=il-lnum
          end if
          if(h<0 .or. k<0 .or. l<0)cycle  !not unique reflections
           S=sqrt( h**2*S1+k**2*S2+l**2*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6 )/Vol
          resolution=1.0/S 
          if(resolution<cutoff_resolution)cycle !!!!
          real_part=real(reciprocal_space_3D(ih,ik,il))
          imag_part=imag(reciprocal_space_3D(ih,ik,il))
          calsf=sqrt(real_part**2+imag_part**2)
          if(real_part>0 .and. imag_part>0)then  
            phase=atan(imag_part/real_part)
          else if(real_part<0 .and. imag_part>0)then  
            phase=atan(imag_part/real_part)+pi
          else if(real_part<0 .and. imag_part<0)then  
            phase=atan(imag_part/real_part)+pi
          else if(real_part>0 .and. imag_part<0)then  
            phase=atan(imag_part/real_part)+2.0*pi
          else if(real_part>=0 .and. imag_part==0)then
            phase=0
          else if(real_part==0 .and. imag_part>0)then     
            phase=0.50*pi
          else if(real_part<0 .and. imag_part==0)then 
            phase=pi
          else if(real_part==0 .and. imag_part<0)then   
            phase=1.50*pi
          end if
          if(abs(calsf)<0.1)phase=0  !attention, not good to find phase
          phase=mod(phase,2.0*pi)
          if(phase<0) phase=phase+2.0*pi
          list_of_uniq_refl(h,k,l,3)=phase
        end do
      end do
    end do
  
    if(iorigin<10)then
      write(filename,'(A30,I1,A4)')'true_phases_for_origin_choice_',iorigin,'.txt'
    else if(iorigin<100)then
      write(filename,'(A30,I2,A4)')'true_phases_for_origin_choice_',iorigin,'.txt'
    end if
    
    open(unit=11,file=filename,status='replace')
    do h=0,hnum/2,1
      do k=0,knum/2,1
        do l=0,lnum/2,1
          S=sqrt(h**2*S1+k**2*S2+l**2*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6)/Vol
          resolution=1.0/S
          if(resolution>=cutoff_resolution)then
            obsstatus=nint(list_of_uniq_refl(h,k,l,1))
            obssf=list_of_uniq_refl(h,k,l,2)
            phase=list_of_uniq_refl(h,k,l,3)
            write(11,'(3I5,I2,F12.3,F8.3)')h,k,l,obsstatus,obssf,phase/pi*180.
          end if
        end do
      end do
    end do  
    close(11)



  end do !iorigin



  stop
end


!----------------------------------------------------------------------------80

subroutine initializ_random_seed()

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

end subroutine initializ_random_seed


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

subroutine convert_Cartesian_coordinates_to_fractional_coordinates(a,b,c,alpha,beta,gamm,x,y,z)

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
  
  !convert fractional coordinates to Cartesian coordinates
  fract=MATMUL(matrix,cart)

  x=fract(1)
  y=fract(2)
  z=fract(3)

  return
end

!----------------------------------------------------------------------------80

subroutine convert_fractional_coordinates_to_Cartesian_coordinates(a,b,c,alpha,beta,gamm,x,y,z)

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
  
  !convert fractional coordinates to Cartesian coordinates
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

  if(h==0 .and.k==0 .and. l==0)then
    sign_of_phase=1
    delta_phase=0
  end if
 
  if(h==0 .and.k==0 .and. l>0)then
    sign_of_phase=1
    delta_phase=0
  else if(h==0 .and.k==0 .and. l<0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(l)),2)
  end if

  if(h==0 .and.k>0 .and. l==0)then
    sign_of_phase=1
    delta_phase=0
  else if(h==0 .and.k<0 .and. l==0)then
    sign_of_phase=-1
    delta_phase=mod((abs(k)+abs(l)),2)
  end if

  if(h>0 .and.k==0 .and. l==0)then
    sign_of_phase=1
    delta_phase=0
  else if(h<0 .and.k==0 .and. l==0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(k)),2)
  end if

  if(h==0 .and.k>0 .and. l>0)then
    sign_of_phase=1
    delta_phase=0
  else if(h==0 .and.k>0 .and. l<0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(l)),2)
  else if(h==0 .and.k<0 .and. l>0)then
    sign_of_phase=-1
    delta_phase=mod((abs(k)+abs(l)),2)
  else if(h==0 .and.k<0 .and. l<0)then
    sign_of_phase=1
    delta_phase=mod(-(abs(h)+abs(k)),2)
  end if

  if(h>0 .and.k==0 .and. l>0)then
    sign_of_phase=1
    delta_phase=0
  else if(h>0 .and.k==0 .and. l<0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(l)),2)
  else if(h<0 .and.k==0 .and. l>0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(k)),2)
  else if(h<0 .and.k==0 .and. l<0)then
    sign_of_phase=1
    delta_phase=mod(-(abs(k)+abs(l)),2)
  end if

  if(h>0 .and. k>0 .and. l==0)then
    sign_of_phase=1
    delta_phase=0
  else if(h>0 .and. k<0 .and. l==0)then
    sign_of_phase=-1
    delta_phase=mod((abs(k)+abs(l)),2)
  else if(h<0 .and. k>0 .and. l==0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(k)),2)
  else if(h<0 .and. k<0 .and. l==0)then
    sign_of_phase=1
    delta_phase=mod(-(abs(h)+abs(l)),2)
  end if

  if(h>0 .and. k>0 .and. l>0)then
    sign_of_phase=1
    delta_phase=0
  else if(h>0 .and. k>0 .and. l<0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(l)),2)
  else if(h>0 .and. k<0 .and. l>0)then
    sign_of_phase=-1
    delta_phase=mod((abs(k)+abs(l)),2)
  else if(h>0 .and. k<0 .and. l<0)then
    sign_of_phase=1
    delta_phase=mod(-(abs(h)+abs(k)),2)
  else if(h<0 .and. k>0 .and. l>0)then
    sign_of_phase=-1
    delta_phase=mod((abs(h)+abs(k)),2)
  else if(h<0 .and. k>0 .and. l<0)then 
    sign_of_phase=1
    delta_phase=mod(-(abs(k)+abs(l)),2)
  else if(h<0 .and. k<0 .and. l>0)then
    sign_of_phase=1
    delta_phase=mod(-(abs(h)+abs(l)),2)
  else if(h<0 .and. k<0 .and. l<0)then
    sign_of_phase=-1
    delta_phase=0
  end if

  return
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
