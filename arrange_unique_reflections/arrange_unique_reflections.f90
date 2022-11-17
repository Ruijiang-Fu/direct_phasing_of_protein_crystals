program main
  implicit none
  real,parameter::pi=dacos(-1.0d+0)
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
  real,parameter::nu=sqrt( 1.0-cos(alpha)**2-cos(beta)**2-cos(gamm)**2  &
                                  +2.0*cos(alpha)*cos(beta)*cos(gamm) )
  real,parameter::Vol=a*b*c*nu
  real,parameter::R_pseudo=0.0
  real,parameter::cutoff_resolution=2.0
  integer,parameter::hnum=38
  integer,parameter::knum=96
  integer,parameter::lnum=104
  integer i,j,ios,obsstatus
  integer h,k,l,h_min,h_max,k_min,k_max,l_min,l_max
  integer irefl,num_of_uniq_refl
  real S,resolution,obssf,radius
  real real_part,imag_part
  real cache(6)
  real random
  real::measured_refl(3,0:hnum/2,0:knum/2,0:lnum/2)
  real,allocatable::list_of_uniq_refl(:,:)
  character(len=80) filename

  !initialization
  num_of_uniq_refl=0
  do l=0,lnum/2,1
    do k=0,knum/2,1
      do h=0,hnum/2,1
        measured_refl(1,h,k,l)=0 !obsstatus
        measured_refl(2,h,k,l)=0 !obssf
        S=sqrt( h**2.0*S1+k**2.0*S2+l**2.0*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6 )/Vol
        resolution=1.0/S
        measured_refl(3,h,k,l)=resolution
        if(resolution>=cutoff_resolution)then
!          if(h==0 .and. k==0 .and. mod(l,2)/=0)then
!            cycle !conditions limiting possible reflections
!          end if
!          if(h==0 .and. mod(k,2)/=0 .and. l==0)then
!            cycle !conditions limiting possible reflections
!          end if
!          if(mod(h,2)/=0 .and. k==0 .and. l==0)then
!            cycle !conditions limiting possible reflections
!          end if
          num_of_uniq_refl=num_of_uniq_refl+1
        end if
      end do
    end do
  end do

  write(*,*)'num_of_uniq_refl=',num_of_uniq_refl

  allocate(list_of_uniq_refl(num_of_uniq_refl,6))

  !calculating all unique reflections incluing both experimental measured and missed reflections
  write(filename,'(A24,A4,A4)')'observed_reflections_of_',PDB_ID,'.txt'
  open(unit=10,file=filename,status='old')
  h_max=-huge(h_max)
  k_max=-huge(k_max)
  l_max=-huge(l_max)
  do
    read(10,'(3I5,I2,F12.3)',iostat=ios)h,k,l,obsstatus,obssf
    if(ios/=0)exit
    if(h_max<abs(h))h_max=abs(h)
    if(k_max<abs(k))k_max=abs(k)
    if(l_max<abs(l))l_max=abs(l)
  end do
  h_min=0 
  k_min=0
  l_min=0
  !double check the gridsize is good 
  write(*,*)h_max,'should <=', hnum/2
  write(*,*)k_max,'should <=', knum/2
  write(*,*)l_max,'should <=', lnum/2
  rewind(10)

  do
    read(10,'(3I5,I2,F12.3)',iostat=ios)h,k,l,obsstatus,obssf
    if(ios/=0)exit
    S=sqrt( h**2.0*S1+k**2.0*S2+l**2.0*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6 )/Vol
    resolution=1.0/S
    measured_refl(1,h,k,l)=obsstatus
    measured_refl(2,h,k,l)=obssf !*exp(-(pi*R_pseudo*S)**2)
    measured_refl(3,h,k,l)=resolution
  end do
  close(10)

  i=0
  irefl=0
  do l=0,lnum/2,1
    do k=0,knum/2,1
      do h=0,hnum/2,1
        S=sqrt( h**2.0*S1+k**2.0*S2+l**2.0*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6 )/Vol
        resolution=1.0/S
        if(resolution>=cutoff_resolution)then 
          obsstatus=measured_refl(1,h,k,l)
          obssf=measured_refl(2,h,k,l)
          if( (h==0 .and. k==0 .and. mod(l,2)/=0) .or. &
              (h==0 .and. mod(k,2)/=0 .and. l==0) .or. &
              (mod(h,2)/=0 .and. k==0 .and. l==0)         )then
            obsstatus=1
            obssf=0
          end if
          irefl=irefl+1
          list_of_uniq_refl(irefl,1)=h
          list_of_uniq_refl(irefl,2)=k
          list_of_uniq_refl(irefl,3)=l
          list_of_uniq_refl(irefl,4)=obsstatus
          list_of_uniq_refl(irefl,5)=obssf
          list_of_uniq_refl(irefl,6)=resolution
        end if
      end do
    end do
  end do 
  num_of_uniq_refl=irefl



  !step 2: arranging unique reflections according to their corresponding 'radius'
  do i=1,num_of_uniq_refl-1,1
    do j=i+1,num_of_uniq_refl,1
      if(list_of_uniq_refl(j,6)>list_of_uniq_refl(i,6))then
        cache(1:6)=list_of_uniq_refl(i,1:6)
        list_of_uniq_refl(i,1:6)=list_of_uniq_refl(j,1:6)
        list_of_uniq_refl(j,1:6)=cache(1:6)
      end if
    end do
  end do

  write(filename,'(A22,A4,A4)')'unique_reflections_of_',PDB_ID,'.txt'
  open(unit=11,file=filename,status='replace')
!  i=0
  do irefl=1,num_of_uniq_refl,1
    h=int(list_of_uniq_refl(irefl,1)) 
    k=int(list_of_uniq_refl(irefl,2)) 
    l=int(list_of_uniq_refl(irefl,3)) 
    obsstatus=int(list_of_uniq_refl(irefl,4))
    obssf=list_of_uniq_refl(irefl,5) 
    resolution=list_of_uniq_refl(irefl,6)
    if(h==0 .and. k==0 .and. l==0)resolution=999.0
!    if(resolution>15.)obsstatus=0
    write(11,'(3I5,I2,F12.3,F8.3)')h,k,l,obsstatus,obssf,resolution
  end do
!  write(*,*)real(i)/real(num_of_uniq_refl)*100.0  !98.447227 
  close(11)

  deallocate(list_of_uniq_refl)
 
  stop
end

!-----------------------------------------------------------------------------80 i=0
