!on cluster: ifort -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
!on PC: ifort -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lm

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
  integer,parameter::anum=40
  integer,parameter::bnum=96
  integer,parameter::cnum=104
  integer,parameter::hnum=anum
  integer,parameter::knum=bnum
  integer,parameter::lnum=cnum
  real,parameter::forward_scale_factor=1.0/Vol
  real,parameter::backward_scale_factor=Vol/real(anum*bnum*cnum)
  real,parameter::cutoff_resolution=2.0
  real,parameter::solvent_content0=0.62
  real,parameter::sigma_weight_initial=4.0
  real,parameter::sigma_weight_final=3.0
  real,parameter::scale_factor=1.0
  integer,parameter::num_of_iterations=10000
  integer,parameter::num_of_iterations_to_limit_den=1000
  integer,parameter::num_of_iterations_for_sol_flat=500
  integer,parameter::num_of_units=300
  integer,parameter::num_of_bins=10000000
  integer,parameter::num_of_origin_choices=16
  real,parameter::limit_density_initial=1.0
  integer i,j,ios
  integer h,k,l,obsstatus
  integer ia,ib,ic,ia0,ib0,ic0
  integer ih,ik,il
  real sign_of_phase,delta_phase
  integer num_of_grid_points_in_protein_region
  integer num_of_grid_points_in_each_histogram_unit
  integer num_of_grid_points_shown_in_ASU
  integer iteration
  integer ibin
  integer iunit
  integer iorigin
  real S,resolution,resolution0
  real sigma_weight,solvent_content
  real x,y,z,x0,y0,z0,cache(3)
  real density,density_min,density_ave,density_max
  real density1,density2,cutoff_density
  real limit_density
  real cutoff_weighted_density
  real obssf,calsf,phase,random
  real real_part,imag_part
  real k_factor,R_value
  real CC(num_of_origin_choices),delta_phi(num_of_origin_choices)
  real sum_of_obssf,sum_of_calsf
  real sum_of_calsf_times_obssf,sum_of_calsf_squ
  real obssfsqu,calsfsqu,obssftimescalsf,obssfsum
  real obssf_times_calsf_times_cosine
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

  real::list_of_uniq_refl(5,0:hnum/2,0:knum/2,0:lnum/2)
  real::density_on_grid_points_in_ASU(0:anum/2-1,0:bnum/2-1,0:cnum-1)
  real::density_on_grid_points_in_ASU_tmp(0:anum/2-1,0:bnum/2-1,0:cnum-1)
  real::density_on_grid_points_in_unitcell(0:anum-1,0:bnum-1,0:cnum-1)
  real::weighted_density_on_grid_points_in_ASU(0:anum/2-1,0:bnum/2-1,0:cnum-1)

  integer::height_of_each_bin(num_of_bins) !frequency, integer
  integer::num_of_grid_points_in_units(num_of_units)
  real::std_boundary_of_each_unit(0:num_of_units)
  real::cal_boundary_of_each_unit(0:num_of_units)
  real::ha(num_of_units),hb(num_of_units)

  real::correct_phases(num_of_origin_choices,0:hnum/2,0:knum/2,0:lnum/2)

  open(unit=9,file='log_of_F000.txt',status='replace')

  open(unit=20,file='log_of_R_work_resolution_above_1A.txt',status='replace')
  open(unit=21,file='log_of_R_work_resolution_above_2A.txt',status='replace')
  open(unit=22,file='log_of_R_work_resolution_above_3A.txt',status='replace')
  open(unit=23,file='log_of_R_work_resolution_above_6A.txt',status='replace')
  open(unit=24,file='log_of_R_work_resolution_above_10A.txt',status='replace')
  open(unit=25,file='log_of_R_work_resolution_above_15A.txt',status='replace')

  open(unit=30,file='log_of_R_free_resolution_above_1A.txt',status='replace')
  open(unit=31,file='log_of_R_free_resolution_above_2A.txt',status='replace')
  open(unit=32,file='log_of_R_free_resolution_above_3A.txt',status='replace')
  open(unit=33,file='log_of_R_free_resolution_above_6A.txt',status='replace')
  open(unit=34,file='log_of_R_free_resolution_above_10A.txt',status='replace')
  open(unit=35,file='log_of_R_free_resolution_above_15A.txt',status='replace')

  open(unit=40,file='log_of_CC_and_deltaphi_resolution_above_1A.txt',status='replace')
  open(unit=41,file='log_of_CC_and_deltaphi_resolution_above_2A.txt',status='replace')
  open(unit=42,file='log_of_CC_and_deltaphi_resolution_above_3A.txt',status='replace')
  open(unit=43,file='log_of_CC_and_deltaphi_resolution_above_6A.txt',status='replace')
  open(unit=44,file='log_of_CC_and_deltaphi_resolution_above_10A.txt',status='replace')
  open(unit=45,file='log_of_CC_and_deltaphi_resolution_above_15A.txt',status='replace')

 
  call initialize_random_seed()
  tab=char(9)

  length_of_each_dimension(1)=anum
  length_of_each_dimension(2)=bnum
  length_of_each_dimension(3)=cnum

  !initialization
  do l=0,lnum/2,1
    do k=0,knum/2,1
      do h=0,hnum/2,1
        list_of_uniq_refl(1,h,k,l)=0 !obsstatus
        list_of_uniq_refl(2,h,k,l)=0 !obssf
        list_of_uniq_refl(3,h,k,l)=0 !calsf
        list_of_uniq_refl(4,h,k,l)=0 !phase
        S=sqrt( h**2*S1+k**2*S2+l**2*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6 )/Vol
        !if use resolution=Vol/sqrt(...), it will give 'NAN' when h,k,l are zero.
        resolution=1.0/S
        list_of_uniq_refl(5,h,k,l)=resolution
      end do
    end do
  end do

  !read experimental observed structure factors
  write(filename,'(A22,2A4)')'unique_reflections_of_',PDB_ID,'.txt'
  open(unit=10,file=filename,status='old')
  do
    read(10,'(3I5,I2,F12.3,F8.3)',iostat=ios)h,k,l,obsstatus,obssf,resolution
    if(ios/=0)exit
    if(resolution>=cutoff_resolution)then
      list_of_uniq_refl(1,h,k,l)=obsstatus
      list_of_uniq_refl(2,h,k,l)=obssf
    end if
  end do
  close(10)

  do l=0,lnum/2,1
    do k=0,knum/2,1
      do h=0,hnum/2,1
        do iorigin=1,num_of_origin_choices,1
          correct_phases(iorigin,h,k,l)=0
        end do
      end do
    end do
  end do

  !read correct phases to prepare for the calculation of CC and deltaphi
  do iorigin=1,num_of_origin_choices,1
    if(iorigin<10)then
      write(filename,'(A30,I1,A4)')'true_phases_for_origin_choice_', &
                         iorigin,'.txt'
    else if(iorigin<100)then
      write(filename,'(A30,I2,A4)')'true_phases_for_origin_choice_', &
                         iorigin,'.txt'
    end if
    open(unit=10,file=filename,status='old')
    do
      read(10,'(3I5,I2,F12.3,F8.3)',iostat=ios)h,k,l,obsstatus,obssf,phase
      if(ios/=0)exit
      S=sqrt( h**2*S1+k**2*S2+l**2*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6 )/Vol
      !if use resolution=Vol/sqrt(...), it will give 'NAN' when h,k,l are zero.
      resolution=1.0/S
      if(resolution>=cutoff_resolution)then
        phase=phase/180.0*pi
        correct_phases(iorigin,h,k,l)=phase
      end if
    end do
    close(10)
  end do

  !read standard histogram as a reference hitogram inside the protein boundary
  write(filename,'(A23,I2,A4)')'std_histogram_with_sol_',  &
                               nint(solvent_content0*100.0),'.txt'
  open(unit=11,file=filename,status='old')
  do
    read(11,*,iostat=ios)iunit,std_boundary_of_each_unit(iunit)
    if(ios/=0)exit
  end do
  close(11)

  !generate random electron density inside the initial protein boundary
  open(unit=13,file='random_number.txt',status='replace')

  do ic=0,cnum-1,1
    do ib=0,bnum/2-1,1
      do ia=0,anum/2-1,1
        call random_number(random)
        write(13,*)random
      end do
    end do
  end do
  rewind(13)

  do ic=0,cnum-1,1
    do ib=0,bnum/2-1,1
      do ia=0,anum/2-1,1
        read(13,*)random
        density_on_grid_points_in_ASU(ia,ib,ic)=random
      end do
    end do
  end do
  close(13)


  !=======================iterations start form here==========================80


  sigma_weight=sigma_weight_initial  !will decrease. It coresponds to average_radius.

  solvent_content=solvent_content0  !solvent content is fixed (can be changed)

  do iteration=1,num_of_iterations,1


    do ic=0,cnum-1,1
      do ib=0,bnum/2-1,1
        do ia=0,anum/2-1,1
          cache(1)=real(ia+0.5)/real(anum)
          cache(2)=real(ib+0.5)/real(bnum)
          cache(3)=real(ic+0.5)/real(cnum)
          density=density_on_grid_points_in_ASU(ia,ib,ic)
          x=cache(1)
          y=cache(2)
          z=cache(3)
          ia0=int(x*anum)
          ib0=int(y*bnum)
          ic0=int(z*cnum)
          density_on_grid_points_in_unitcell(ia0,ib0,ic0)=density
          x=-cache(1)+0.5
          y=-cache(2)+1.0
          z= cache(3)+0.5
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinate_into_unitcell(x,y,z)
          end if
          ia0=int(x*anum)
          ib0=int(y*bnum)
          ic0=int(z*cnum)
          density_on_grid_points_in_unitcell(ia0,ib0,ic0)=density
          x=-cache(1)+1.0
          y= cache(2)+0.5
          z=-cache(3)+0.5
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinate_into_unitcell(x,y,z)
          end if
          ia0=int(x*anum)
          ib0=int(y*bnum)
          ic0=int(z*cnum)
          density_on_grid_points_in_unitcell(ia0,ib0,ic0)=density
          x= cache(1)+0.5
          y=-cache(2)+0.5
          z=-cache(3)+1.0
          if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)then
            continue
          else
            call adjust_fractional_coordinate_into_unitcell(x,y,z)
          end if
          ia0=int(x*anum)
          ib0=int(y*bnum)
          ic0=int(z*cnum)
          density_on_grid_points_in_unitcell(ia0,ib0,ic0)=density
        end do
      end do
    end do

    !read real density into complex 3D matrix
    do ic=0,cnum-1,1
      do ib=0,bnum-1,1
        do ia=0,anum-1,1
          real_space_3D(ia,ib,ic)=density_on_grid_points_in_unitcell(ia,ib,ic)
        end do
      end do
    end do

    !Perform a complex to complex backword Fourier transform (definition of FFT)
    Status = DftiCreateDescriptor( my_descriptor_handle, DFTI_SINGLE,  &
                                   DFTI_COMPLEX, 3, length_of_each_dimension )
    Status = DftiSetValue( my_descriptor_handle, DFTI_PLACEMENT,  &
                           DFTI_NOT_INPLACE )
    Status = DftiSetValue( my_descriptor_handle, DFTI_BACKWARD_SCALE,  &
                           backward_scale_factor )
    Status = DftiCommitDescriptor( my_descriptor_handle )
    Status = DftiComputeBackward( my_descriptor_handle, real_space,  &
                                  reciprocal_space )
    Status = DftiFreeDescriptor( my_descriptor_handle )

    !When we put the density into the complex matrix before FFT, the density in 
    !real space are located in the middle of each grid point. However, the FFT 
    !program supposes the density are located on each grid point. After FFT 
    !calculation, we need to time a half_grid phase to get correct answer.

    do il=0,lnum-1,1
      do ik=0,knum-1,1
        do ih=0,hnum-1,1
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
          * ( cos(2.0*pi*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum))) &
          +IM*sin(2.0*pi*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum))) )

          !record the unique calculated phases and the unique calculated magnitudes
          if(h<0 .or. k<0 .or. l<0)cycle  !not unique reflections
          resolution=list_of_uniq_refl(5,h,k,l)
          if(resolution<cutoff_resolution)cycle
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
          if(abs(calsf)<0.1)then  !not good to find phase
            calsf=0
            phase=0
          end if
          phase=mod(phase,2.0*pi)
          if(phase<0) phase=phase+2.0*pi
          list_of_uniq_refl(3,h,k,l)=calsf
          list_of_uniq_refl(4,h,k,l)=phase
        end do
      end do
    end do

    !calculate R_work above different resolutions
   if(iteration==1 .or. mod(iteration,10)==0)then

    resolution0=15
    do while (resolution0>=int(cutoff_resolution))
      if(int(resolution0)/=2 .and. int(resolution0)/=3 .and. &
         int(resolution0)/=6 .and. int(resolution0)/=10 .and. &
         int(resolution0)/=15 ) then
        resolution0=resolution0-1.0
        cycle
      end if
      sum_of_obssf=0
      sum_of_calsf=0
      sum_of_calsf_times_obssf=0
      sum_of_calsf_squ=0
      do l=0,lnum/2-1,1
        do k=0,knum/2-1,1  
          do h=0,hnum/2-1,1
            if(h==0 .and. k==0 .and. l==0)cycle  !special reflection F000
            resolution=list_of_uniq_refl(5,h,k,l)
            if(resolution<cutoff_resolution)cycle
            if(resolution<resolution0) cycle
            obsstatus=nint(list_of_uniq_refl(1,h,k,l))
            obssf=list_of_uniq_refl(2,h,k,l)
            calsf=list_of_uniq_refl(3,h,k,l)
            if(obsstatus==1)then
              sum_of_obssf=sum_of_obssf+obssf
              sum_of_calsf=sum_of_calsf+calsf
              sum_of_calsf_times_obssf=sum_of_calsf_times_obssf+calsf*obssf
              sum_of_calsf_squ=sum_of_calsf_squ+calsf**2
            end if
          end do
        end do
      end do
      k_factor=sum_of_calsf_times_obssf/sum_of_calsf_squ
  
      R_value=0
      do l=0,lnum/2-1,1
        do k=0,knum/2-1,1
          do h=0,hnum/2-1,1
            if(h==0 .and. k==0 .and. l==0)cycle  !special reflection F000
            resolution=list_of_uniq_refl(5,h,k,l)
            if(resolution<cutoff_resolution)cycle
            if(resolution<resolution0) cycle
            obsstatus=nint(list_of_uniq_refl(1,h,k,l))
            obssf=list_of_uniq_refl(2,h,k,l)
            calsf=list_of_uniq_refl(3,h,k,l)
            if(obsstatus==1)then
              R_value=R_value+abs( obssf-k_factor*calsf )
            end if
          end do
        end do
      end do
      R_value=R_value/sum_of_obssf
      if(int(resolution0)==1)write(20,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==2)write(21,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==3)write(22,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==6)write(23,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==10)write(24,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==15)write(25,*)iteration,tab,k_factor,tab,R_value
      resolution0=resolution0-1.0
    end do

   end if

    !calculate R_free above different resolutions
   if(iteration==1 .or. mod(iteration,10)==0)then

    resolution0=15
    do while (resolution0>=int(cutoff_resolution))
      if(int(resolution0)/=2 .and. int(resolution0)/=3 .and. &
         int(resolution0)/=6 .and. int(resolution0)/=10 .and. &
         int(resolution0)/=15 ) then
        resolution0=resolution0-1.0
        cycle
      end if
      sum_of_obssf=0
      sum_of_calsf=0
      sum_of_calsf_times_obssf=0
      sum_of_calsf_squ=0
      do l=0,lnum/2-1,1
        do k=0,knum/2-1,1  
          do h=0,hnum/2-1,1
            if(h==0 .and. k==0 .and. l==0)cycle  !special reflection F000
            resolution=list_of_uniq_refl(5,h,k,l)
            if(resolution<cutoff_resolution)cycle
            if(resolution<resolution0) cycle
            obsstatus=nint(list_of_uniq_refl(1,h,k,l))
            obssf=list_of_uniq_refl(2,h,k,l)
            calsf=list_of_uniq_refl(3,h,k,l)
            if(obsstatus==2)then
              sum_of_obssf=sum_of_obssf+obssf
              sum_of_calsf=sum_of_calsf+calsf
              sum_of_calsf_times_obssf=sum_of_calsf_times_obssf+calsf*obssf
              sum_of_calsf_squ=sum_of_calsf_squ+calsf**2
            end if
          end do
        end do
      end do
      k_factor=sum_of_calsf_times_obssf/sum_of_calsf_squ
  
      R_value=0
      do l=0,lnum/2-1,1
        do k=0,knum/2-1,1
          do h=0,hnum/2-1,1
            if(h==0 .and. k==0 .and. l==0)cycle  !special reflection F000
            resolution=list_of_uniq_refl(5,h,k,l)
            if(resolution<cutoff_resolution)cycle
            if(resolution<resolution0) cycle
            obsstatus=nint(list_of_uniq_refl(1,h,k,l))
            obssf=list_of_uniq_refl(2,h,k,l)
            calsf=list_of_uniq_refl(3,h,k,l)
            if(obsstatus==2)then
              R_value=R_value+abs( obssf-k_factor*calsf )
            end if
          end do
        end do
      end do
      R_value=R_value/sum_of_obssf
      if(int(resolution0)==1)write(30,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==2)write(31,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==3)write(32,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==6)write(33,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==10)write(34,*)iteration,tab,k_factor,tab,R_value
      if(int(resolution0)==15)write(35,*)iteration,tab,k_factor,tab,R_value
      resolution0=resolution0-1.0
    end do

   end if


    !for missing reflections, use scaled calculated ones as a substitude
    sum_of_obssf=0
    sum_of_calsf=0
    sum_of_calsf_times_obssf=0
    sum_of_calsf_squ=0
    do l=0,lnum/2-1,1
      do k=0,knum/2-1,1  
        do h=0,hnum/2-1,1
          if(h==0 .and. k==0 .and. l==0)cycle  !special reflection F000
          resolution=list_of_uniq_refl(5,h,k,l)
          if(resolution<cutoff_resolution) cycle 
          obsstatus=nint(list_of_uniq_refl(1,h,k,l))
          obssf=list_of_uniq_refl(2,h,k,l)
          calsf=list_of_uniq_refl(3,h,k,l)
          if(obsstatus==1)then
            sum_of_obssf=sum_of_obssf+obssf
            sum_of_calsf=sum_of_calsf+calsf
            sum_of_calsf_times_obssf=sum_of_calsf_times_obssf+calsf*obssf
            sum_of_calsf_squ=sum_of_calsf_squ+calsf**2
          end if
        end do
      end do
    end do
    k_factor=sum_of_calsf_times_obssf/sum_of_calsf_squ

    do l=0,lnum/2-1,1
      do k=0,knum/2-1,1
        do h=0,hnum/2-1,1
          if(h==0 .and. k==0 .and. l==0)cycle  !special reflection F000
          resolution=list_of_uniq_refl(5,h,k,l)
          if(resolution<cutoff_resolution) cycle
          obsstatus=nint(list_of_uniq_refl(1,h,k,l))
          obssf=list_of_uniq_refl(2,h,k,l)
          calsf=list_of_uniq_refl(3,h,k,l)
          if(obsstatus/=1)then
            list_of_uniq_refl(3,h,k,l)=sum_of_obssf/sum_of_calsf*calsf  !!!!!!
          end if
        end do
      end do
    end do

    write(9,*)iteration,tab,list_of_uniq_refl(3,0,0,0)

    !monitor deltaphi
   if(iteration==1 .or. mod(iteration,10)==0)then

    resolution0=15
    do while (resolution0>=int(cutoff_resolution))
      if(int(resolution0)/=2 .and. int(resolution0)/=3 .and. &
         int(resolution0)/=6 .and. int(resolution0)/=10 .and. &
         int(resolution0)/=15 ) then
        resolution0=resolution0-1.0
        cycle
      end if
      do iorigin=1,num_of_origin_choices,1
        delta_phi(iorigin)=0
      end do
      i=0
      do l=0,lnum/2-1,1
        do k=0,knum/2-1,1
          do h=0,hnum/2-1,1
            if(h==0 .and. k==0 .and. l==0)cycle  !special reflection F000
            resolution=list_of_uniq_refl(5,h,k,l)
            if(resolution<resolution0) cycle
            if(resolution<cutoff_resolution) cycle
            obsstatus=nint(list_of_uniq_refl(1,h,k,l))
            phase=list_of_uniq_refl(4,h,k,l)
            if(obsstatus/=0)then
              i=i+1
              do iorigin=1,num_of_origin_choices,1
                delta_phi(iorigin)=delta_phi(iorigin) &
                    +acos(cos(correct_phases(iorigin,h,k,l)-phase))
              end do
            end if
          end do
        end do
      end do
      do iorigin=1,num_of_origin_choices,1
        delta_phi(iorigin)=delta_phi(iorigin)/real(i)/pi*180.0
      end do
      if(int(resolution0)==1)write(40,'(I8,16(A1,F8.4))')iteration, &
       ',',delta_phi(1),',',delta_phi(2),',',delta_phi(3),',',delta_phi(4), &
       ',',delta_phi(5),',',delta_phi(6),',',delta_phi(7),',',delta_phi(8), &
       ',',delta_phi(9),',',delta_phi(10),',',delta_phi(11),',',delta_phi(12), &
       ',',delta_phi(13),',',delta_phi(14),',',delta_phi(15),',',delta_phi(16)
      if(int(resolution0)==2)write(41,'(I8,16(A1,F8.4))')iteration, &
       ',',delta_phi(1),',',delta_phi(2),',',delta_phi(3),',',delta_phi(4), &
       ',',delta_phi(5),',',delta_phi(6),',',delta_phi(7),',',delta_phi(8), &
       ',',delta_phi(9),',',delta_phi(10),',',delta_phi(11),',',delta_phi(12), &
       ',',delta_phi(13),',',delta_phi(14),',',delta_phi(15),',',delta_phi(16)
      if(int(resolution0)==3)write(42,'(I8,16(A1,F8.4))')iteration, &
       ',',delta_phi(1),',',delta_phi(2),',',delta_phi(3),',',delta_phi(4), &
       ',',delta_phi(5),',',delta_phi(6),',',delta_phi(7),',',delta_phi(8), &
       ',',delta_phi(9),',',delta_phi(10),',',delta_phi(11),',',delta_phi(12), &
       ',',delta_phi(13),',',delta_phi(14),',',delta_phi(15),',',delta_phi(16)
      if(int(resolution0)==6)write(43,'(I8,16(A1,F8.4))')iteration, &
       ',',delta_phi(1),',',delta_phi(2),',',delta_phi(3),',',delta_phi(4), &
       ',',delta_phi(5),',',delta_phi(6),',',delta_phi(7),',',delta_phi(8), &
       ',',delta_phi(9),',',delta_phi(10),',',delta_phi(11),',',delta_phi(12), &
       ',',delta_phi(13),',',delta_phi(14),',',delta_phi(15),',',delta_phi(16)
      if(int(resolution0)==10)write(44,'(I8,16(A1,F8.4))')iteration, &
       ',',delta_phi(1),',',delta_phi(2),',',delta_phi(3),',',delta_phi(4), &
       ',',delta_phi(5),',',delta_phi(6),',',delta_phi(7),',',delta_phi(8), &
       ',',delta_phi(9),',',delta_phi(10),',',delta_phi(11),',',delta_phi(12), &
       ',',delta_phi(13),',',delta_phi(14),',',delta_phi(15),',',delta_phi(16)
      if(int(resolution0)==15)write(45,'(I8,16(A1,F8.4))')iteration, &
       ',',delta_phi(1),',',delta_phi(2),',',delta_phi(3),',',delta_phi(4), &
       ',',delta_phi(5),',',delta_phi(6),',',delta_phi(7),',',delta_phi(8), &
       ',',delta_phi(9),',',delta_phi(10),',',delta_phi(11),',',delta_phi(12), &
       ',',delta_phi(13),',',delta_phi(14),',',delta_phi(15),',',delta_phi(16)
      resolution0=resolution0-1.0
    end do

   end if

    !We want the density in real space be located in the middle of each grid 
    !point. However, the FFT program supposes the density are located on each 
    !grid point. Before FFT calculation, we need to time a half_grid phase 
    !to satisfy the default assumption of FFT program.

    do il=0,lnum-1,1
      do ik=0,knum-1,1
        do ih=0,hnum-1,1
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
          obsstatus=nint(list_of_uniq_refl(1,abs(h),abs(k),abs(l)))
          if(obsstatus==1)then
            obssf=list_of_uniq_refl(2,abs(h),abs(k),abs(l))
          else
            obssf=list_of_uniq_refl(3,abs(h),abs(k),abs(l))
          end if
          phase=list_of_uniq_refl(4,abs(h),abs(k),abs(l))
          !expand reflections with positive hkl to the whole reciprocal space
          call expand_reflections(h,k,l,sign_of_phase,delta_phase)
          phase=sign_of_phase*phase+delta_phase*pi
          reciprocal_space_3D(ih,ik,il)=obssf*(cos(phase)+IM*sin(phase))  &
          * ( cos(2.0*pi*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum))) &
          -IM*sin(2.0*pi*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum))) )
        end do
      end do
    end do

    !Perform a complex to complex transform
    Status = DftiCreateDescriptor( my_descriptor_handle, DFTI_SINGLE,  &
                                   DFTI_COMPLEX, 3, length_of_each_dimension )
    Status = DftiSetValue( my_descriptor_handle, DFTI_PLACEMENT,  &
                           DFTI_NOT_INPLACE )
    Status = DftiSetValue( my_descriptor_handle, DFTI_FORWARD_SCALE,  &
                           forward_scale_factor )
    Status = DftiCommitDescriptor( my_descriptor_handle )
    Status = DftiComputeForward( my_descriptor_handle, reciprocal_space,  &
                                 real_space )
    Status = DftiFreeDescriptor( my_descriptor_handle )

    !read density from the complex 3D matrix; the density are already located in
    !the middle of each grid point

    do ic=0,cnum-1,1
      do ib=0,bnum/2-1,1
        do ia=0,anum/2-1,1
          density_on_grid_points_in_ASU_tmp(ia,ib,ic)  &
          =real(real_space_3D(ia,ib,ic))
        end do
      end do
    end do

    if(sigma_weight>sigma_weight_final)then
      sigma_weight=sigma_weight_initial-(sigma_weight_initial-sigma_weight_final)  &
       *real(iteration)/real(num_of_iterations-num_of_iterations_for_sol_flat)
    end if

    !use density convolution in reciprocal space to find the weighted average 
    !electron density in real space
    !time a half-grid phase factor and a convolution factor
    do il=0,lnum-1,1
      do ik=0,knum-1,1
        do ih=0,hnum-1,1
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
          obsstatus=nint(list_of_uniq_refl(1,abs(h),abs(k),abs(l)))
          if(obsstatus==1)then
            obssf=list_of_uniq_refl(2,abs(h),abs(k),abs(l))
          else
            obssf=list_of_uniq_refl(3,abs(h),abs(k),abs(l))
          end if
          phase=list_of_uniq_refl(4,abs(h),abs(k),abs(l))
          !expand reflections with positive hkl to the whole reciprocal space
          call expand_reflections(h,k,l,sign_of_phase,delta_phase)
          phase=sign_of_phase*phase+delta_phase*pi
          S=sqrt(h**2*S1+k**2*S2+l**2*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6)/Vol
          reciprocal_space_3D(ih,ik,il)=obssf*(cos(phase)+IM*sin(phase))  &
          *sqrt(2.0*pi*sigma_weight**2)*exp(-2.0*(pi*sigma_weight*S)**2) &
          * ( cos(2.0*pi*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum))) &
          -IM*sin(2.0*pi*(h*0.5/real(anum)+k*0.5/real(bnum)+l*0.5/real(cnum))) )
        end do
      end do
    end do

    !Perform a complex to complex transform
    Status = DftiCreateDescriptor( my_descriptor_handle, DFTI_SINGLE,  &
                                   DFTI_COMPLEX, 3, length_of_each_dimension )
    Status = DftiSetValue( my_descriptor_handle, DFTI_PLACEMENT,  &
                           DFTI_NOT_INPLACE )
    Status = DftiSetValue( my_descriptor_handle, DFTI_FORWARD_SCALE,  &
                           forward_scale_factor )
    Status = DftiCommitDescriptor( my_descriptor_handle )
    Status = DftiComputeForward( my_descriptor_handle, reciprocal_space,  &
                                 real_space )
    Status = DftiFreeDescriptor( my_descriptor_handle )

    !gradually update the weighted average density
    do ic=0,cnum-1,1
      do ib=0,bnum/2-1,1
        do ia=0,anum/2-1,1  !A1+A2=0.5+0.5=1.0+0=0.9+0.1=1, all are OK
          weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
          =0.5*real(real_space_3D(ia,ib,ic))  &
          +0.5*weighted_density_on_grid_points_in_ASU(ia,ib,ic)  !previous one
        end do
      end do
    end do

    num_of_grid_points_in_protein_region  &
    =int(anum/2*bnum/2*cnum*(1.0-solvent_content))

    density_min=huge(density_min)
    density_max=-huge(density_max)
    do ic=0,cnum-1,1
      do ib=0,bnum/2-1,1
        do ia=0,anum/2-1,1   
          density=weighted_density_on_grid_points_in_ASU(ia,ib,ic)
          if(density<density_min)density_min=density
          if(density>density_max)density_max=density
        end do
      end do
    end do

    !find a proper cutoff density on the weighted average density map to 
    !approximately locate the protein boundary  
    density1=density_min
    density2=density_max
    do
      cutoff_weighted_density=(density2+density1)/2.0
      i=0
      do ic=0,cnum-1,1
        do ib=0,bnum/2-1,1
          do ia=0,anum/2-1,1
            if(weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
               >=cutoff_weighted_density)then
              i=i+1
            end if
          end do
        end do
      end do
      if(i>num_of_grid_points_in_protein_region)then
        density1=cutoff_weighted_density
        if(abs(density2-density1)<10.0**(-5))exit  !exit the loop
      else if(i<num_of_grid_points_in_protein_region)then
        density2=cutoff_weighted_density
        if(abs(density2-density1)<10.0**(-5))exit  !exit the loop
      else if(i==num_of_grid_points_in_protein_region)then
        exit !exit the loop
      end if
    end do
    num_of_grid_points_in_protein_region=i

    !use hybrid input-output (HIO) method to modify the electron density inside
    !the calculated protein boundary
    if(iteration>=( num_of_iterations-num_of_iterations_to_limit_den &
                             -num_of_iterations_for_sol_flat) &
       .and. iteration<(num_of_iterations-num_of_iterations_for_sol_flat))then
      limit_density=limit_density_initial-limit_density_initial  &
           *real( iteration-(num_of_iterations-num_of_iterations_to_limit_den &
                             -num_of_iterations_for_sol_flat) ) &
           /real(num_of_iterations_to_limit_den)
    end if

    if(iteration<(num_of_iterations-num_of_iterations_for_sol_flat))then
      do ic=0,cnum-1,1
        do ib=0,bnum/2-1,1
          do ia=0,anum/2-1,1
            if(weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
               >=cutoff_weighted_density)then  !inside boundary
              density_on_grid_points_in_ASU(ia,ib,ic)  &
              =density_on_grid_points_in_ASU_tmp(ia,ib,ic)
            else  !outside boundary
              density_on_grid_points_in_ASU(ia,ib,ic)  &
              =density_on_grid_points_in_ASU(ia,ib,ic)  &
              -0.9*density_on_grid_points_in_ASU_tmp(ia,ib,ic)
              if(iteration>=( num_of_iterations-num_of_iterations_to_limit_den &
                             -num_of_iterations_for_sol_flat))then
                if(density_on_grid_points_in_ASU(ia,ib,ic)>limit_density)then
                  density_on_grid_points_in_ASU(ia,ib,ic)=limit_density
                else if(density_on_grid_points_in_ASU(ia,ib,ic)<-limit_density)then
                  density_on_grid_points_in_ASU(ia,ib,ic)=-limit_density
                end if   
              end if         
            end if
          end do
        end do
      end do
    end if

    !apply solvent flattening at the last num_of_iterations_for_sol_flat iterations
    if(iteration>=(num_of_iterations-num_of_iterations_for_sol_flat))then
      do ic=0,cnum-1,1
        do ib=0,bnum/2-1,1
          do ia=0,anum/2-1,1
            if(weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
               >=cutoff_weighted_density)then  !inside boundary
              density_on_grid_points_in_ASU(ia,ib,ic)  &
              =density_on_grid_points_in_ASU_tmp(ia,ib,ic)
            else  !outside boundary
              density_on_grid_points_in_ASU(ia,ib,ic)=0
            end if
          end do
        end do
      end do
    end if


    !---histogram matching begins------------
    density_min=huge(density_min)
    density_max=-huge(density_max)
    do ic=0,cnum-1,1
      do ib=0,bnum/2-1,1
        do ia=0,anum/2-1,1   
          if(weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
             >=cutoff_weighted_density)then
            density=density_on_grid_points_in_ASU(ia,ib,ic)
            if(density<density_min)density_min=density
            if(density>density_max)density_max=density
          end if
        end do
      end do
    end do

    do ibin=1,num_of_bins,1
      height_of_each_bin(ibin)=0
    end do

    do ic=0,cnum-1,1
      do ib=0,bnum/2-1,1
        do ia=0,anum/2-1,1
          if(weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
             >=cutoff_weighted_density)then
            density=density_on_grid_points_in_ASU(ia,ib,ic)
            ibin=int( (density-density_min)/(density_max-density_min)  &
                      *num_of_bins )+1
            if(ibin>num_of_bins) ibin=num_of_bins
            height_of_each_bin(ibin)=height_of_each_bin(ibin)+1
          end if
        end do
      end do
    end do

    num_of_grid_points_in_each_histogram_unit  &
    =nint( real(anum/2*bnum/2*cnum*(1.0-solvent_content))/real(num_of_units) )
  
    cal_boundary_of_each_unit(0)=density_min
    iunit=1
    i=0
    do ibin=1,num_of_bins,1
      i=i+height_of_each_bin(ibin)
      if(i>=num_of_grid_points_in_each_histogram_unit)then
        num_of_grid_points_in_units(iunit)=i
        cal_boundary_of_each_unit(iunit)=cal_boundary_of_each_unit(0)  &
           +real(ibin)/real(num_of_bins)*(density_max-density_min)
        iunit=iunit+1
        i=0
      end if
    end do
    cal_boundary_of_each_unit(num_of_units)=density_max

    do iunit=1,num_of_units,1
      ha(iunit)=( std_boundary_of_each_unit(iunit)  &
                  -std_boundary_of_each_unit(iunit-1) )  &
               /( cal_boundary_of_each_unit(iunit)  &
                  -cal_boundary_of_each_unit(iunit-1) )
      hb(iunit)=( std_boundary_of_each_unit(iunit-1)  &
                  *cal_boundary_of_each_unit(iunit)  &
                  -std_boundary_of_each_unit(iunit)  &
                  *cal_boundary_of_each_unit(iunit-1) ) &
               /( cal_boundary_of_each_unit(iunit)  &
                  -cal_boundary_of_each_unit(iunit-1) )
    end do

    do ic=0,cnum-1,1
      do ib=0,bnum/2-1,1  
        do ia=0,anum/2-1,1
          if(weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
             >=cutoff_weighted_density)then
            density=density_on_grid_points_in_ASU(ia,ib,ic)
            do iunit=1,num_of_units,1
              if(density>=cal_boundary_of_each_unit(iunit-1)   &
                .and. density<=cal_boundary_of_each_unit(iunit))then
                density=ha(iunit)*density+hb(iunit)
                density_on_grid_points_in_ASU(ia,ib,ic)=density   
                exit
              end if
            end do
          end if
        end do
      end do  
    end do
    !---histogram matching ends------------


    !record results at the end of eath iteration; this can be deleted---

    if(iteration==num_of_iterations)then
  
      if(iteration<10)then
        write(filename,'(I1,A4)')iteration,'.txt'
      else if(iteration<100)then
        write(filename,'(I2,A4)')iteration,'.txt'
      else if(iteration<1000)then
        write(filename,'(I3,A4)')iteration,'.txt'
      else if(iteration<10000)then
        write(filename,'(I4,A4)')iteration,'.txt'
      else if(iteration<100000)then
        write(filename,'(I5,A4)')iteration,'.txt'
      end if
    
      open(unit=14,file=filename,status='replace')
      do l=0,lnum/2,1
        do k=0,knum/2,1
          do h=0,hnum/2,1
            S=sqrt(h**2*S1+k**2*S2+l**2*S3+2.0*h*k*S4+2.0*h*l*S5+2.0*k*l*S6)/Vol
            resolution=1.0/S
            if(resolution>=cutoff_resolution)then
              obsstatus=nint(list_of_uniq_refl(1,h,k,l))
              obssf=list_of_uniq_refl(2,h,k,l)
              calsf=list_of_uniq_refl(3,h,k,l)
              phase=list_of_uniq_refl(4,h,k,l)
              write(14,'(3I5,I2,F12.3,F8.3)')h,k,l,obsstatus,calsf,phase/pi*180.
            end if
          end do
        end do
      end do  
      close(14)

    end if
  
!    if(iteration==num_of_iterations)then
!
!      if(iteration<10)then
!        write(filename,'(I1,A4)')iteration,'.pdb'
!      else if(iteration<100)then
!        write(filename,'(I2,A4)')iteration,'.pdb'
!      else if(iteration<1000)then
!        write(filename,'(I3,A4)')iteration,'.pdb'
!      else if(iteration<10000)then
!        write(filename,'(I4,A4)')iteration,'.pdb'
!      else if(iteration<100000)then
!        write(filename,'(I5,A4)')iteration,'.pdb'
!      end if
!      open(unit=15,file=filename,status='replace')
!      write(15,'(A6,3F9.3,3F7.2,A11,4X,A2,9X)')  &
!           'CRYST1',a,b,c,alpha/pi*180.0,beta/pi*180.0,gamm/pi*180.0, &
!           ' P 43 21 2 ','40'
!      i=0
!      do ic=0,cnum-1,1
!        do ib=0,bnum/2-1,1
!          do ia=0,anum/2-1,1
!            if(weighted_density_on_grid_points_in_ASU(ia,ib,ic)  &
!               >=cutoff_weighted_density)then
!              cache(1)=real(ia+0.5)/real(anum)
!              cache(2)=real(ib+0.5)/real(bnum)
!              cache(3)=real(ic+0.5)/real(cnum)
!              i=i+1
!              x=cache(1)
!              y=cache(2)
!              z=cache(3)
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                   (a,b,c,alpha,beta,gamm,x,y,z)
!              write(15,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!              x=-cache(1)+0.5
!              y=-cache(2)+1.0
!              z= cache(3)+0.5
!              if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and.  &
!                 z<1.0)then
!                continue
!              else
!                call adjust_fractional_coordinate_into_unitcell(x,y,z)
!              end if
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                  (a,b,c,alpha,beta,gamm,x,y,z)
!              i=i+1
!              write(15,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!              x=-cache(1)+1.0
!              y= cache(2)+0.5
!              z=-cache(3)+0.5
!              if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and.  &
!                 z<1.0)then
!                continue
!              else
!                call adjust_fractional_coordinate_into_unitcell(x,y,z)
!              end if
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                   (a,b,c,alpha,beta,gamm,x,y,z)
!              i=i+1
!              write(15,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!              x= cache(1)+0.5
!              y=-cache(2)+0.5
!              z=-cache(3)+1.0
!              if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and.  &
!                 z<1.0)then
!                continue
!              else
!                call adjust_fractional_coordinate_into_unitcell(x,y,z)
!              end if
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                   (a,b,c,alpha,beta,gamm,x,y,z)
!              i=i+1
!              write(15,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!            end if
!          end do
!        end do
!      end do
!      close(15)
!
!    end if
!
!
!    if(iteration==num_of_iterations)then
!
!      density_min=huge(density_min)
!      density_max=-huge(density_max)
!      do ic=0,cnum-1,1
!        do ib=0,bnum/2-1,1
!          do ia=0,anum/2-1,1
!            density=density_on_grid_points_in_ASU(ia,ib,ic)
!            if(density<density_min) density_min=density
!            if(density>density_max) density_max=density
!          end do
!        end do
!      end do
!
!      num_of_grid_points_shown_in_ASU=int(anum/2*bnum/2*cnum*0.05)  !!!
!
!      density1=density_min
!      density2=density_max
!      do
!        cutoff_density=(density2+density1)/2.0
!        i=0
!        do ic=0,cnum-1,1
!          do ib=0,bnum/2-1,1
!            do ia=0,anum/2-1,1
!              if(density_on_grid_points_in_ASU(ia,ib,ic)>=cutoff_density)then
!                i=i+1
!              end if
!            end do
!          end do
!        end do
!        if(i>num_of_grid_points_shown_in_ASU)then
!          density1=cutoff_density
!          if(abs(density2-density1)<10.0**(-3))exit  !exit loop
!        else if(i<num_of_grid_points_shown_in_ASU)then
!          density2=cutoff_density
!          if(abs(density2-density1)<10.0**(-3))exit  !exit loop
!        else if(i==num_of_grid_points_shown_in_ASU)then
!          exit  !exit loop
!        end if
!      end do
!      num_of_grid_points_shown_in_ASU=i
!
!      if(iteration<10)then
!        write(filename,'(I1,A5)')iteration,'_.pdb'
!      else if(iteration<100)then
!        write(filename,'(I2,A5)')iteration,'_.pdb'
!      else if(iteration<1000)then
!        write(filename,'(I3,A5)')iteration,'_.pdb'
!      else if(iteration<10000)then
!        write(filename,'(I4,A5)')iteration,'_.pdb'
!      else if(iteration<100000)then
!        write(filename,'(I5,A5)')iteration,'_.pdb'
!      end if
!      open(unit=16,file=filename,status='replace')
!      write(16,'(A6,3F9.3,3F7.2,A11,4X,A2,9X)')  &
!           'CRYST1',a,b,c,alpha/pi*180.0,beta/pi*180.0,gamm/pi*180.0, &
!           ' P 43 21 2 ','40'
!      i=0
!      do ic=0,cnum-1,1
!        do ib=0,bnum/2-1,1
!          do ia=0,anum/2-1,1
!            if(density_on_grid_points_in_ASU(ia,ib,ic)>=cutoff_density)then
!              cache(1)=real(ia+0.5)/real(anum)
!              cache(2)=real(ib+0.5)/real(bnum)
!              cache(3)=real(ic+0.5)/real(cnum)
!              i=i+1
!              x=cache(1)
!              y=cache(2)
!              z=cache(3)
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                   (a,b,c,alpha,beta,gamm,x,y,z)
!              write(16,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!              x=-cache(1)+0.5
!              y=-cache(2)+1.0
!              z= cache(3)+0.5
!              if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and.  &
!                 z<1.0)then
!                continue
!              else
!                call adjust_fractional_coordinate_into_unitcell(x,y,z)
!              end if
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                  (a,b,c,alpha,beta,gamm,x,y,z)
!              i=i+1
!              write(16,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!              x=-cache(1)+1.0
!              y= cache(2)+0.5
!              z=-cache(3)+0.5
!              if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and.  &
!                 z<1.0)then
!                continue
!              else
!                call adjust_fractional_coordinate_into_unitcell(x,y,z)
!              end if
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                   (a,b,c,alpha,beta,gamm,x,y,z)
!              i=i+1
!              write(16,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!              x= cache(1)+0.5
!              y=-cache(2)+0.5
!              z=-cache(3)+1.0
!              if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and.  &
!                 z<1.0)then
!                continue
!              else
!                call adjust_fractional_coordinate_into_unitcell(x,y,z)
!              end if
!              call convert_fractional_coordinate_to_orthogonal_coordinate  &
!                   (a,b,c,alpha,beta,gamm,x,y,z)
!              i=i+1
!              write(16,'(A4,I7,A3,16X,3F8.3,26X)')'ATOM',i,'  H',x,y,z
!            end if
!          end do
!        end do
!      end do
!      write(filename,'(A4,A14)')PDB_ID,'_unitcell1.pdb'
!      open(unit=17,file=filename,status='old')
!      do
!        read(17,'(A80)',iostat=ios)string
!        if(ios/=0)exit
!        write(16,'(A80)')string
!      end do
!      close(17)
!      close(16)

!    end if

  end do  !iteration  

  close(9)
  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  close(25)
  close(30)
  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(40)
  close(41)
  close(42)
  close(43)
  close(44)
  close(45)
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

subroutine adjust_fractional_coordinate_into_unitcell(x,y,z)

  implicit none
  integer i
  real x,y,z
  real fract(3)

  if(x>=0 .and. x<1.0 .and. y>=0 .and. y<1.0 .and. z>=0 .and. z<1.0)return

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

subroutine convert_orthogonal_coordinate_to_fractional_coordinate(a,b,c,alpha,beta,gamm,x,y,z)

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
  
  !convert fractional coordinate to orthogonal coordinate
  fract=MATMUL(matrix,cart)

  x=fract(1)
  y=fract(2)
  z=fract(3)

  return
end

!----------------------------------------------------------------------------80

subroutine convert_fractional_coordinate_to_orthogonal_coordinate(a,b,c,alpha,beta,gamm,x,y,z)

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
  
  !convert fractional coordinate to orthogonal coordinate
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
