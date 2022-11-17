program main
  implicit none
  integer i,j,n,ios,h,k,l
  integer obsstatus
  integer num_of_uniq_refl,num_of_free_refl
  real F_meas_au,F_meas_sigma_au,I_meas_au,I_meas_sigma_au,F_cal,phase_cal
  character(len=4)::PDB_ID_uppercase='1UII'
  character(len=4) PDB_ID_lowercase
  character(len=80) filename,linestring
  character(len=1) refl_status
  real random

  !converting file name from uppercase to lowercase
  do i=1,4,1
    n=iachar(PDB_ID_uppercase(i:i))
    if(n>=48 .and. n<=57)then
      PDB_ID_lowercase(i:i)=achar(n)
    else if( n>=65 .and. n<=90)then
      PDB_ID_lowercase(i:i)=achar(n+32)
    end if
  end do

  open(unit=10,file='1uii-sf.cif',status='old')

  write(linestring,'(A24,A4,A4)')'observed_reflections_of_',PDB_ID_uppercase,'.txt'
  open(unit=11,file=linestring,status='replace')

  num_of_uniq_refl=0
  num_of_free_refl=0
  do  
    read(10,'(A80)',iostat=ios)linestring
    if(ios/=0)exit

    if( linestring(1:5)=='1 1 1' .and. linestring(22:22)=='o')then

      read(linestring,'(5X,3I5,3X,F8.1,1X,F6.1)')h,k,l,F_meas_au,F_meas_sigma_au
      if(F_meas_au>1.0*F_meas_sigma_au)then
        num_of_uniq_refl=num_of_uniq_refl+1
        obsstatus=1
        write(11,'(3I5,I2,F12.3)')h,k,l,obsstatus,F_meas_au
      end if

    else if( linestring(1:5)=='1 1 1' .and. linestring(22:22)=='f')then

      read(linestring,'(5X,3I5,3X,F8.1,1X,F6.1)')h,k,l,F_meas_au,F_meas_sigma_au
      if(F_meas_au>1.0*F_meas_sigma_au)then
        num_of_uniq_refl=num_of_uniq_refl+1
        call random_number(random)
        if (random < 0.1) then
            obsstatus=2
            num_of_free_refl=num_of_free_refl+1
        else
            obsstatus=1
        end if
        write(11,'(3I5,I2,F12.3)')h,k,l,obsstatus,F_meas_au
      end if

    end if

  end do
  write(*,*)'num_of_uniq_refl=',num_of_uniq_refl
  close(10)
  close(11)
  write(*,*)num_of_free_refl,real(num_of_free_refl)/real(num_of_uniq_refl)
  stop
end

!-----------------------------------------------------------------------------80
