program main
  implicit none
  integer i,j,ios,ifile
  integer ishell
  integer iter
  integer,parameter::numShell=9
  integer,parameter::numFile=200
  character(len=80) fileName0,fileName
  real meanPhaseError(0:numShell-1)
  real corrCoef(0:numShell-1)
  real protMaskMatch
  real Rfree
  real Rwork

  open(unit=11,file='2_log_of_final_deltaPhi_each_run.txt',status='replace')
  open(unit=12,file='2_log_of_final_corrCoef_each_run.txt',status='replace')
  open(unit=13,file='2_log_of_final_protMaskMatch_each_run.txt',status='replace')
  open(unit=14,file='2_log_of_final_Rfree_each_run.txt',status='replace')
  open(unit=15,file='2_log_of_final_Rwork_each_run.txt',status='replace')

  do ifile=1,numFile,1

    if(ifile<10)then
      write(fileName0,'(A27,I1)')'directPhasing/all_runs/run_',ifile
    else if(ifile<100)then
      write(fileName0,'(A27,I2)')'directPhasing/all_runs/run_',ifile
    else if(ifile<1000)then
      write(fileName0,'(A27,I3)')'directPhasing/all_runs/run_',ifile
    end if

	fileName = trim(filename0) // '/mean_phase_error_reso_sphere.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,9(1X,F8.1))',iostat=ios) iter, (meanPhaseError(i), i=0, 8)
      if(ios/=0)exit
    end do
    close(10)

    write(11,'(I8,9(A1,F8.1))')ifile, (',', meanPhaseError(i), i=0, 8)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/corr_coef.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,9(1X,F8.3))',iostat=ios) iter, (corrCoef(i), i=0, 8)
      if(ios/=0)exit
    end do
    close(10)

    write(12,'(I8,9(A1,F8.3))')ifile, (',',corrCoef(i), i=0, 8)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/prot_mask_match.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,protMaskMatch
      if(ios/=0)exit
    end do
    close(10)

    write(13,'(I8,A1,F8.3)')ifile,',',protMaskMatch

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/R_free.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,Rfree
      if(ios/=0)exit
    end do
    close(10)

    write(14,'(I8,A1,F8.3)')ifile,',',Rfree

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/R_work.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,Rwork
      if(ios/=0)exit
    end do
    close(10)

    write(15,'(I8,A1,F8.3)')ifile,',',Rwork

  end do
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)


  stop
end 
