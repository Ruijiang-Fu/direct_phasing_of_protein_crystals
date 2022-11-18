program main
  implicit none
  integer i,j,ios,ifile0,ifile
  integer flag
  integer iter0,iter,iorig0
  integer numIter0
  integer,parameter::numFile=200
  integer,parameter::numIter=30000
  integer,parameter::numResoShell=9
  character(len=80) fileName0, fileName
  real meanPhaseErrorFile(numFile,0:numIter-1),meanPhaseError(0:numResoShell-1),meanPhaseError0
  real corrCoefFile(numFile,0:numIter-1),corrCoef(0:numResoShell-1)
  real protMaskMatchFile(numFile,0:numIter-1),protMaskMatch
  real RfreeFile(numFile,0:numIter-1),Rfree
  real RworkFile(numFile,0:numIter-1),Rwork
  real F000File(numFile,0:numIter-1),F000
  real scaleFactFile(numFile,0:numIter-1),scaleFact
  real centMassDeviDistFile(numFile,0:numIter-1),centMassDeviDist
  real ncsAxisDeviAnglFile(numFile,0:numIter-1),ncsAxisDeviAngl
  character(len=20) :: myfmt
  character(len=80) :: string

  open(unit=11,file='0_log_of_ncs_cent_mass_all_runs.pdb',status='replace')

  ifile0=0
  do ifile=1,numFile,1

    if(ifile<10)then
      write(fileName0,'(A27,I1)')'directPhasing/all_runs/run_',ifile
    else if(ifile<100)then
      write(fileName0,'(A27,I2)')'directPhasing/all_runs/run_',ifile
    else if(ifile<1000)then
      write(fileName0,'(A27,I3)')'directPhasing/all_runs/run_',ifile
    end if

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/mean_phase_error_reso_sphere.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

	ifile0=ifile0+1

    do 
      read(10,'(I8,9(1X,F8.1))',iostat=ios) iter,(meanPhaseError(i),i=0,8)
      if(ios/=0)exit
      meanPhaseErrorFile(ifile0,iter)=meanPhaseError(0)
    end do
    close(10)

	numIter0 = iter

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/corr_coef.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,9(1X,F8.3))',iostat=ios) iter,(corrCoef(i),i=0,8)
      if(ios/=0)exit
      corrCoefFile(ifile0,iter)=corrCoef(0)
    end do
    close(10)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/prot_mask_match.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,protMaskMatch
      if(ios/=0)exit
      protMaskMatchFile(ifile0,iter)=protMaskMatch
    end do
    close(10)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/R_free.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle
    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,Rfree
      if(ios/=0)exit
      RfreeFile(ifile0,iter)=Rfree
    end do
    close(10)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/R_work.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle
    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,Rwork
      if(ios/=0)exit
      RworkFile(ifile0,iter)=Rwork
    end do
    close(10)

	




  end do

  close(11)

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  open(unit=11,file='0_log_of_deltaPhi_all_runs.txt',status='replace')
  open(unit=12,file='0_log_of_corrCoef_all_runs.txt',status='replace')
  open(unit=13,file='0_log_of_protMaskMatch_all_runs.txt',status='replace')
  open(unit=14,file='0_log_of_Rfree_all_runs.txt',status='replace')
  open(unit=15,file='0_log_of_Rwork_all_runs.txt',status='replace')
  do iter=0,numIter0-1,1
	if(mod(iter,20)/=0)cycle

    write(myfmt, '("(I8,",I0,"(A1,F6.1))")') ifile0
    write(11, fmt=myfmt ) iter, (',',meanPhaseErrorFile(i,iter), i=1,ifile0)

    write(myfmt, '("(I8,",I0,"(A1,F6.3))")') ifile0
    write(12, fmt=myfmt ) iter, (',',corrCoefFile(i,iter), i=1,ifile0)

    write(13, fmt=myfmt ) iter, (',',protMaskMatchFile(i,iter), i=1,ifile0)

    write(14, fmt=myfmt ) iter, (',',RfreeFile(i,iter), i=1,ifile0)

    write(15, fmt=myfmt ) iter, (',',RworkFile(i,iter), i=1,ifile0)

  end do
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  stop
end 
