program main
  implicit none
  integer i,j,ios,ifile0,ifile
  integer flag
  integer cnt
  integer iter0,iter,iorig0
  integer numIter0
  integer,parameter::numFile=300
  integer,parameter::numIter=30000
  character(len=80) fileName0, fileName
  real RfreeFile(numFile,0:numIter-1),Rfree
  real RworkFile(numFile,0:numIter-1),Rwork
  real DensCutoffmidFile(numFile,0:numIter-1),DensCutoffmid
  character(len=20) :: myfmt

  ifile0=0
  cnt=0
  do ifile=1,numFile,1

    if(ifile<10)then
      write(fileName0,'(A27,I1)')'directPhasing/all_runs/run_',ifile
    else if(ifile<100)then
      write(fileName0,'(A27,I2)')'directPhasing/all_runs/run_',ifile
    else if(ifile<1000)then
      write(fileName0,'(A27,I3)')'directPhasing/all_runs/run_',ifile
    end if

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '/R_free.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    ifile0=ifile0+1
    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,Rfree
      if(ios/=0)exit
      RfreeFile(ifile0,iter)=Rfree
    end do
    if(Rfree<0.25) cnt=cnt+1
    close(10)
    numIter0 = iter

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

    !---------------------------------------------------------------------------

	fileName = trim(filename0) // '/DensCutoffmid.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle
    do 
      read(10,'(I8,1X,F15.8)',iostat=ios) iter,DensCutoffmid
      if(ios/=0)exit
      DensCutoffmidFile(ifile0,iter)=DensCutoffmid
    end do
    close(10)

	!---------------------------------------------------------------------------

  end do
  write(*,*)cnt
  !-----------------------------------------------------------------------------

  open(unit=14,file='0_log_of_Rfree_all_runs.txt',status='replace')
  open(unit=15,file='0_log_of_Rwork_all_runs.txt',status='replace')
  open(unit=16,file='0_log_of_DensCutoffmid_all_runs.txt',status='replace')
  do iter=0,numIter0-1,1
	if(mod(iter,20)/=0)cycle

    write(myfmt, '("(I8,",I0,"(A1,F6.3))")') ifile0
    write(14, fmt=myfmt ) iter, (',',RfreeFile(i,iter), i=1,ifile0)
    write(15, fmt=myfmt ) iter, (',',RworkFile(i,iter), i=1,ifile0)
    write(myfmt, '("(I8,",I0,"(A1,F15.8))")') ifile0
    write(16, fmt=myfmt ) iter, (',',DensCutoffmidFile(i,iter), i=1,ifile0)


  end do
  close(14)
  close(15)
  close(16)

  stop
end 
