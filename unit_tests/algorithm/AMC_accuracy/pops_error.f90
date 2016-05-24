program compare_error
  implicit none
  ! model : cosine /shu1d
  character*32,parameter::model='cosine'
  ! mode : FullyRandom /QuasiRandom
  character*32,parameter::mode='FullyRandom'
  integer,parameter::nlev=21,ncell=64,Nsnr=9
  integer,dimension(Nsnr)::snr=(/ 10,15,20,30,40,50,60,80,100 /)
  
  real*8,dimension(ncell,nlev)::pops_FMC,pops_QMC,pops_ref
  character*32 filenameFR,filenameQR,filenameOut,case_name
  real*8 Rf,Rq ! Radius of Quasi-random data and Fully-random data 
  real*8 L2,L_inf ! L_2-norm error and L_inf-norm error
  real*8 diff ! difference of pops_FMC and pops_QMC
  real*8 FL2,FL_inf,QL2,QL_inf
  integer nL2,nFL2,nQL2
  integer i,j,k ! iterator
  
  write (case_name, "(A,A,A5)") adjustl(trim(model)),adjustl(trim(mode)),"_snr_"
  
  ! read reference solution : FMC snr 100
  write (filenameFR, "(A,I3,A4)") case_name, 800, '.dat'
  open(unit=8,FILE=filenameFR,STATUS='OLD',FORM='FORMATTED')
  do j=1,ncell
    ! Read from the population file
    read(8,*) Rf,(pops_ref(j,k),k=1,nlev)
  enddo
  close(8)
  
  write (filenameOut, "(A,A)") adjustl(trim(model)),'_error.dat'
  open(unit=10,FILE=filenameOut,STATUS='REPLACE',FORM='FORMATTED')
  do i=1,Nsnr
    
    ! set the target file name of population data
    if (snr(i) < 100 )then
      write (filenameFR, "(A,I3,A4)") case_name, snr(i), '.dat'
      write (filenameQR, "(A,I3,A4)") case_name, snr(i), '.dat'
    else
      write (filenameFR, "(A,I3,A4)") case_name, snr(i), '.dat'
      write (filenameQR, "(A,I3,A4)") case_name, snr(i), '.dat'
    endif
    ! Read from the population file
    open(unit=8,FILE=filenameFR,STATUS='OLD',FORM='FORMATTED')
    open(unit=9,FILE=filenameQR,STATUS='OLD',FORM='FORMATTED')
    do j=1,ncell
      read(8,*) Rf,(pops_FMC(j,k),k=1,nlev)
      read(9,*) Rq,(pops_QMC(j,k),k=1,nlev)
      ! check the grid consistency
      if (Rf /= Rq) then
        ! the grid is not consistent
        write(*,*) 'Attention!!! radius mismatched'
      endif
    enddo
    close(8)
    close(9)
    
    ! reset error
    FL2=0.0
    FL_inf=0.0
    nFL2=0
    
    QL2=0.0
    QL_inf=0.0
    nQL2=0
    
    do j=1,ncell      
      
      do k=1,nlev
        if(pops_ref(j,k)>=1e-6)then
          ! error of FMC to the reference
          diff=abs(pops_FMC(j,k)-pops_ref(j,k))/pops_ref(j,k)
          FL2=FL2+diff*diff
          nFL2=nFL2+1
          if(diff>FL_inf) then
            FL_inf=diff
          endif
          ! error of QMC to the reference
          diff=abs(pops_QMC(j,k)-pops_ref(j,k))/pops_ref(j,k)
          QL2=QL2+diff*diff
          nQL2=nQL2+1
          if(diff>QL_inf) then
            QL_inf=diff
          endif
        endif
      enddo
      
      
    enddo

    
    L2=sqrt(L2/nL2)
    FL2=sqrt(FL2/nFL2)
    QL2=sqrt(QL2/nQL2)
    
    write(10,*) snr(i),FL2,FL_inf,QL2,QL_inf
  enddo
  close(10)

end program