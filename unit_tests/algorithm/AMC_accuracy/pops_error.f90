program compare_error
  implicit none
  ! model : Cosine /Shu1d
  character*32,parameter::model='Cosine'
  integer,parameter::nlev=21,ncell=64,Nsnr=12
  integer,dimension(Nsnr)::snr=(/ 10,15,20,30,40,50,60,80,100,200,400,800 /)
  
  real*8,dimension(ncell,nlev)::pops_FMC,pops_QMC,pops_ref
  character*32 filenameFR,filenameQR,filenameOut,case_name
  real*8 Rf,Rq ! Radius of Quasi-random data and Fully-random data 
  real*8 L2,L_inf ! L_2-norm error and L_inf-norm error
  real*8 diff ! difference of pops_FMC and pops_QMC
  real*8 QL2,QL_inf
  real*8,dimension(12)::FL2,FL_inf
  integer nL2,nFL2,nQL2
  integer i,j,k,l ! iterator
  integer,parameter::Nrun=10

  
  
  ! read reference solution 
  write (filenameFR, "(A,A)") adjustl(trim(model)),'Quasi_1_snr_800.dat'
  open(unit=8,FILE=filenameFR,STATUS='OLD',FORM='FORMATTED')
  do j=1,ncell
    ! Read from the population file
    read(8,*) Rf,(pops_ref(j,k),k=1,nlev)
  enddo
  close(8)
  
  write (filenameOut, "(A,A)") adjustl(trim(model)),'_error.dat'
  open(unit=10,FILE=filenameOut,STATUS='REPLACE',FORM='FORMATTED')
  
  
  
  do i=1,Nsnr
    ! Read Quasi Error
    write (case_name, "(A,A5)") adjustl(trim(model)),'Quasi'
    write (filenameQR, "(A,A1,I0,A5,I0,A4)") adjustl(trim(case_name)),'_',1,'_snr_',snr(i),'.dat'
    open(unit=9,FILE=filenameQR,STATUS='OLD',FORM='FORMATTED')
    do j=1,ncell
      read(9,*) Rq,(pops_QMC(j,k),k=1,nlev)
    enddo
    close(9)
    
    QL2=0.0
    QL_inf=0.0
    nQL2=0
    
    ! error of QMC to the reference
    do j=1,ncell; do k=1,nlev
        if(pops_ref(j,k)>=1e-6)then
          diff=abs(pops_QMC(j,k)-pops_ref(j,k))/pops_ref(j,k)
          QL2=QL2+diff*diff
          nQL2=nQL2+1
          if(diff>QL_inf) then
            QL_inf=diff
          endif
        endif
    enddo; enddo
    QL2=sqrt(QL2/nQL2)
    
    
    
    FL2(11)=0.0
    FL_inf(11)=0.0
    FL2(12)=0.0
    FL_inf(12)=0.0
    
    if(snr(i)<=100)then
    
        write (case_name, "(A,A6)") adjustl(trim(model)),'Pseudo'
        do l=1,Nrun 
          ! set the target file name of population data
          write (filenameFR, "(A,A1,I0,A5,I0,A4)") adjustl(trim(case_name)),'_',l,'_snr_',snr(i),'.dat'
          ! Read from the population file
          open(unit=8,FILE=filenameFR,STATUS='OLD',FORM='FORMATTED')
          do j=1,ncell
            read(8,*) Rf,(pops_FMC(j,k),k=1,nlev)
          enddo
          close(8)
          
          ! reset error
          FL2(l)=0.0
          FL_inf(l)=0.0
          nFL2=0
        
          do j=1,ncell; do k=1,nlev
              if(pops_ref(j,k)>=1e-6)then
                ! error of FMC to the reference
                diff=abs(pops_FMC(j,k)-pops_ref(j,k))/pops_ref(j,k)
                FL2(l)=FL2(l)+diff*diff
                nFL2=nFL2+1
                if(diff>FL_inf(l)) then
                  FL_inf(l)=diff
                endif
              endif
          enddo; enddo
          FL2(l)=sqrt(FL2(l)/nFL2)
          
          FL2(11)=FL2(11)+FL2(l)
          FL_inf(11)=FL_inf(11)+FL_inf(l)
          FL2(12)=FL2(12)+FL2(l)*FL2(l)
          FL_inf(12)=FL_inf(12)+FL_inf(l)*FL_inf(l)
        enddo
        FL2(11)=FL2(11)/Nrun
        FL_inf(11)=FL_inf(11)/Nrun
        FL2(12) = sqrt( FL2(12)/Nrun - FL2(11)*FL2(11) )
        FL_inf(12) = sqrt( FL_inf(12)/Nrun - FL_inf(11)*FL_inf(11) )
    
    endif
    
    write(10,*) snr(i),FL2(11),FL2(12),FL_inf(11),FL_inf(12),QL2,QL_inf
  enddo
  close(10)

end program