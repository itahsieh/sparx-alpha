program compare_error
  implicit none
  ! model : Cosine /Shu1d
  character*32,parameter::model='pops__uniform'
  integer,parameter::nlev=21,Nres=6
  integer,parameter,dimension(Nres)::ncell=(/ 3,9,27,81,243,729 /)
  
  real*8,dimension(ncell(Nres),nlev)::pops_ref
  real*8,ALLOCATABLE,dimension(:,:)::pops
  character*32 filename,filename_tmp,filenameOut
  real*8,dimension(ncell(Nres))::R_ref
  real*8,ALLOCATABLE,dimension(:)::R_tmp 
  real*8 L2,L_inf ! L_2-norm error and L_inf-norm error
  real*8 diff ! difference of pops_FMC and pops_QMC
  integer nL2
  integer i,j,k,l ! iterator

  
  
  ! read reference solution 
  write(filename,"(A,A)") 'pops__interp','_nr729.dat'
  open(unit=8,FILE=filename,STATUS='OLD',FORM='FORMATTED')
  do j=1,ncell(Nres)
    ! Read from the population file
    read(8,*) R_ref(j),(pops_ref(j,k),k=1,nlev)
  enddo
  close(8)
  
  write(filenameOut,"(A,A)") adjustl(trim(model)),'_error.dat'
  open(unit=10,FILE=filenameOut,STATUS='REPLACE',FORM='FORMATTED')
  
  do i=1,Nres-1
    write(filename_tmp,"(A,A,I0,A)") adjustl(trim(model)),'_nr',ncell(i),'.dat'
    open(unit=9,FILE=filename_tmp,STATUS='OLD',FORM='FORMATTED')
    allocate(R_tmp(ncell(i)))
    allocate(pops(ncell(i),nlev))
    do j=1,ncell(i)
      read(9,*) R_tmp(j),(pops(j,k),k=1,nlev)
    enddo
    close(9)

    L2=0.0
    L_inf=0.0
    nL2=0
    
    ! error of QMC to the reference
    do j=1,ncell(i); do l=1,ncell(Nres)
        if ( R_tmp(j) .eq. R_ref(l) )then
            do k=1,nlev
                if ( pops_ref(l,k) >= 1e-6 ) then
                    diff = abs( pops(j,k) - pops_ref(l,k) ) / pops_ref(l,k)
                    L2 = L2 + diff * diff
                    nL2 = nL2 + 1
                    if ( diff > L_inf ) then
                            L_inf = diff
                    endif
                endif
            enddo
        endif
    enddo;enddo
    L2=sqrt(L2/nL2)
    write(10,*) ncell(i),L2,L_inf
    deallocate(R_tmp)
    deallocate(pops)
  enddo
  close(10)

end program