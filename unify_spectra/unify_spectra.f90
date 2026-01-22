! This program merges the spectra.
!   It's a modification of the original unify_stats for riblets.
! It's fed with a list of sp files and produces a single file with the sum

program unify_spectra

  use declaration
  implicit none

  integer ivar,nvar

  ! Variables to be merged
  nvar = 5
  allocate(extV(nvar))
  extV(1) = 'U'
  extV(2) = 'V'
  extV(3) = 'W'
  extV(4) = 'UV'
  extV(5) = 'P'

  ! Reads the name of files, size of the box, etc 
  call start

  ! Loop over the variables (U,V,W,UV and P)
  do ivar = 1,nvar
    istat = 0
    spU = 0d0
    spV = 0d0
    ! Variable to be merged 
    if(ivar==1.or.ivar==3.or.ivar==4.or.ivar==5)then !C! Not sure about 5
      ny0=nyu0
      nyF=nyuF
    else !C! ivar=2
      ny0=nyv0
      nyF=nyvF
    endif
    
    filout = extV(ivar)
    ! Common part of the filename
    filout = fildir(1:index(fildir,' ')-1)//'sp'//filout(1:index(filout,' ')-1)//'_'//ext1//'x'//ext2//'x'//ext3//'_'
    ! Loop over the list of files
    ! Reads a new file and add the new spectrum to the previous ones.
    do isamp = 1,nsamp
      if(ivar==1.or.ivar==3.or.ivar==4)then
	call readspXu
      elseif(ivar==2)then
	call readspXv
      elseif(ivar==5)then
	call readspXp
      endif  
    end do
    filout = extV(ivar)
    ! Save the accumulated spectrum
    if(ivar==1.or.ivar==3.or.ivar==4)then
      call writeSPu
    elseif(ivar==2)then
      call writeSPv
    elseif(ivar==5)then
      call writeSPp
    endif 
    
  end do

  write(*,*) 'DONE'

33 FORMAT(A80)

end program

subroutine readspXu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  READ SPEC  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reads the spectra

  use declaration
  implicit none
  integer distat
  integer i,k,j

  fnameima = filout(1:index(filout,' ')-1)//extL(isamp)//'.dat'

  open(20,file=fnameima,form='unformatted')
  read(20) 
  read(20) 
  read(20) 
  read(20) 
  read(20) distat
  read(20) spXu
  close(20)

  ! Iterations between consecutive outputs
  istat = istat + distat
  ! Accumulates the spectra
  do j = nyu0,nyuF
    do k = 1,nz
      do i = 0,nx
        spU(i,k,j) = spU(i,k,j) + spXu(i,k,j)
      end do
    end do
  end do

end subroutine


subroutine readspXv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  READ SPEC  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reads the spectra

  use declaration
  implicit none
  integer distat
  integer i,k,j

  fnameima = filout(1:index(filout,' ')-1)//extL(isamp)//'.dat'
  
  open(20,file=fnameima,form='unformatted')
  read(20) 
  read(20) 
  read(20) 
  read(20) 
  read(20) distat
  read(20) spXv
  close(20)

  ! Iterations between consecutive outputs
  istat = istat + distat
  ! Accumulates the spectra
  do j = nyv0,nyvF
    do k = 1,nz
      do i = 0,nx
        spV(i,k,j) = spV(i,k,j) + spXv(i,k,j)
      end do
    end do
  end do

end subroutine

subroutine readspXp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  READ SPEC  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reads the spectra

  use declaration
  implicit none
  integer distat
  integer i,k,j

  fnameima = filout(1:index(filout,' ')-1)//extL(isamp)//'.dat'

  open(20,file=fnameima,form='unformatted')
  read(20) 
  read(20) 
  read(20) 
  read(20) 
  read(20) distat
  read(20) spXp
  close(20)

  ! Iterations between consecutive outputs
  istat = istat + distat
  ! Accumulates the spectra
  do j = nyu0+1,nyuF-1
    do k = 1,nz
      do i = 0,nx
        spP(i,k,j) = spP(i,k,j) + spXp(i,k,j)
      end do
    end do
  end do

end subroutine


subroutine start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!    START    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  pi = 4d0*datan(1d0)
  allocate(dummint(88))

  ! Reads the input file
  open(40,file='unify_spectra.input',form='formatted')
  read(40,31) filelist
  read(40,20) firstgood ! Not used but still here to make the input file more general
  read(40,31) fildir    ! Directory
  read(40,20) nx
  read(40,20) nz
  read(40,20) ny0
  read(40,20) ptstilex  ! Points per tile in x, i.e. width of the tiles
  read(40,20) ptstilez  ! Points per tile in z, i.e. width of the tiles
  close(40)

  open(40,file=filelist,form='formatted')
  read(40,37) nsamp     ! Number of samples
  allocate(extL(nsamp))
  do isamp = 1,nsamp
    read(40,34) extL(isamp) ! Iteration (last number in files)
  end do
  close(40)

  write(ext1,'(i4.4)') nx
  write(ext2,'(i4.4)') nz
  write(ext3,'(i4.4)') ny0

  fnameima = fildir(1:index(fildir,' ')-1)//'spU_'//ext1//'x'//ext2//'x'//ext3//'_'//extL(1)//'.dat'
  open(20,file=fnameima,form='unformatted')
  read(20) t,Re,alp,bet,mpgx,nband,iter
  allocate(N(4,0:nband+1))
  read(20) N
  write(*,*) 'N'
  do dnz=1,4
    write(*,*) N(dnz,:)
  end do
  write(*,*) ''
  nx  = N(1,(nband+1)/2)/2
  nz  = N(2,(nband+1)/2)/2+1 !The plus 1 is due to a rouge +1 in spectra.f90...
  nyu0 = N(4,0)
  nyuF = N(4,nband)+1
  nyv0 = N(3,0)
  nyvF = N(3,nband)+1
  allocate(yu(nyu0:nyuF),dthdyu(nyu0:nyuF))
  allocate(yv(nyv0:nyvF),dthdyv(nyv0:nyvF))
  read(20) yu,dthetai,dthdyu
  read(20) yv,dthetai,dthdyv
  close(20)

  ntilex = N(1,1)/ptstilex  ! Ngal(1,1)/(points per tile in x)
  ntilez = N(2,1)/ptstilez  ! Ngal(2,1)/(points per tile in z)

  write(*,*) 'ntilez',ntilez
  write(*,*) 'nsamp',nsamp
  write(*,*) ''
  write(extC,'(i3.3)') ntilex
  write(extD,'(i3.3)') ntilez

  Lx = 2d0*pi/alp
  Lz = 2d0*pi/bet
  Ly = 2d0
  
  dnyu = N(4,1)-N(4,0)+1
  dnyv = N(3,1)-N(3,0)+1
  
  allocate(spU(0:nx,1:nz,nyu0:nyuF),spXu(0:nx,1:nz,nyu0:nyuF))
  allocate(spV(0:nx,1:nz,nyv0:nyvF),spXv(0:nx,1:nz,nyv0:nyvF))
  allocate(spP(0:nx,1:nz,nyu0+1:nyuF-1),spXp(0:nx,1:nz,nyu0+1:nyuF-1))

  istat = 0

  write(*,*) 'Lx',Lx
  write(*,*) 'Ly',Ly
  write(*,*) 'Lz',Lz
  write(*,*) 'Nx',N(1,1:nband)
  write(*,*) 'Nz',N(2,1:nband)
  write(*,*) 'Nyu',N(4,0:nband)!-N(3,0:nband-1)
  write(*,*) 'Nyv',N(3,0:nband)!-N(3,0:nband-1)
  write(*,*) ''

20 FORMAT(10X,I10)
31 FORMAT(10X,A60)
34 FORMAT(A5)
37 FORMAT(I5)

end subroutine

subroutine writeSPu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  WRITE OUT  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  fnameima='sp'//filout(1:index(filout,' ')-1)//extD//'x'//extC//'.dat'
  write(*,*) 'writing out ',fnameima(1:index(fnameima,' ')-1)

  open(10,file=fnameima,form='unformatted')
  write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
  write(10) N
  write(10) yu,dthetai,dthdyu
  write(10) yv,dthetai,dthdyv
  write(10) istat
  write(10) spU
  close(10)

end subroutine

subroutine writeSPv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  WRITE OUT  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  fnameima='sp'//filout(1:index(filout,' ')-1)//extD//'x'//extC//'.dat'
  write(*,*) 'writing out ',fnameima(1:index(fnameima,' ')-1)

  open(10,file=fnameima,form='unformatted')
  write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
  write(10) N
  write(10) yu,dthetai,dthdyu
  write(10) yv,dthetai,dthdyv
  write(10) istat
  write(10) spV
  close(10)

end subroutine

subroutine writeSPp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  WRITE OUT  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use declaration
  implicit none

  fnameima='sp'//filout(1:index(filout,' ')-1)//extD//'x'//extC//'.dat'
  write(*,*) 'writing out ',fnameima(1:index(fnameima,' ')-1)

  open(10,file=fnameima,form='unformatted')
  write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
  write(10) N
  write(10) yu,dthetai,dthdyu
  write(10) yv,dthetai,dthdyv
  write(10) istat
  write(10) spP
  close(10)

end subroutine
