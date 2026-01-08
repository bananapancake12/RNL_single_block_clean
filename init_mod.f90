
module init_mod
  use declaration
  use transpose
  implicit none
contains

  subroutine getini(u1,u2,u3,p,div,myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!     GETINI  NEW   !!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Grab the initial condition and initialize some variables

    use declaration
    use littleharsh_mod

    implicit none

    include 'mpif.h'                                  ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables
    ! type(cfield)  u1,u2,u3
    ! type(cfield)  p 
    complex(8), intent(inout) :: u1(jlim(1,ugrid):,:), u2(jlim(1,vgrid):,:), u3(jlim(1,ugrid):,:), p(:,:)
    type(cfield)  div

    ! write(*,*) "call shape(u1)=", shape(u1), "lb=", lbound(u1), "ub=", ubound(u1)


    if (myid==0) then
      write(*,*) 'Launching...'
    end if

    if (flag_init==1) then       ! Initial conditions borrowed from another 
      if (myid==0) then
        write(*,*) 'starting from multi-block, flat channel?'
        write(*,*) 'start file: ',trim(dirin)//trim(fnameimb)
      end if
      call mblock_ini(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
      if (myid==0) then
        ! call flowrateIm(Qx,u1(nyu_LB,1))        !check why flowrate used midband ! Column 1 of proc 0 is mode (0,1) [the meeeean]
        call flowrateIm(Qx,u1(:,1))
        !write(6,*) "flowrateIm", Qx
        if (flag_ctpress/=1) then
          QxT = Qx
        end if
      end if
      iter0 = 0
      mpgz  = 0d0
      t     = 0d0
    else if (flag_init==2) then  ! Continuing simulation
      if (myid==0) then
        write(*,*) 'continuing simulation'
        write(*,*) 'start file:',fnameimb
      end if
      call mblock_ini(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
      if (myid==0) then
        !call flowrateIm(Qx,u1(nyu_LB,1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
        call flowrateIm(Qx,u1(:,1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
        ! TODO change this condition for 'flag_ctpress == 0'
        if (flag_ctpress/=1) then  ! Constant flow rate (flag_ctpress == 1, constant pressure gradient)
          QxT = Qx
        end if
      end if
      mpgz = 0d0
    else if (flag_init==3) then  ! Parabolic profile
      if (myid==0) then
        write(*,*) 'parabolic profile'
      end if
      call mblock_ini_parabolic_profile(u1,u2,u3,p,myid,status,ierr)         ! Initializes u1, u2, u3, p, u1PL, u2PL, u3PL, ppPL
      if (myid==0) then
        call flowrateIm(Qx,u1(:,1))        ! Column 1 of proc 0 is mode (0,1) [the meeeean]
        ! TODO change this condition for 'flag_ctpress == 0'
        if (flag_ctpress/=1) then  ! Constant flow rate (flag_ctpress == 1, constant pressure gradient)
          QxT = Qx
        end if
      end if
      iter0 = 0
      mpgz  = 0d0
      t     = 0d0
    else 
      write(*,*) 'INITIAL CONDITIONS NOT IMPLEMENTED YET'
      call MPI_FINALIZE(ierr)
      stop
    end if

    ! write(6,*) "check1 ", myid

    ! Broadcast the time step and time.
    ! If it's initialized from a previous simulation this value is already known, otherwise it's set to 0
    call MPI_BCAST(iter0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(t    ,1,MPI_REAL8  ,0,MPI_COMM_WORLD,ierr)
    iter   = iter0
    !  iter0=iter-nstat
    iwrite = iter

    ! 'Probably' this is used to initialize the divergence in the case of a new simulation
    ! write(6,*) "call divergence ", myid
    call divergence(div%f,u1,u2,u3,myid)

    ! write(6,*) "call init stats", myid
    call init_stats(myid)
    ! write(6,*) "finished init_stats", myid
    
    ! write(6,*) "call init_sl stats", myid
    ! call init_sl_stats(myid)

    if (myid==0) then
      write(*,*) 'Lx    ',Lx
      write(*,*) 'Ly    ',Ly
      write(*,*) 'Lz    ',Lz
      write(*,*) 'Nx    ',Nspec_x
      write(*,*) 'Nz    ',Nspec_z
      write(*,*) 'Nyv   ',Nyv
      write(*,*) 'Nyu   ',nyu
      write(*,*) 'Ngalx ',Ngal_x
      write(*,*) 'Ngalz ',Ngal_z
      write(*,*) ''
      write(*,*) 'dymin ',yu(1)-yu(0)
      write(*,*) 'dymax ',yu((nyu+1)/2+1)-yu((nyu+1)/2)
      write(*,*) ''
      write(*,*) 'Re ',Re
      write(*,*) 't  ',t
    end if

  end subroutine

  subroutine mblock_ini(u1,u2,u3,p,myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!   MBLOCK INI   !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none

    include 'mpif.h'                                  ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr,myid         ! MPI variables

    integer j,i,k,column
    ! real(8) sigmaz,z,sigmax,x,fact
    complex(8), intent(inout) :: u1(:,:), u2(:,:), u3(:,:), p(:,:)

    u1PL = 0d0
    u2PL = 0d0
    u3PL = 0d0
    ppPL = 0d0

    ! write(*,*) "shape(u1)=", shape(u1), "lb=", lbound(u1), "ub=", ubound(u1)



    u1 = 0d0
    u2 = 0d0
    u3 = 0d0
    p = 0d0

    write(6,*) " Calling read in"
    call read_in(myid)
    
    ! write(6,*) " start planes to modes"
    call planes_to_modes_UVP(u1,u1PL,2,nyu,nyu_LB,myid,status,ierr)
    call planes_to_modes_UVP(u2,u2PL,1,nyv,nyv_LB,myid,status,ierr)
    call planes_to_modes_UVP(u3,u3PL,2,nyu,nyu_LB,myid,status,ierr)
    call planes_to_modes_UVP(p ,ppPL,3,nyp,nyp_LB,myid,status,ierr)
    ! write(6,*) " finished pplanes to modes"

    ! if(myid ==0) then 
    !   do j= 1, 22
    !     write(6,*) p%f(j,1)
    !   end do
    ! end if 

  end subroutine
end module init_mod