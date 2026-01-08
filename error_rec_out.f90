module error_rec_out
  use declaration
  use transpose
  implicit none

  contains

    subroutine error(A,myid,ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!     ERROR      !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none

    include 'mpif.h'             ! MPI variables
    integer ierr

    integer j,column,myid
    type(cfield) A
    real(8) erri,errband

    erri = 0d0
    errband = 0d0
    do column = 1,columns_num(myid)
        do j = jlim(1,pgrid),jlim(2,pgrid)
        err     = abs(A%f(j,column))
        errband = max(errband,err)

        ! write(6,*) "errband", errband
        end do
    end do
    erri = max(erri,errband)


    !if (myid == 0) then
    ! write(6,*) 'error(): size(A%f,1:2)=', size(A%f,1), size(A%f,2)
    ! write(6,*) 'error(): jlim(1:2,pgrid)=', jlim(1,pgrid), jlim(2,pgrid)
    ! write(6,*) 'error(): columns_num(myid)=', columns_num(myid), ' myid=', myid
    !end if


    call MPI_ALLREDUCE(erri,err,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)



    end subroutine

    subroutine record_out(u1,myid)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!   RECORD OUT   !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    ! use littleharsh_mod
    implicit none

    include 'mpif.h'             ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr,myid

    ! type(cfield) u1
    complex(8), intent(in) :: u1(jlim(1,ugrid):,:)
    integer nx,nz, i 
    integer j,jmax,iproc
    real(8), allocatable:: buffSR(:,:)
    integer, allocatable:: dummint(:)
    real(8) Uslip  

    ! Rebuilding N 
    if (.not. allocated(N)) allocate(N(4,0:4))


    N= 0 
    N(1,1:3) = Nspec_x
    N(1,4) = -2
    N(2,1:3) = Nspec_z
    N(3,3) = nyv
    N(4,3) = nyu

    if (myid == 0) then
        write(6,*) "N:"
        do i = 1,4
        write(6,*) N(i,0:4)
        end do
    end if 

    if (myid/=0) then
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
        call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,123*myid,MPI_COMM_WORLD,ierr)
        end do
        do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
        call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,124*myid,MPI_COMM_WORLD,ierr)
        end do
        do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
        call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,125*myid,MPI_COMM_WORLD,ierr)
        end do
        do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
        call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
        call MPI_SEND(buffSR,nx*nz,MPI_REAL8,0,126*myid,MPI_COMM_WORLD,ierr)
        end do

        deallocate(buffSR)
    else
        write(ext4,'(i5.5)') int(10d0*(t))!int(t)!
        allocate(dummint(88))
        dummint = 0
        !!!!!!!!!!!!!    u1    !!!!!!!!!!!!!
        fnameima = 'output/u1_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
        open(10,file=fnameima,form='unformatted')
        write(10) t,Re,alp,bet,mpgx,nband,iter,dummint 
        write(10) N
        write(10) yu,dthetavi,dthdyu
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
        call u_to_buff(buffSR,u1PL(1,1,j),nx,nz,igal,kgal)
        write(10) j,1,nx,nz,yu(j),buffSR
        end do
        deallocate(buffSR)
        do iproc = 1,np-1
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
            call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,123*iproc,MPI_COMM_WORLD,status,ierr)
            write(10) j,1,nx,nz,yu(j),buffSR
        end do
        deallocate(buffSR)
        end do
        close(10)
        !!!!!!!!!!!!!    u2    !!!!!!!!!!!!!
        fnameima='output/u2_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
        open(10,file=fnameima,form='unformatted')
        write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
        write(10) N
        write(10) yv,dthetavi,dthdyv
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(vgrid,1,myid),limPL_incw(vgrid,2,myid)
        call u_to_buff(buffSR,u2PL(1,1,j),nx,nz,igal,kgal)
        write(10) j,2,nx,nz,yv(j),buffSR
        end do
        deallocate(buffSR)
        do iproc = 1,np-1
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(vgrid,1,iproc),limPL_incw(vgrid,2,iproc)
            call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,124*iproc,MPI_COMM_WORLD,status,ierr)
            write(10) j,2,nx,nz,yv(j),buffSR
        end do
        deallocate(buffSR)
        end do
        close(10)
        !!!!!!!!!!!!!    u3    !!!!!!!!!!!!!
        fnameima = 'output/u3_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
        open(10,file=fnameima,form='unformatted')
        write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
        write(10) N
        write(10) yu,dthetavi,dthdyu
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(ugrid,1,myid),limPL_incw(ugrid,2,myid)
        call u_to_buff(buffSR,u3PL(1,1,j),nx,nz,igal,kgal)
        write(10) j,3,nx,nz,yu(j),buffSR
        end do
        deallocate(buffSR)
        do iproc = 1,np-1
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(ugrid,1,iproc),limPL_incw(ugrid,2,iproc)
            call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,125*iproc,MPI_COMM_WORLD,status,ierr)
            write(10) j,3,nx,nz,yu(j),buffSR
        end do
        deallocate(buffSR)
        end do
        close(10)
        !!!!!!!!!!!!!    p     !!!!!!!!!!!!!
        fnameima = 'output/p_'//ext1//'x'//ext2//'x'//ext3//'_t'//ext4//'.dat'
        open(10,file=fnameima,form='unformatted')
        write(10) t,Re,alp,bet,mpgx,nband,iter,dummint
        write(10) N
        write(10) yu,dthetavi,dthdyu
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(pgrid,1,myid),limPL_incw(pgrid,2,myid)
        call u_to_buff(buffSR,ppPL(1,1,j),nx,nz,igal,kgal)
        write(10) j,4,nx,nz,yu(j),buffSR
        end do
        deallocate(buffSR)
        do iproc = 1,np-1
        nx = Nspec_x+2
        nz = Nspec_z
        allocate(buffSR(nx,nz))
        do j = limPL_incw(pgrid,1,iproc),limPL_incw(pgrid,2,iproc)
            call MPI_RECV(buffSR,nx*nz,MPI_REAL8,iproc,126*iproc,MPI_COMM_WORLD,status,ierr)
            write(10) j,4,nx,nz,yu(j),buffSR
        end do
        deallocate(buffSR)
        end do
        close(10)
        deallocate(dummint)
    end if

    if (myid==0) then

        Uslip = ((u1(1,1))*(-1d0-yu(0))+(u1(0,1))*(yu(1)+1d0))/(yu(1)-yu(0))

        ! call flowrateIm(Qx,u1(nyu_LB,1))
        call flowrateIm(Qx,u1(:,1))
        ! call maxvel(u1(nyu_LB,1))
        call maxvel(u1(:,1))
        write(*,*) ''
        write(*,*) 'iter',iter
        write(*,*) 't   ',t
        write(*,*) 'dtv ',dtv
        write(*,*) 'dtc ',dtc
        write(*,*) 'dt  ',dt
        write(*,*) 'err ',err
        write(*,*) 'Qx  ',Qx
        if (flag_ctpress==0) then
        write(*,*) 'QxT ',QxT
        write(*,*) 'mpgx',mpgx
        write(*,*) 'dpgx',dgx
        else
        write(*,*) 'mpgx',mpgx
        end if
        write(*,*) 'Umax',Umax
        write(*,*) 'Uslp',Uslip
    !    write(*,*) 'utau',utau
    ! Save to history file
        if (flag_ctpress==0) then
        write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,QxT,mpgx,dgx,Umax,Uslip
        else
        write(30) flag_ctpress,iter,t,dtv,dtc,dt,err,Qx,    mpgx,    Umax,Uslip
        end if
        flush(30)
    end if

    end subroutine

end module error_rec_out