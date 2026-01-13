module transpose
  use declaration
  implicit none
contains

! transposes that need to be compilled first as its called by other modules


  subroutine modes_to_planes_UVP (xPL,x,grid,nygrid,nygrid_LB,myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! MODES TO PLANES !!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use declaration
    implicit none

    include 'mpif.h'             ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr,myid

    integer i,k,j,jminS,jmaxS,jminR,jmaxR,dki,grid
    integer column
    integer inode,yourid
    integer msizeR,msizeS
    complex(8), intent(in) :: x(jlim(1,grid):,:)
    real(8)      xPL(igal,kgal,jgal(grid,1)-1:jgal(grid,2)+1)
    integer, intent(in) :: nygrid, nygrid_LB
    complex(8), allocatable :: buffS(:,:),buffR(:,:)

    ! Loop for itself
    ! Transpose the cube that it already owns

    yourid = myid

    jminR = max(planelim(grid,1,  myid),jlim(1,grid)+1)
    jmaxR = min(planelim(grid,2,  myid),jlim(2,grid)-1)

    if (jminR==nygrid_LB+1 .and. jmaxR>=jminR) then
      jminR = jminR-1
    end if
    if (jmaxR==nygrid   .and. jmaxR>=jminR) then
      jmaxR = jmaxR+1
    end if
    do j = jminR,jmaxR
      do column = 1,columns_num(yourid)
        i = columns_i(column,yourid)
        k = columns_k(column,yourid) - dk(column,yourid)
        xPL(2*i+1,k,j) = dreal(x(j,column))
        xPL(2*i+2,k,j) = dimag(x(j,column))
      end do
    end do


    do inode = 1,pnodes-1
      yourid = ieor(myid,inode)
      if (yourid<np) then

        jminS = max(planelim(grid,1,yourid),jlim(1,grid)+1)
        jmaxS = min(planelim(grid,2,yourid),jlim(2,grid)-1)
        jminR = max(planelim(grid,1,  myid),jlim(1,grid)+1)
        jmaxR = min(planelim(grid,2,  myid),jlim(2,grid)-1)
        if (jminS==nygrid_LB+1  ) then
          jminS = jminS-1
        end if
        if (jmaxS==nygrid) then
          jmaxS = jmaxS+1
        end if
        if (jminR==nygrid_LB+1  ) then
          jminR = jminR-1
        end if
        if (jmaxR==nygrid) then
          jmaxR = jmaxR+1
        end if
        allocate(buffS(jminS:jmaxS,columns_num(  myid)))
        allocate(buffR(jminR:jmaxR,columns_num(yourid)))
        msizeS = 2*(columns_num(  myid)*(jmaxS-jminS+1))  ! 2 times because it's complex
        msizeR = 2*(columns_num(yourid)*(jmaxR-jminR+1))
        msizeS = max(msizeS,0)
        msizeR = max(msizeR,0)

        do j = jminS,jmaxS
          do column = 1,columns_num(myid)
            buffS(j,column) = x(j,column)

          end do
        end do

        call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid, &
  &                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid, &
  &                         MPI_COMM_WORLD,status,ierr)
        do j = jminR,jmaxR
          do column = 1,columns_num(yourid)
            i = columns_i(column,yourid)
            k = columns_k(column,yourid) - dk(column,yourid)
            xPL(2*i+1,k,j) = dreal(buffR(j,column))
            xPL(2*i+2,k,j) = dimag(buffR(j,column))
          end do
        end do

        deallocate(buffR,buffS)

        ! end do
      end if
    end do

  end subroutine

  subroutine planes_to_modes_UVP (x,xPL,grid,nygrid,nygrid_LB,myid,status,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!! PLANES TO MODES  NEW !!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Prepare the vectors for the Fourier transform
  ! The procs broadcast the data they have of a plane, and receive the data of a pencil,
  !  they also transpose the data for the Fourier transform XZY -> YXZ

    use declaration
    implicit none

    include 'mpif.h'             ! MPI variables
    integer status(MPI_STATUS_SIZE),ierr,myid

    real :: t_start, t_end

    integer :: i,k,j,jminS,jmaxS,jminR,jmaxR,grid
    integer :: column, nprocs
    integer inode,yourid
    integer msizeR,msizeS
    complex(8), intent(inout) :: x(jlim(1,grid):,:) 
    real(8)      xPL(igal,kgal,jgal(grid,1)-1:jgal(grid,2)+1)
    complex(8), allocatable:: buffS(:,:),buffR(:,:)
    integer, intent(in) :: nygrid,nygrid_LB

    ! write(6,*) "starting self transpose"

    ! Loop for itself
    ! Transpose the cube that it already owns


    ! call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    ! write(6,*) "DEBUG myid=", myid, "MPI nprocs=", nprocs, "np var=", np
    ! call flush(6)



    jminR = max(limPL_excw(grid,1,myid),jlim(1,grid)+1)  ! Select the planes to transpose 
    jmaxR = min(limPL_excw(grid,2,myid),jlim(2,grid)-1)

    ! write(*,*) "rank", myid, "x bounds:", lbound(x,1), ubound(x,1), "cols:", lbound(x,2), ubound(x,2)
    

    
    if (jminR==nygrid_LB+1 .and. jmaxR>=jminR) then   ! Special cases: walls
      jminR = jminR-1
    end if
    if (jmaxR==nygrid   .and. jmaxR>=jminR) then
      jmaxR = jmaxR+1
    end if
    ! write(*,*) "rank", myid, "jminR/jmaxR:", jminR, jmaxR, "grid", grid

    do j = jminR,jmaxR
      do column = 1,columns_num(myid)
        i = columns_i(column,myid)
        k = columns_k(column,myid) - dk(column,myid)
        x(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j)) ! Transposition: Reordering from XZY to YC
      end do
    end do

    !end do

    ! write(6,*) "finished self transpose"

    do inode = 1,pnodes-1
      yourid = ieor(myid,inode)   ! XOR. It's used to pair procs 1-to-1
      if (yourid<np) then
          jminS = max(limPL_excw(grid,1,  myid),jlim(1,grid)+1)  ! Select the planes to be SENT.
          jmaxS = min(limPL_excw(grid,2,  myid),jlim(2,grid)-1)  ! max and min because maybe this proc needs less planes that the other proc has
          jminR = max(limPL_excw(grid,1,yourid),jlim(1,grid)+1)  ! Select the planes to be RECEIVED
          jmaxR = min(limPL_excw(grid,2,yourid),jlim(2,grid)-1)

          ! Adding the walls =

          if (jminS==nygrid_LB+1  ) then
            jminS=jminS-1
          end if
          if (jmaxS==nygrid) then
            jmaxS=jmaxS+1
          end if
          if (jminR==nygrid_LB+1  ) then
            jminR=jminR-1
          end if
          if (jmaxR==nygrid) then
            jmaxR=jmaxR+1
          end if
          allocate(buffS(jminS:jmaxS,columns_num(yourid)))
          allocate(buffR(jminR:jmaxR,columns_num(  myid)))

          ! if (myid ==0 .and. yourid == 4 .and. jband == 2) then
          !   write(6,*) "buffs", size(buffS,1), size(buffS,2)
          ! end if 

          msizeS = 2*(columns_num(yourid)*(jmaxS-jminS+1))     ! Size of the data to be SENDER (times 2, because it is complex)
          msizeR = 2*(columns_num(  myid)*(jmaxR-jminR+1))     ! Size of the data to be RECEIVED
          msizeS = max(msizeS,0)                                     ! The size has to be 0 or positive. 
          msizeR = max(msizeR,0)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do j=jminS,jmaxS
            do column = 1,columns_num(yourid)
              i = columns_i(column,yourid)
              k = columns_k(column,yourid) - dk(column,yourid)
              buffS(j,column) = dcmplx(xPL(2*i+1,k,j),xPL(2*i+2,k,j))     ! The data is transposed and stored in a buffer
            end do
          end do

          ! call cpu_time(t_start)
          
          call MPI_SENDRECV(buffS,msizeS,MPI_REAL8,yourid,77*yourid+53*myid, &   ! SEND_RECV so it can send and receive at the same time
  &                         buffR,msizeR,MPI_REAL8,yourid,53*yourid+77*myid, &
  &                         MPI_COMM_WORLD,status,ierr)

          ! call cpu_time(t_end)

          ! if (myid ==0 .and. yourid == 4 .and. jband == 2) then
          !   write(6,*) "cpu time", t_start, t_end
          ! end if

          do j=jminR,jmaxR
            do column = 1,columns_num(myid)
              x(j,column) = buffR(j,column)                         ! Store the data received
            end do
          end do
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          deallocate(buffR,buffS)
        ! end do
      end if
    end do

  end subroutine
end module transpose