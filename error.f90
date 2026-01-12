module error_mod
  use declaration
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
    complex(8), intent(in) :: A(:,:)
    real(8) erri,errband

    erri = 0d0
    errband = 0d0
    do column = 1,columns_num(myid)
        do j = jlim(1,pgrid),jlim(2,pgrid)
        err     = abs(A(j,column))
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


end module error_mod