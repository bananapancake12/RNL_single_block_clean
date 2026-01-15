        !COMPILER-GENERATED INTERFACE MODULE: Thu Jan 15 17:25:35 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RECVX__genmod
          INTERFACE 
            SUBROUTINE RECVX(X,XM,GRID,NYGRID,NYGRID_LB,STATUS,IERR)
              INTEGER(KIND=4) :: NYGRID_LB
              INTEGER(KIND=4) :: NYGRID
              INTEGER(KIND=4) :: GRID
              REAL(KIND=8) :: X(NYGRID_LB:NYGRID+1)
              REAL(KIND=8) :: XM(LIMPL_INCW((GRID,1,0)):LIMPL_INCW((GRID&
     &,2,0)))
              INTEGER(KIND=4) :: STATUS(5)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE RECVX
          END INTERFACE 
        END MODULE RECVX__genmod
