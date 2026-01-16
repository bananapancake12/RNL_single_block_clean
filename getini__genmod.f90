        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 16 10:54:43 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETINI__genmod
          INTERFACE 
            SUBROUTINE GETINI(U1,U2,U3,P,DIV,MYID,STATUS,IERR)
              INTEGER(KIND=4) :: MYID
              COMPLEX(KIND=8) :: U1(JLIM((1,2)):JLIM((2,2)),COLUMNS_NUM(&
     &(MYID)))
              COMPLEX(KIND=8) :: U2(JLIM((1,1)):JLIM((2,1)),COLUMNS_NUM(&
     &(MYID)))
              COMPLEX(KIND=8) :: U3(JLIM((1,2)):JLIM((2,2)),COLUMNS_NUM(&
     &(MYID)))
              COMPLEX(KIND=8) :: P(JLIM((1,3)):JLIM((2,3)),COLUMNS_NUM((&
     &MYID)))
              COMPLEX(KIND=8) :: DIV(JLIM((1,3)):JLIM((2,3)),COLUMNS_NUM&
     &((MYID)))
              INTEGER(KIND=4) :: STATUS(5)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE GETINI
          END INTERFACE 
        END MODULE GETINI__genmod
