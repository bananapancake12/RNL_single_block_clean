        !COMPILER-GENERATED INTERFACE MODULE: Fri Jan 16 10:54:43 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MBLOCK_INI_PARABOLIC_PROFILE__genmod
          INTERFACE 
            SUBROUTINE MBLOCK_INI_PARABOLIC_PROFILE(U1,U2,U3,P,MYID,    &
     &STATUS,IERR)
              COMPLEX(KIND=8), INTENT(IN) :: U1(JLIM((1,2)):,:)
              COMPLEX(KIND=8), INTENT(IN) :: U2(JLIM((1,1)):,:)
              COMPLEX(KIND=8), INTENT(IN) :: U3(JLIM((1,2)):,:)
              COMPLEX(KIND=8), INTENT(IN) :: P(:,:)
              INTEGER(KIND=4) :: MYID
              INTEGER(KIND=4) :: STATUS(5)
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE MBLOCK_INI_PARABOLIC_PROFILE
          END INTERFACE 
        END MODULE MBLOCK_INI_PARABOLIC_PROFILE__genmod
