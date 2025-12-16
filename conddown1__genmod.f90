        !COMPILER-GENERATED INTERFACE MODULE: Mon Dec  8 16:28:23 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CONDDOWN1__genmod
          INTERFACE 
            SUBROUTINE CONDDOWN1(YMC,Y2MC,X,GRID,MYID)
              USE DECLARATION
              INTEGER(KIND=4) :: MYID
              INTEGER(KIND=4) :: GRID
              REAL(KIND=8) :: YMC(DNX,DNZ,LIMPL_INCW((GRID,1,MYID))-NY((&
     &GRID,0)):LIMPL_INCW((GRID,2,MYID))-NY((GRID,0)))
              REAL(KIND=8) :: Y2MC(DNX,DNZ,LIMPL_INCW((GRID,1,MYID))-NY(&
     &(GRID,0)):LIMPL_INCW((GRID,2,MYID))-NY((GRID,0)))
              REAL(KIND=8) :: X(IGAL,KGAL,LIMPL_INCW((GRID,1,MYID)):    &
     &LIMPL_INCW((GRID,2,MYID)))
            END SUBROUTINE CONDDOWN1
          END INTERFACE 
        END MODULE CONDDOWN1__genmod
