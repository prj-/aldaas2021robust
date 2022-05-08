!
!  Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
!       Date: 2021-07-11
!
!  Copyright (C) 2021-     Centre National de la Recherche Scientifique
!
!  This script is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!  If you use this script, you are kindly asked to cite the following article:
!
!  "A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations",
!  H. Al Daas, P. Jolivet, and J. A. Scott (2022).
!

  program main
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
  Vec                            x,b,c
  Mat ::                         A,Bf,aux,Cf = PETSC_NULL_MAT,perm,Neumann(2)
  MatPartitioning                mpart
  KSP                            ksp
  PC                             pc
  KSPType                        ktype
  PCType                         ptype
  PCHPDDMCoarseCorrectionType    ctype
  IS                             is,rows,cols(1)
  PetscMPIInt                    rank
  PetscReal                      norm(2)
  PetscInt                       m,n
  PetscBool ::                   normal = PETSC_FALSE,flg
  PetscViewer                    viewer
  character*(PETSC_MAX_PATH_LEN) name
  PetscErrorCode                 ierr

  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
  ! begin loading phase (section 4.1.1)
  PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
  PetscCallA(PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-mat_name',name,flg,ierr))
  if (.not. flg) then
    print *,'Must specify matrix name with -mat_name'
    stop
  endif
  PetscCallA(PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_READ,viewer,ierr))
  PetscCallA(MatLoad(A,viewer,ierr))
  PetscCallA(PetscViewerDestroy(viewer,ierr))
  ! end loading phase
  PetscCallA(KSPCreate(PETSC_COMM_WORLD,ksp,ierr))
  PetscCallA(KSPSetFromOptions(ksp,ierr))
  PetscCallA(KSPGetPC(ksp,pc,ierr))
  PetscCallA(PCGetType(pc,ptype,ierr))
  if (ptype == PCQR) then
    PetscCallA(KSPSetOperators(ksp,A,A,ierr))
    PetscCallA(KSPSetType(ksp,KSPPREONLY,ierr))
  else
    ! begin partitioning phase (section 4.1.1)
    PetscCallA(MatTransposeMatMult(A,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,Bf,ierr))
    PetscCallA(MatPartitioningCreate(PETSC_COMM_WORLD,mpart,ierr))
    PetscCallA(MatPartitioningSetAdjacency(mpart,Bf,ierr))
    PetscCallA(MatPartitioningSetFromOptions(mpart,ierr))
    PetscCallA(MatPartitioningApply(mpart,is,ierr))
    PetscCallA(MatPartitioningDestroy(mpart,ierr))
    PetscCallA(ISBuildTwoSided(is,PETSC_NULL_IS,rows,ierr))
    PetscCallA(ISDestroy(is,ierr))
    PetscCallA(MatCreateSubMatrix(Bf,rows,rows,MAT_INITIAL_MATRIX,perm,ierr))
    PetscCallA(MatDestroy(Bf,ierr))
    Bf = perm
    PetscCallA(MatGetOwnershipRange(A,m,n,ierr))
    PetscCallA(ISCreateStride(PETSC_COMM_WORLD,n-m,m,1,is,ierr))
    PetscCallA(MatCreateSubMatrix(A,is,rows,MAT_INITIAL_MATRIX,perm,ierr))
    PetscCallA(MatDestroy(A,ierr))
    PetscCallA(ISDestroy(is,ierr))
    PetscCallA(ISDestroy(rows,ierr))
    A = perm
    ! end partitioning phase
    ! begin setup phase (section 4.1.2)
    PetscCallA(MatNorm(Bf,NORM_FROBENIUS,norm(1),ierr))
    PetscCallA(MatSetOption(Bf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr))
    PetscCallA(MatSetOption(Bf,MAT_SYMMETRIC,PETSC_TRUE,ierr))
    PetscCallA(MatShift(Bf,norm(1) * 1.0e-10,ierr))
    PetscCallA(KSPGetType(ksp,ktype,ierr))
    normal = ktype /= KSPLSQR
    if (normal) then
      PetscCallA(MatCreateNormal(A,Cf,ierr))
    endif
    PetscCallA(KSPSetOperators(ksp,merge(Cf,A,normal),Bf,ierr))
    if (ptype == PCHPDDM) then
      PetscCallA(MatGetOwnershipRange(Bf,m,n,ierr))
      PetscCallA(ISCreateStride(PETSC_COMM_SELF,n-m,m,1,cols(1),ierr))
      PetscCallA(MatGetSize(A,m,PETSC_NULL_INTEGER,ierr))
      PetscCallA(ISCreateStride(PETSC_COMM_SELF,m,0,1,rows,ierr))
      PetscCallA(MatSetOption(A,MAT_SUBMAT_SINGLEIS,PETSC_TRUE,ierr))
      PetscCallA(MatCreateSubMatrices(A,1,rows,cols(1),MAT_INITIAL_MATRIX,Neumann,ierr))
      PetscCallA(MatFindZeroRows(Neumann(1),is,ierr))
      PetscCallA(MatDestroySubMatrices(1,Neumann,ierr))
      PetscCallA(MatIncreaseOverlap(Bf,1,cols,1,ierr))
      PetscCallA(MatCreateSubMatrices(A,1,rows,cols(1),MAT_INITIAL_MATRIX,Neumann,ierr))
      PetscCallA(ISDestroy(rows,ierr))
      PetscCallA(MatZeroRowsIS(Neumann(1),is,0.0D0,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr))
      PetscCallA(ISDestroy(is,ierr))
      PetscCallA(MatTransposeMatMult(Neumann(1),Neumann(1),MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,aux,ierr))
      PetscCallA(MatDestroySubMatrices(1,Neumann,ierr))
      PetscCallA(MatNorm(aux,NORM_FROBENIUS,norm(1),ierr))
      PetscCallA(MatSetOption(aux,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr))
      PetscCallA(MatShift(aux,norm(1) * 1.0e-8,ierr))
      PetscCallA(PCHPDDMSetAuxiliaryMat(pc,cols(1),aux,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ierr))
      PetscCallA(PCHPDDMHasNeumannMat(pc,PETSC_TRUE,ierr))
      if (.not. normal) then
        PetscCallA(PCHPDDMGetCoarseCorrectionType(pc,ctype,ierr))
        if (ctype == PC_HPDDM_COARSE_CORRECTION_DEFLATED) then
          PetscCallA(PCHPDDMSetCoarseCorrectionType(pc,PC_HPDDM_COARSE_CORRECTION_BALANCED,ierr))
        endif
      endif
      PetscCallA(ISDestroy(cols(1),ierr))
      PetscCallA(MatDestroy(aux,ierr))
    endif
    ! end setup phase
    PetscCallA(MatDestroy(Bf,ierr))
    PetscCallA(MatDestroy(Cf,ierr))
  endif
  ! begin solution phase (section 4.1.3)
  PetscCallA(MatCreateVecs(A,x,b,ierr))
  PetscCallA(VecSetRandom(b,PETSC_NULL_RANDOM,ierr))
  PetscCallA(VecDuplicate(x,c,ierr))
  PetscCallA(MatMultTranspose(A,b,c,ierr))
  PetscCallA(KSPSolve(ksp,merge(c,b,normal),x,ierr))
  ! end solution phase
  PetscCallA(VecScale(b,-1.0D0,ierr))
  PetscCallA(MatMultAdd(A,x,b,b,ierr))
  PetscCallA(VecNorm(b,NORM_2,norm(1),ierr))
  PetscCallA(MatMultTranspose(A,b,c,ierr))
  PetscCallA(VecNorm(c,NORM_2,norm(2),ierr))
  if (rank .eq. 0) then
    PetscCallA(PetscPrintf(PETSC_COMM_SELF,"||A^T(Ax-b)|| / ||Ax-b|| = ",ierr))
    write(name,*) norm(2)
    PetscCallA(PetscPrintf(PETSC_COMM_SELF,trim(name)//" / ",ierr))
    write(name,*) norm(1)
    PetscCallA(PetscPrintf(PETSC_COMM_SELF,trim(name)//" = ",ierr))
    norm(1) = norm(2) / norm(1)
    write(name,*) norm(1)
    PetscCallA(PetscPrintf(PETSC_COMM_SELF,trim(name)//"\n",ierr))
  endif
  PetscCallA(KSPDestroy(ksp,ierr))
  PetscCallA(VecDestroy(c,ierr))
  PetscCallA(VecDestroy(x,ierr))
  PetscCallA(VecDestroy(b,ierr))
  PetscCallA(MatDestroy(A,ierr))
  PetscCallA(PetscFinalize(ierr))
  end
