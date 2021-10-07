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
!  H. Al Daas, P. Jolivet, and J. A. Scott (2021).
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

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print *,'Unable to initialize PETSc'
    stop
  endif
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
  ! begin loading phase (section 4.1.1)
  call MatCreate(PETSC_COMM_WORLD,A,ierr);CHKERRA(ierr)
  call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-mat_name',name,flg,ierr);CHKERRA(ierr)
  if (.not. flg) then
    print *,'Must specify matrix name with -mat_name'
    stop
  endif
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_READ,viewer,ierr);CHKERRA(ierr)
  call MatLoad(A,viewer,ierr);CHKERRA(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRA(ierr)
  ! end loading phase
  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr);CHKERRA(ierr)
  call KSPSetFromOptions(ksp,ierr);CHKERRA(ierr)
  call KSPGetPC(ksp,pc,ierr);CHKERRA(ierr)
  call PCGetType(pc,ptype,ierr);CHKERRA(ierr)
  if (ptype == PCQR) then
    call KSPSetOperators(ksp,A,A,ierr);CHKERRA(ierr)
    call KSPSetType(ksp,KSPPREONLY,ierr);CHKERRA(ierr)
  else
    ! begin partitioning phase (section 4.1.1)
    call MatTransposeMatMult(A,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,Bf,ierr);CHKERRA(ierr)
    call MatPartitioningCreate(PETSC_COMM_WORLD,mpart,ierr);CHKERRA(ierr)
    call MatPartitioningSetAdjacency(mpart,Bf,ierr);CHKERRA(ierr)
    call MatPartitioningSetFromOptions(mpart,ierr);CHKERRA(ierr)
    call MatPartitioningApply(mpart,is,ierr);CHKERRA(ierr)
    call MatPartitioningDestroy(mpart,ierr);CHKERRA(ierr)
    call ISBuildTwoSided(is,PETSC_NULL_IS,rows,ierr);CHKERRA(ierr)
    call ISDestroy(is,ierr);CHKERRA(ierr)
    call MatCreateSubMatrix(Bf,rows,rows,MAT_INITIAL_MATRIX,perm,ierr);CHKERRA(ierr)
    call MatDestroy(Bf,ierr);CHKERRA(ierr)
    Bf = perm
    call MatGetOwnershipRange(A,m,n,ierr);CHKERRA(ierr)
    call ISCreateStride(PETSC_COMM_WORLD,n-m,m,1,is,ierr);CHKERRA(ierr)
    call MatCreateSubMatrix(A,is,rows,MAT_INITIAL_MATRIX,perm,ierr);CHKERRA(ierr)
    call MatDestroy(A,ierr);CHKERRA(ierr)
    call ISDestroy(is,ierr);CHKERRA(ierr)
    call ISDestroy(rows,ierr);CHKERRA(ierr)
    A = perm
    ! end partitioning phase
    ! begin setup phase (section 4.1.2)
    call MatNorm(Bf,NORM_FROBENIUS,norm(1),ierr);CHKERRA(ierr)
    call MatSetOption(Bf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr);CHKERRA(ierr)
    call MatSetOption(Bf,MAT_SYMMETRIC,PETSC_TRUE,ierr);CHKERRA(ierr)
    call MatShift(Bf,norm(1) * 1.0e-10,ierr);CHKERRA(ierr)
    call KSPGetType(ksp,ktype,ierr);CHKERRA(ierr)
    normal = ktype /= KSPLSQR
    if (normal) then
      call MatCreateNormal(A,Cf,ierr);CHKERRA(ierr)
    endif
    call KSPSetOperators(ksp,merge(Cf,A,normal),Bf,ierr);CHKERRA(ierr)
    if (ptype == PCHPDDM) then
      call MatGetOwnershipRange(Bf,m,n,ierr);CHKERRA(ierr)
      call ISCreateStride(PETSC_COMM_SELF,n-m,m,1,cols(1),ierr);CHKERRA(ierr)
      call MatGetSize(A,m,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      call ISCreateStride(PETSC_COMM_SELF,m,0,1,rows,ierr);CHKERRA(ierr)
      call MatSetOption(A,MAT_SUBMAT_SINGLEIS,PETSC_TRUE,ierr);CHKERRA(ierr)
      call MatCreateSubMatrices(A,1,rows,cols(1),MAT_INITIAL_MATRIX,Neumann,ierr);CHKERRA(ierr)
      call MatFindZeroRows(Neumann(1),is,ierr);CHKERRA(ierr)
      call MatDestroySubMatrices(1,Neumann,ierr);CHKERRA(ierr)
      call MatIncreaseOverlap(Bf,1,cols,1,ierr);CHKERRA(ierr)
      call MatCreateSubMatrices(A,1,rows,cols(1),MAT_INITIAL_MATRIX,Neumann,ierr);CHKERRA(ierr)
      call ISDestroy(rows,ierr);CHKERRA(ierr)
      call MatZeroRowsIS(Neumann(1),is,0.0D0,PETSC_NULL_VEC,PETSC_NULL_VEC,ierr);CHKERRA(ierr)
      call ISDestroy(is,ierr);CHKERRA(ierr)
      call MatTransposeMatMult(Neumann(1),Neumann(1),MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,aux,ierr);CHKERRA(ierr)
      call MatDestroySubMatrices(1,Neumann,ierr);CHKERRA(ierr)
      call MatNorm(aux,NORM_FROBENIUS,norm(1),ierr);CHKERRA(ierr)
      call MatSetOption(aux,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr);CHKERRA(ierr)
      call MatShift(aux,norm(1) * 1.0e-8,ierr);CHKERRA(ierr)
      call PCHPDDMSetAuxiliaryMat(pc,cols(1),aux,PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER,ierr);CHKERRA(ierr)
      call PCHPDDMHasNeumannMat(pc,PETSC_TRUE,ierr);CHKERRA(ierr)
      if (.not. normal) then
        call PCHPDDMGetCoarseCorrectionType(pc,ctype,ierr);CHKERRA(ierr)
        if (ctype == PC_HPDDM_COARSE_CORRECTION_DEFLATED) then
          call PCHPDDMSetCoarseCorrectionType(pc,PC_HPDDM_COARSE_CORRECTION_BALANCED,ierr);CHKERRA(ierr)
        endif
      endif
      call ISDestroy(cols(1),ierr);CHKERRA(ierr)
      call MatDestroy(aux,ierr);CHKERRA(ierr)
    endif
    ! end setup phase
    call MatDestroy(Bf,ierr);CHKERRA(ierr)
    call MatDestroy(Cf,ierr);CHKERRA(ierr)
  endif
  ! begin solution phase (section 4.1.3)
  call MatCreateVecs(A,x,b,ierr);CHKERRA(ierr)
  call VecSetRandom(b,PETSC_NULL_RANDOM,ierr);CHKERRA(ierr)
  call VecDuplicate(x,c,ierr);CHKERRA(ierr)
  call MatMultTranspose(A,b,c,ierr);CHKERRA(ierr)
  call KSPSolve(ksp,merge(c,b,normal),x,ierr);CHKERRA(ierr)
  ! end solution phase
  call VecScale(b,-1.0D0,ierr);CHKERRA(ierr)
  call MatMultAdd(A,x,b,b,ierr);CHKERRA(ierr)
  call VecNorm(b,NORM_2,norm(1),ierr);CHKERRA(ierr)
  call MatMultTranspose(A,b,c,ierr);CHKERRA(ierr)
  call VecNorm(c,NORM_2,norm(2),ierr);CHKERRA(ierr)
  if (rank .eq. 0) then
    call PetscPrintf(PETSC_COMM_SELF,"||A^T(Ax-b)|| / ||Ax-b|| = ",ierr);CHKERRA(ierr)
    write(name,*) norm(2)
    call PetscPrintf(PETSC_COMM_SELF,trim(name)//" / ",ierr);CHKERRA(ierr)
    write(name,*) norm(1)
    call PetscPrintf(PETSC_COMM_SELF,trim(name)//" = ",ierr);CHKERRA(ierr)
    norm(1) = norm(2) / norm(1)
    write(name,*) norm(1)
    call PetscPrintf(PETSC_COMM_SELF,trim(name)//"\n",ierr);CHKERRA(ierr)
  endif
  call KSPDestroy(ksp,ierr);CHKERRA(ierr)
  call VecDestroy(c,ierr);CHKERRA(ierr)
  call VecDestroy(x,ierr);CHKERRA(ierr)
  call VecDestroy(b,ierr);CHKERRA(ierr)
  call MatDestroy(A,ierr);CHKERRA(ierr)
  call PetscFinalize(ierr)
  end
