/*
   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2021-07-11

   Copyright (C) 2021-     Centre National de la Recherche Scientifique

   This script is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

   If you use this script, you are kindly asked to cite the following article:

   "A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations",
   H. Al Daas, P. Jolivet, and J. A. Scott (2021).
 */

#include <petsc.h>

static char help[] = "Solves a linear least squares problem system using PCHPDDM.\n\n";

int main(int argc,char **args)
{
  Vec            x,b,c;
  Mat            A,B,*Neumann,aux,C = NULL;
  KSP            ksp;
  PC             pc;
  IS             is,rows,cols;
  PetscMPIInt    rank;
  PetscReal      norm[2];
  PetscInt       m,n;
  PetscBool      normal = PETSC_FALSE,flg;
  PetscViewer    viewer;
  char           name[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc,&args,NULL,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRMPI(ierr);
  // begin loading phase (section 4.1.1)
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(NULL,NULL,"-mat_name",name,sizeof(name),&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify matrix name with -mat_name");
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = MatLoad(A,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  // end loading phase
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)pc,PCQR,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  } else {
    MatPartitioning mpart;
    Mat             perm;
    // begin partitioning phase (section 4.1.1)
    ierr = MatTransposeMatMult(A,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&B);CHKERRQ(ierr);
    ierr = MatPartitioningCreate(PETSC_COMM_WORLD,&mpart);CHKERRQ(ierr);
    ierr = MatPartitioningSetAdjacency(mpart,B);CHKERRQ(ierr);
    ierr = MatPartitioningSetFromOptions(mpart);CHKERRQ(ierr);
    ierr = MatPartitioningApply(mpart,&is);CHKERRQ(ierr);
    ierr = MatPartitioningDestroy(&mpart);CHKERRQ(ierr);
    ierr = ISBuildTwoSided(is,NULL,&rows);CHKERRQ(ierr);
    ierr = ISDestroy(&is);CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(B,rows,rows,MAT_INITIAL_MATRIX,&perm);CHKERRQ(ierr);
    ierr = MatDestroy(&B);CHKERRQ(ierr);
    B = perm;
    ierr = MatGetOwnershipRange(A,&m,&n);CHKERRQ(ierr);
    ierr = ISCreateStride(PETSC_COMM_WORLD,n-m,m,1,&is);CHKERRQ(ierr);
    ierr = MatCreateSubMatrix(A,is,rows,MAT_INITIAL_MATRIX,&perm);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = ISDestroy(&is);CHKERRQ(ierr);
    ierr = ISDestroy(&rows);CHKERRQ(ierr);
    A = perm;
    // end partitioning phase
    // begin setup phase (section 4.1.2)
    ierr = MatNorm(B,NORM_FROBENIUS,norm);CHKERRQ(ierr);
    ierr = MatSetOption(B,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
    ierr = MatSetOption(B,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatShift(B,*norm * 1.0e-10);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)ksp,KSPLSQR,&normal);CHKERRQ(ierr);
    normal = (PetscBool)!normal;
    if (normal) {
      ierr = MatCreateNormal(A,&C);CHKERRQ(ierr);
    }
    ierr = KSPSetOperators(ksp,normal?C:A,B);CHKERRQ(ierr);
    ierr = PetscObjectTypeCompare((PetscObject)pc,PCHPDDM,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = MatGetOwnershipRange(B,&m,&n);CHKERRQ(ierr);
      ierr = ISCreateStride(PETSC_COMM_SELF,n-m,m,1,&cols);CHKERRQ(ierr);
      ierr = MatGetSize(A,&m,NULL);CHKERRQ(ierr);
      ierr = ISCreateStride(PETSC_COMM_SELF,m,0,1,&rows);CHKERRQ(ierr);
      ierr = MatSetOption(A,MAT_SUBMAT_SINGLEIS,PETSC_TRUE);CHKERRQ(ierr);
      ierr = MatCreateSubMatrices(A,1,&rows,&cols,MAT_INITIAL_MATRIX,&Neumann);CHKERRQ(ierr);
      ierr = MatFindZeroRows(*Neumann,&is);CHKERRQ(ierr);
      ierr = MatDestroySubMatrices(1,&Neumann);CHKERRQ(ierr);
      ierr = MatIncreaseOverlap(B,1,&cols,1);CHKERRQ(ierr);
      ierr = MatCreateSubMatrices(A,1,&rows,&cols,MAT_INITIAL_MATRIX,&Neumann);CHKERRQ(ierr);
      ierr = ISDestroy(&rows);CHKERRQ(ierr);
      ierr = MatZeroRowsIS(*Neumann,is,0.0,NULL,NULL);CHKERRQ(ierr);
      ierr = ISDestroy(&is);CHKERRQ(ierr);
      ierr = MatTransposeMatMult(*Neumann,*Neumann,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&aux);CHKERRQ(ierr);
      ierr = MatDestroySubMatrices(1,&Neumann);CHKERRQ(ierr);
      ierr = MatNorm(aux,NORM_FROBENIUS,norm);CHKERRQ(ierr);
      ierr = MatSetOption(aux,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);
      ierr = MatShift(aux,*norm * 1.0e-8);CHKERRQ(ierr);
      ierr = PCHPDDMSetAuxiliaryMat(pc,cols,aux,NULL,NULL);CHKERRQ(ierr);
      ierr = PCHPDDMHasNeumannMat(pc,PETSC_TRUE);CHKERRQ(ierr);
      ierr = ISDestroy(&cols);CHKERRQ(ierr);
      ierr = MatDestroy(&aux);CHKERRQ(ierr);
    }
    // end setup phase
    ierr = MatDestroy(&B);CHKERRQ(ierr);
    ierr = MatDestroy(&C);CHKERRQ(ierr);
  }
  // begin solution phase (section 4.1.3)
  ierr = MatCreateVecs(A,&x,&b);CHKERRQ(ierr);
  ierr = VecSetRandom(b,NULL);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&c);CHKERRQ(ierr);
  ierr = MatMultTranspose(A,b,c);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,normal?c:b,x);CHKERRQ(ierr);
  // end solution phase
  ierr = VecScale(b,-1.0);CHKERRQ(ierr);
  ierr = MatMultAdd(A,x,b,b);CHKERRQ(ierr);
  ierr = VecNorm(b,NORM_2,norm);CHKERRQ(ierr);
  ierr = MatMultTranspose(A,b,c);CHKERRQ(ierr);
  ierr = VecNorm(c,NORM_2,norm+1);CHKERRQ(ierr);
  if (!rank) {
    ierr = PetscPrintf(PETSC_COMM_SELF,"||A^T(Ax-b)|| / ||Ax-b|| = %f / %f = %f\n",norm[1],norm[0],norm[1]/norm[0]);CHKERRQ(ierr);
  }
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&c);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}
