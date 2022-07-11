/*
   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2021-07-11

   Copyright (C) 2021-     Centre National de la Recherche Scientifique

   This script is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

   If you use this script, you are kindly asked to cite the following article:

   "A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations",
   H. Al Daas, P. Jolivet, and J. A. Scott (2022).
 */

#include <petsc.h>

static char help[] = "Solves a linear least squares problem system using PCHPDDM.\n\n";

int main(int argc,char **args)
{
  Vec         x,b,c;
  Mat         A,B,*Neumann,aux,C = NULL;
  KSP         ksp;
  PC          pc;
  IS          is,rows,cols;
  PetscMPIInt rank;
  PetscReal   norm[2];
  PetscInt    m,n;
  PetscBool   normal = PETSC_FALSE,flg,qr = PETSC_FALSE;
  PetscViewer viewer;
  char        name[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc,&args,NULL,help));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
  // begin loading phase (section 4.1.1)
  PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
  PetscCall(PetscOptionsGetString(NULL,NULL,"-mat_name",name,sizeof(name),&flg));
  PetscCheck(flg,PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify matrix name with -mat_name");
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_READ,&viewer));
  PetscCall(MatLoad(A,viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  // end loading phase
  PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PetscObjectTypeCompare((PetscObject)pc,PCQR,&flg));
  if (flg) {
    PetscCall(KSPSetOperators(ksp,A,A));
    PetscCall(KSPSetType(ksp,KSPPREONLY));
  } else {
    MatPartitioning mpart;
    Mat             perm;
    // begin partitioning phase (section 4.1.1)
    PetscCall(MatProductCreate(A,A,NULL,&B));
    PetscCall(MatProductSetType(B,MATPRODUCT_AtB));
    PetscCall(MatProductSetFromOptions(B));
    PetscCall(MatProductSymbolic(B));
    PetscCall(MatPartitioningCreate(PETSC_COMM_WORLD,&mpart));
    PetscCall(MatPartitioningSetAdjacency(mpart,B));
    PetscCall(MatPartitioningSetFromOptions(mpart));
    PetscCall(MatPartitioningApply(mpart,&is));
    PetscCall(MatPartitioningDestroy(&mpart));
    PetscCall(ISBuildTwoSided(is,NULL,&rows));
    PetscCall(ISDestroy(&is));
    PetscCall(PetscOptionsGetBool(NULL,NULL,"-pc_use_qr",&qr,NULL));
    if (!qr) PetscCall(MatProductNumeric(B));
    PetscCall(MatCreateSubMatrix(B,rows,rows,MAT_INITIAL_MATRIX,&perm));
    PetscCall(MatHeaderReplace(B,&perm));
    if (!qr) {
      PetscCall(MatNorm(B,NORM_FROBENIUS,norm));
      PetscCall(MatSetOption(B,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE));
      PetscCall(MatSetOption(B,MAT_SYMMETRIC,PETSC_TRUE));
      PetscCall(MatShift(B,*norm * 1.0e-10));
    }
    PetscCall(MatGetOwnershipRange(A,&m,&n));
    PetscCall(ISCreateStride(PETSC_COMM_WORLD,n-m,m,1,&is));
    PetscCall(MatCreateSubMatrix(A,is,rows,MAT_INITIAL_MATRIX,&perm));
    PetscCall(MatHeaderReplace(A,&perm));
    PetscCall(ISDestroy(&is));
    PetscCall(ISDestroy(&rows));
    // end partitioning phase
    // begin setup phase (section 4.1.2)
    PetscCall(PetscObjectTypeCompare((PetscObject)ksp,KSPLSQR,&normal));
    normal = (PetscBool)!normal;
    if (normal || qr) PetscCall(MatCreateNormal(A,&C));
    PetscCall(KSPSetOperators(ksp,normal?C:A,qr?C:B));
    PetscCall(PetscObjectTypeCompare((PetscObject)pc,PCHPDDM,&flg));
    if (flg) {
      flg = PETSC_FALSE;
      PetscCall(PetscOptionsGetBool(NULL,NULL,"-pc_hidden_setup",&flg,NULL));
      if (!flg) {
        PetscCall(MatGetOwnershipRangeColumn(A,&m,&n));
        PetscCall(ISCreateStride(PETSC_COMM_SELF,n-m,m,1,&cols));
        PetscCall(MatGetSize(A,&m,NULL));
        PetscCall(ISCreateStride(PETSC_COMM_SELF,m,0,1,&rows));
        PetscCall(MatSetOption(A,MAT_SUBMAT_SINGLEIS,PETSC_TRUE));
        PetscCall(MatCreateSubMatrices(A,1,&rows,&cols,MAT_INITIAL_MATRIX,&Neumann));
        PetscCall(MatFindZeroRows(*Neumann,&is));
        PetscCall(MatDestroySubMatrices(1,&Neumann));
        PetscCall(MatIncreaseOverlap(B,1,&cols,1));
        PetscCall(MatCreateSubMatrices(A,1,&rows,&cols,MAT_INITIAL_MATRIX,&Neumann));
        PetscCall(ISDestroy(&rows));
        PetscCall(MatSetOption(*Neumann,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE));
        PetscCall(MatZeroRowsIS(*Neumann,is,0.0,NULL,NULL));
        PetscCall(ISDestroy(&is));
        if (!qr) {
          PetscCall(MatTransposeMatMult(*Neumann,*Neumann,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&aux));
          PetscCall(MatNorm(aux,NORM_FROBENIUS,norm));
          PetscCall(MatSetOption(aux,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE));
          PetscCall(MatShift(aux,*norm * 1.0e-8));
        } else PetscCall(MatCreateNormal(Neumann[0],&aux));
        PetscCall(MatDestroySubMatrices(1,&Neumann));
        PetscCall(PCHPDDMSetAuxiliaryMat(pc,cols,aux,NULL,NULL));
        PetscCall(PCHPDDMHasNeumannMat(pc,PETSC_TRUE));
        if (!normal) {
          PCHPDDMCoarseCorrectionType type;
          PetscCall(PCHPDDMGetCoarseCorrectionType(pc,&type));
          if (type == PC_HPDDM_COARSE_CORRECTION_DEFLATED) PetscCall(PCHPDDMSetCoarseCorrectionType(pc,PC_HPDDM_COARSE_CORRECTION_BALANCED));
        }
        PetscCall(ISDestroy(&cols));
        PetscCall(MatDestroy(&aux));
      }
    }
    // end setup phase
    PetscCall(MatDestroy(&B));
    PetscCall(MatDestroy(&C));
  }
  // begin solution phase (section 4.1.3)
  PetscCall(MatCreateVecs(A,&x,&b));
  PetscCall(VecSetRandom(b,NULL));
  PetscCall(VecDuplicate(x,&c));
  PetscCall(MatMultTranspose(A,b,c));
  PetscCall(KSPSolve(ksp,normal?c:b,x));
  // end solution phase
  PetscCall(VecScale(b,-1.0));
  PetscCall(MatMultAdd(A,x,b,b));
  PetscCall(VecNorm(b,NORM_2,norm));
  PetscCall(MatMultTranspose(A,b,c));
  PetscCall(VecNorm(c,NORM_2,norm+1));
  if (!rank) PetscCall(PetscPrintf(PETSC_COMM_SELF,"||A^T(Ax-b)|| / ||Ax-b|| = %f / %f = %f\n",(double)norm[1],(double)norm[0],(double)(norm[1]/norm[0])));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&c));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
