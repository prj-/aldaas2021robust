#
#  Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
#       Date: 2021-07-11
#
#  Copyright (C) 2021-     Centre National de la Recherche Scientifique
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  If you use this script, you are kindly asked to cite the following article:
#
#  "A Robust Algebraic Domain Decomposition Preconditioner for Sparse Normal Equations",
#  H. Al Daas, P. Jolivet, and J. A. Scott (2021).
#

import petsc4py,sys
petsc4py.init(sys.argv)
from petsc4py import PETSc

A = PETSc.Mat().create(PETSc.COMM_WORLD)
# begin loading phase (section 4.1.1)
viewer = PETSc.Viewer().createBinary(PETSc.Options().getString('-mat_name'))
A.load(viewer)
viewer.destroy()
# end loading phase
x = PETSc.Vec()
b = PETSc.Vec()
c = PETSc.Vec()
ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
ksp.setFromOptions()
pc = ksp.getPC()
normal = False
if pc.getType() == PETSc.PC.Type.QR:
  ksp.setOperators(A,A)
  ksp.setType(PETSc.KSP.Type.PREONLY)
else:
  # begin partitioning phase (section 4.1.1)
  B = A.transposeMatMult(A)
  mpart = PETSc.MatPartitioning().create(PETSc.COMM_WORLD)
  mpart.setFromOptions()
  mpart.setAdjacency(B)
  isp = PETSc.IS()
  mpart.apply(isp)
  rows = PETSc.IS()
  rows = isp.buildTwoSided()
  perm = B.createSubMatrix(rows,rows)
  B = perm
  (m,n) = A.getOwnershipRange()
  isp = PETSc.IS().createStride(n-m,m,1,PETSc.COMM_WORLD)
  perm = A.createSubMatrix(isp,rows)
  A = perm
  for obj in [mpart,isp,rows]:
    obj.destroy()
  # end partitioning phase
  # begin setup phase (section 4.1.2)
  B.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,False)
  B.setOption(PETSc.Mat.Option.SYMMETRIC,True)
  B.shift(B.norm() * 1.0e-10)
  C = PETSc.Mat()
  if not ksp.getType() == PETSc.KSP.Type.LSQR:
    normal = True
    C.createNormal(A)
  ksp.setOperators(C if normal else A,B)
  if pc.getType() == PETSc.PC.Type.HPDDM:
    (m,n) = B.getOwnershipRange()
    cols = PETSc.IS().createStride(n-m,m,1,PETSc.COMM_SELF)
    (m,n) = A.getSize()
    rows.createStride(m,0,1,PETSc.COMM_SELF)
    A.setOption(PETSc.Mat.Option.SUBMAT_SINGLEIS,True)
    Neumann = A.createSubMatrices(rows,cols)
    isp = Neumann[0].findZeroRows()
    B.increaseOverlap(cols)
    Neumann = A.createSubMatrices(rows,cols)
    Neumann[0].zeroRows(isp,0.0)
    aux = PETSc.Mat()
    aux = Neumann[0].transposeMatMult(Neumann[0])
    aux.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    aux.shift(aux.norm() * 1.0e-8)
    pc.setHPDDMAuxiliaryMat(cols,aux)
    pc.setHPDDMHasNeumannMat(True)
  # end setup phase
# begin solution phase (section 4.1.3)
(x,b) = A.createVecs()
b.setRandom()
c = x.duplicate()
A.multTranspose(b,c)
ksp.solve(c if normal else b,x)
# end solution phase (section 4.1.3)
b.scale(-1.0)
A.multAdd(x,b,b)
A.multTranspose(b,c)
norm = [b.norm(),c.norm()]
if not PETSc.COMM_WORLD.Get_rank():
  print("||A^T(Ax-b)|| / ||Ax-b|| = "+str(norm[1])+" / "+str(norm[0])+" = "+str(norm[1]/norm[0]))
