-include ${PETSC_DIR}/petscdir.mk
EXAMPLESC = sparse_ls.c
EXAMPLESF = sparse_ls_f.F90
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

all: sparse_ls sparse_ls_f

clean::
	rm -rf sparse_ls sparse_ls_f output/*tmp.out

check: all sparse_ls.py
	@for BinName in "./sparse_ls" "./sparse_ls_f" "${PYTHON} ./sparse_ls.py"; do \
		for PCType in hpddm gamg asm hypre qr; do \
			for MatName in lp_stocfor3 mesh_deform; do \
				if [ "$${PCType}" = "qr" ]; then NP=1; else NP=4; fi; \
				cmd="${MPIEXEC} -np $${NP} $${BinName} -mat_name datafiles/$${MatName}.dat -pc_type $${PCType} -options_file default.rc -ksp_view"; \
				echo "$${cmd}"; \
				$${cmd} > output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out || exit; \
				${PETSC_DIR}/lib/petsc/bin/petscdiff output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.out output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out; \
				unlink output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out 2> /dev/null; \
			done \
		done \
	done; \
	for PCType in hpddm gamg asm hypre; do \
		for MatName in lp_stocfor3 mesh_deform; do \
			cmd="${MPIEXEC} -np 4 ./sparse_ls -mat_name datafiles/$${MatName}.dat -pc_type $${PCType} -options_file gmres.rc -ksp_view"; \
			echo "$${cmd}"; \
			$${cmd} > output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out || exit; \
			${PETSC_DIR}/lib/petsc/bin/petscdiff output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.out output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out; \
			unlink output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out 2> /dev/null; \
		done \
	done
