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
	@if [ "${PETSC4PY}" = "yes" ]; then PyBin="${PYTHON} ./sparse_ls.py"; fi; \
	for BinName in "./sparse_ls" "./sparse_ls_f" "$${PyBin}"; do \
		if [ -z "$${BinName}" ]; then continue; fi; \
		for PCType in hpddm gamg asm hypre qr; do \
			for MatName in lp_stocfor3 mesh_deform; do \
				if [ "$${PCType}" = "qr" ]; then NP=1; else NP=4; fi; \
				cmd="${MPIEXEC} -np $${NP} $${BinName} -mat_name datafiles/$${MatName}.dat -pc_type $${PCType} -options_file default.rc -ksp_view"; \
				echo "$${cmd}"; \
				$${cmd} > output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out || exit; \
				breakdown=`grep -e DIVERGED_ITS -e DIVERGED_BREAKDOWN output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.out | wc -l | xargs echo`; \
				if [ "$${breakdown}" != "2" ]; then ${PETSC_DIR}/lib/petsc/bin/petscdiff output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.out output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out; fi; \
				unlink output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out 2> /dev/null; \
			done \
		done \
	done; \
	for PCType in hpddm gamg asm hypre; do \
		for MatName in lp_stocfor3 mesh_deform; do \
			cmd="${MPIEXEC} -np 4 ./sparse_ls -mat_name datafiles/$${MatName}.dat -pc_type $${PCType} -options_file gmres.rc -ksp_view"; \
			echo "$${cmd}"; \
			$${cmd} > output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out || exit; \
			breakdown=`grep -e DIVERGED_ITS -e DIVERGED_BREAKDOWN output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.out | wc -l | xargs echo`; \
			if [ "$${breakdown}" != "2" ]; then  ${PETSC_DIR}/lib/petsc/bin/petscdiff output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.out output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out; fi; \
			unlink output/sparse_ls_ksp_type-gmres_mat_name-$${MatName}_pc_type-$${PCType}.tmp.out 2> /dev/null; \
		done \
	done; \
	for PCType in hpddm; do \
		for MatName in mesh_deform; do \
			for UseQR in true false; do \
				for HiddenSetup in true; do \
					if [ "$${UseQR}" = "true" ]; then suffix="_pc_use_qr"; options="$${UseQR} -pc_hpddm_levels_1_sub_pc_type qr"; else suffix=; options=false; fi; \
					cmd="${MPIEXEC} -np 4 ./sparse_ls -mat_name datafiles/$${MatName}.dat -pc_type $${PCType} -options_file default.rc -pc_use_qr $${options} -ksp_view -pc_hidden_setup $${HiddenSetup}"; \
					echo "$${cmd}"; \
					$${cmd} > output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}$${suffix}.tmp.out || exit; \
					${PETSC_DIR}/lib/petsc/bin/petscdiff output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}$${suffix}.out output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}$${suffix}.tmp.out; \
					unlink output/sparse_ls_ksp_type-lsqr_mat_name-$${MatName}_pc_type-$${PCType}$${suffix}.tmp.out 2> /dev/null; \
				done \
			done \
		done \
	done
