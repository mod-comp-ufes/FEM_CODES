CC=gcc
#'CC=tau_cc.sh
IDIR=
CFLAGS= -DEulerEquations2D -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DEulerEquations2D -lm -lrt -Wall -Ofast -march=native
LDFLAGS=
SOURCES= main.c \
Preprocess.c \
Process.c \
Postprocess.c \
printf_rho_Residue.c \
Fill_ID.c \
Fill_LM.c \
setProblem.c \
setMatrixVectorProductType.c \
setSolver.c \
setStabilizationForm.c \
setzeros.c \
setPreconditioner.c \
setScaling.c \
setDimensionlessness.c \
ebe_assembly.c \
ede_Initialization.c \
ede_assembly.c \
csr_Initialization.c \
csr_assembly.c \
Build_M_K_F_SUPG.c \
Build_M_K_F_DD_Transiente.c \
Build_M_F_DD_Transiente.c \
eval_U_dU.c \
Ax_Ay_calculations.c \
no_penetrability.c \
NO_BC_no_penetrability.c \
BC_General_theta.c \
set_BC_no_penetrability.c \
Delta_CAU.c \
Delta_DD.c \
Delta_YZBeta.c \
Delta_YZBetaNMV.c \
Paraview_Output.c \
Paraview_Output_3D.c \
Paraview_Output_DeltaT.c \
Paraview_Output_3D_DeltaT.c \
calculate_DaB.c \
uB_InitialSolution.c \
rotation.c \
../EXPLOSION/EXPLOSION_InitialSolution.c \
../EXPLOSION/EXPLOSION_gamma.c \
../EXPLOSION/EXPLOSION_cv.c \
../EXPLOSION/EXPLOSION_rhopresc.c \
../EXPLOSION/EXPLOSION_v1presc.c \
../EXPLOSION/EXPLOSION_v2presc.c \
../EXPLOSION/EXPLOSION_epresc.c \
../OBLIQUO/OBLIQUO_InitialSolution.c \
../OBLIQUO/OBLIQUO_gamma.c \
../OBLIQUO/OBLIQUO_cv.c \
../OBLIQUO/OBLIQUO_rhopresc.c \
../OBLIQUO/OBLIQUO_v1presc.c \
../OBLIQUO/OBLIQUO_v2presc.c \
../OBLIQUO/OBLIQUO_epresc.c \
../REFLETIDO/REFLETIDO_InitialSolution.c \
../REFLETIDO/REFLETIDO_gamma.c \
../REFLETIDO/REFLETIDO_cv.c \
../REFLETIDO/REFLETIDO_rhopresc.c \
../REFLETIDO/REFLETIDO_v1presc.c \
../REFLETIDO/REFLETIDO_v2presc.c \
../REFLETIDO/REFLETIDO_epresc.c \
../TUNEL/TUNEL_InitialSolution.c \
../TUNEL/TUNEL_gamma.c \
../TUNEL/TUNEL_cv.c \
../TUNEL/TUNEL_rhopresc.c \
../TUNEL/TUNEL_v1presc.c \
../TUNEL/TUNEL_v2presc.c \
../TUNEL/TUNEL_epresc.c \
../SOD/SOD_InitialSolution.c \
../SOD/SOD_gamma.c \
../SOD/SOD_cv.c \
../SOD/SOD_rhopresc.c \
../SOD/SOD_v1presc.c \
../SOD/SOD_v2presc.c \
../SOD/SOD_epresc.c \
../NACA0012/NACA0012_InitialSolution.c \
../NACA0012/NACA0012_gamma.c \
../NACA0012/NACA0012_cv.c \
../NACA0012/NACA0012_rhopresc.c \
../NACA0012/NACA0012_v1presc.c \
../NACA0012/NACA0012_v2presc.c \
../NACA0012/NACA0012_epresc.c \
../NACA0012/NACA0012_BC_no_penetrability.c \
../NACA0012/NACA0012_theta.c \
../CYLINDER/CYLINDER_InitialSolution.c \
../CYLINDER/CYLINDER_gamma.c \
../CYLINDER/CYLINDER_cv.c \
../CYLINDER/CYLINDER_rhopresc.c \
../CYLINDER/CYLINDER_v1presc.c \
../CYLINDER/CYLINDER_v2presc.c \
../CYLINDER/CYLINDER_epresc.c \
../CYLINDER/CYLINDER_BC_no_penetrability.c \
../CYLINDER/CYLINDER_theta.c \
../BAROCLINIC/BAROCLINIC_InitialSolution.c \
../BAROCLINIC/BAROCLINIC_gamma.c \
../BAROCLINIC/BAROCLINIC_cv.c \
../BAROCLINIC/BAROCLINIC_rhopresc.c \
../BAROCLINIC/BAROCLINIC_v1presc.c \
../BAROCLINIC/BAROCLINIC_v2presc.c \
../BAROCLINIC/BAROCLINIC_epresc.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/amg_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_dpa.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_dpa_counting.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_fib_heap.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_list.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_matrix.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_solvers.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_util.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/bcm_linmatch.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/Diag_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/Diag_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/Diag_precond_CSR_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/BlockDiag_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/NOBlockDiag2_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/NOBlockDiag2_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/BlockDiag_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/JacobiDOF4_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/JacobiBlockDOF4_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSORDOF4_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSORBlockDOF4_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/BlockDiag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/BlockDiag2_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NOBlockDiag2_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/JacobiDOF4_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/JacobiBlockDOF4_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSORDOF4_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSORBlockDOF4_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGSBlock_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGSBlock_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/LUBlock_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/LUBlock_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Block_scaling.c \
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/BLAS_Operations/izero.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemvNDOF4.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebe2mvNDOF4.c \
../../../00_CommonFiles/MatrixVector_Operations/EDE/edemvNDOF4.c \
../../../00_CommonFiles/MatrixVector_Operations/EDE/ede2mvNDOF4.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/C_Operations/dmemcpy.c \
../../../00_CommonFiles/C_Operations/memsetzero.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c \
../../../00_CommonFiles/Time_Integration/Criteria.c \
../../../00_CommonFiles/Time_Integration/Predictor_New.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old.c \
../../../00_CommonFiles/Time_Integration/Predictor_New_BDF.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old_BDF.c \
../../../00_CommonFiles/Time_Integration/Predictor_New_TRBDF2.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old_TRBDF2.c \

#Build_M_K_F_DD_Estatico.c \
#Delta_SHOC1.c \
#Delta_SHOC2.c \
#Delta_SHOC3.c \
#Delta_SHOC4.c \
#ebe_assembly2.c \
#Reverse_Matrix.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/gmresBD.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=EulerEquations2D

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DEulerEquations2D -g -Ddebug -Dcheck_memory_leak -c -lm -lrt -Wall
debug: GCFLAGS= -DEulerEquations2D -g -Ddebug -Dcheck_memory_leak -lm -lrt -Wall 
#../../../00_CommonFiles/Reordering/hsl_mc73/hsl_mc73d.o 
debug: all

gprof: CFLAGS+=-pg
gprof: LDFLAGS+=-pg
gprof: all


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
