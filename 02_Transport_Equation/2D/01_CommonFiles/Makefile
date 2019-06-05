CC=gcc
#CC=tau_cc.sh
IDIR=
CFLAGS= -DTranspEquation2D -Ddebug -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DTransEquation2D -Ddebug -lm -lrt -Wall -Ofast -march=native
LDFLAGS=
SOURCES= main.c \
Preprocess.c \
Process.c \
Postprocess.c \
Fill_LM.c \
ede_Initialization.c \
csr_Initialization.c \
Build_M_K_R_SUPG.c \
Build_M_K_R_DD.c \
Build_M_F_DD.c \
eval_U_dU.c \
ebe_assembly.c \
ede_assembly.c \
csr_assembly.c \
Paraview_Output.c \
CAU_ShockCapture.c \
CAU_DD_ShockCapture.c \
YZBeta_ShockCapture.c \
h_shock_2sqrtArea.c \
h_shock_Option1.c \
h_shock_Option2.c \
setStabilizationForm.c \
setSolver.c \
setPreconditioner.c \
setScaling.c \
setProblem.c \
setMatrixVectorProductType.c \
setzeros.c \
calculate_DaB.c \
uB_InitialSolution.c \
../CARTOLA/CARTOLA_Condutivity.c \
../CARTOLA/CARTOLA_Font.c \
../CARTOLA/CARTOLA_Reaction.c \
../CARTOLA/CARTOLA_Velocity.c \
../CARTOLA/CARTOLA_upresc.c \
../CARTOLA/CARTOLA_InitialSolution.c \
../CONE/CONE_Condutivity.c \
../CONE/CONE_Font.c \
../CONE/CONE_Reaction.c \
../CONE/CONE_Velocity.c \
../CONE/CONE_upresc.c \
../CONE/CONE_InitialSolution.c \
../CONE2/CONE2_Condutivity.c \
../CONE2/CONE2_Font.c \
../CONE2/CONE2_Reaction.c \
../CONE2/CONE2_Velocity.c \
../CONE2/CONE2_upresc.c \
../CONE2/CONE2_InitialSolution.c \
../PUDIM/PUDIM_Condutivity.c \
../PUDIM/PUDIM_Font.c \
../PUDIM/PUDIM_Reaction.c \
../PUDIM/PUDIM_Velocity.c \
../PUDIM/PUDIM_upresc.c \
../PUDIM/PUDIM_InitialSolution.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/C_Operations/dmemcpy.c \
../../../00_CommonFiles/C_Operations/memsetzero.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemv.c \
../../../00_CommonFiles/MatrixVector_Operations/EDE/edemv.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pcg.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/Diag_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/Diag_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/Diag_precond_CSR_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGS_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SGS_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSOR_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SSOR_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/LU_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/LU_precond_EDE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SGS_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SGS_precond_EDE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/SSOR_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/SSOR_precond_EDE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EBE/LU_precond_EBE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/EDE/LU_precond_EDE.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/scaling.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c \
../../../00_CommonFiles/Time_Integration/Criteria.c \
../../../00_CommonFiles/Time_Integration/Predictor_New.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old.c \
../../../00_CommonFiles/Time_Integration/Predictor_New_BDF.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old_BDF.c \
../../../00_CommonFiles/Time_Integration/Predictor_New_TRBDF2.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old_TRBDF2.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=TranspEquation2D

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DTranspEquation2D -Ddebug -c -lm -lrt -Wall -g
debug: GCFLAGS= -DTransEquation2D -Ddebug -lm -lrt -Wall -g
debug: all

gprof: CFLAGS+=-pg -Ddebug
gprof: LDFLAGS+=-pg -Ddebug
gprof: all

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
