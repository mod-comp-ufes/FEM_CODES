CC=gcc
IDIR= 
CFLAGS= -DSSNavierStokesEquations2D -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DSSNavierStokesEquations2D -std=gnu99 -lm -lrt -Wall \
-lrt -lm -lgfortran -Wall -lblas -llapack -lquadmath -Ofast -march=native
LDFLAGS= 
SOURCES= main.c \
Preprocess.c \
Process.c \
Fill_ID.c \
Fill_LM.c \
setProblem.c \
setSolver.c \
setStabilizationForm.c \
setPreconditioner.c \
setMatrixVectorProductType.c \
csr_Initialization.c \
ebe_assembly.c \
csr_assembly.c \
Build_K_F_SUPG.c \
Build_K_F_MS.c \
Build_K_F_MS_Time.c \
eta_newton.c\
Coef_Lift_Drag.c\
eval_U.c \
eval_U_dU.c \
Postprocess.c \
Paraview_Output.c \
setzeros.c\
../CAVITY/CAVITY_v1presc.c \
../CAVITY/CAVITY_v2presc.c \
../CAVITY/CAVITY_ppresc.c \
../CAVITY/CAVITY_f1ext.c \
../CAVITY/CAVITY_f2ext.c \
../CHANNEL/CHANNEL_v1presc.c \
../CHANNEL/CHANNEL_v2presc.c \
../CHANNEL/CHANNEL_ppresc.c \
../CHANNEL/CHANNEL_f1ext.c \
../CHANNEL/CHANNEL_f2ext.c \
../EXATA/EXATA_v1presc.c \
../EXATA/EXATA_v2presc.c \
../EXATA/EXATA_ppresc.c \
../EXATA/EXATA_f1ext.c \
../EXATA/EXATA_f2ext.c \
../CYLINDER/CYLINDER_v1presc.c \
../CYLINDER/CYLINDER_v2presc.c \
../CYLINDER/CYLINDER_ppresc.c \
../CYLINDER/CYLINDER_f1ext.c \
../CYLINDER/CYLINDER_f2ext.c \
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/BLAS_Operations/izero.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/amg_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_dpa.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_dpa_counting.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_fib_heap.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_list.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_matrix.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_solvers.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_util.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/bcm_linmatch.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/BlockDiagDOF3_precond_EBE_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/BlockDiagDOF3_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemvNDOF3.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c 

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=SSNavierStokesEquations2D 

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DSSNavierStokesEquations2D -g -Ddebug -c -lm -lrt -Wall
debug: GCFLAGS= -DSSNavierStokesEquations2D -g -std=gnu99 -Ddebug -lm -lrt -Wall \
-lrt -lm -lgfortran -Wall -lblas -llapack -lquadmath 
debug: all

fast: CFLAGS= -DSSNavierStokesEquations2D -c -std=gnu99 -lm -lrt -Wall -Ofast -march=native
fast: GCFLAGS= -DSSNavierStokesEquations2D -lm -lrt -Wall \
-lrt -lm -lgfortran -Wall -lblas -llapack -lquadmath -Ofast -march=native
fast: all

gprof: CFLAGS+=-pg  
gprof: LDFLAGS+=-pg
gprof: all
 
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)



