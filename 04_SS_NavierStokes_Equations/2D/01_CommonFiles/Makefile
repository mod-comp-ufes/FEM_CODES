CC=gcc
IDIR= 
CFLAGS= -DSSNavierStokesEquations2D -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DSSNavierStokesEquations2D -lm -lrt -Wall -Ofast -march=native
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
setzeros.c \
setMatrixVectorProductType.c \
csr_Initialization.c \
ebe_assembly.c \
csr_assembly.c \
Build_K_F_SUPG_PSPG.c \
Build_K_F_VMS_DCDD.c \
eval_U.c \
Postprocess.c \
Paraview_Output.c \
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
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/BLAS_Operations/izero.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
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
debug: GCFLAGS= -DSSNavierStokesEquations2D -g -Ddebug -lm -lrt -Wall
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



