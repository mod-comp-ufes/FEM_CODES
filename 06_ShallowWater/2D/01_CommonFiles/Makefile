EXE = ShalowWater
CFLAGS = -DShalowWater -lm -g# -lrt -Wall -Ofast -march=native

SRC = $(wildcard *.c) \
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
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/Diag_precond_CSR_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/Diag_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/BLAS_Operations/dtrsvUP.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/BLAS_Operations/izero.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemvNDOF3.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/C_Operations/dmemcpy.c \
../../../00_CommonFiles/C_Operations/memsetzero.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c #\
../../../00_CommonFiles/Time_Integration/Criteria.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old_BDF.c \
../../../00_CommonFiles/Time_Integration/Predictor_Old_TRBDF2.c \

all: $(SRC)
	gcc -o $(EXE) $^ $(CFLAGS) -DonlyCompile


1DDamBreak: $(wildcard ../1D_dam_break/*.c) $(SRC) # ../1D_dam_break/ is the path do implemented functions
	gcc -o $(EXE) $^ $(CFLAGS)

transcritical_flow: $(wildcard ../transcritical_flow/*.c) $(SRC) # ../transcritical_flow/ is the path do implemented functions
	gcc -o $(EXE) $^ $(CFLAGS)

DamBreakBumps: $(wildcard ../DamBreakBumps/*.c) $(SRC) # ../DamBreakBumps/ is the path do implemented functions
	gcc -o $(EXE) $^ $(CFLAGS)
