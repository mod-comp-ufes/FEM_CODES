CC=gcc
#CC=tau_cc.sh
IDIR=
CFLAGS= -DSSTransportEquation3D -g -c -lm -lrt -Wall -Ofast -march=native
GCFLAGS= -DSSTransportEquation3D -g -lm -lrt -Wall -lrt -lm -Wall -Ofast -march=native
LDFLAGS=
SOURCES= main.c \
Preprocess.c \
Process.c \
Postprocess.c \
Fill_ID.c \
Fill_LM.c \
setProblem.c \
setMatrixVectorProductType.c \
setSolver.c \
setPreconditioner.c \
setStabilizationForm.c \
setZeros.c \
AssemblyGlobalMatrix.c \
Element_Configuration.c \
ElementsType.c \
LengthMesh.c \
eval_U_Space.c \
ebe_assembly.c \
csr_Initialization.c \
csr_assembly.c \
font_assembly.c \
Space_Algorithm_NonLinear.c \
Space_Algorithm_Damping.c \
Space_Algorithm_Galerkin.c \
Build.c \
Build_SUPG.c \
Build_CAU.c \
Build_VMS.c \
Build_NSGS.c \
Build_DD.c \
AproximationInEachNode.c \
Paraview_Output.c \
Paraview_Exact.c \
Gnuplot_Output.c \
Norm_L2.c \
Norm_H1.c \
Calculating_Errors.c \
Compute_Residual.c \
TauMacro.c \
TauMicro.c \
../APLI4/APLI4_upresc.c \
../APLI4/APLI4_f.c \
../APLI4/APLI4_ExactSolution.c \
../APLI4/APLI4_kappa.c \
../APLI4/APLI4_beta.c \
../APLI4/APLI4_sigma.c \
../APLI6/APLI6_upresc.c \
../APLI6/APLI6_f.c \
../APLI6/APLI6_ExactSolution.c \
../APLI6/APLI6_kappa.c \
../APLI6/APLI6_beta.c \
../APLI6/APLI6_sigma.c \
../DIFCONV/DIFCONV_upresc.c \
../DIFCONV/DIFCONV_f.c \
../DIFCONV/DIFCONV_ExactSolution.c \
../DIFCONV/DIFCONV_kappa.c \
../DIFCONV/DIFCONV_beta.c \
../DIFCONV/DIFCONV_sigma.c \
../DIFCONV/DIFCONV_DuDx.c \
../DIFCONV/DIFCONV_DuDy.c \
../DIFCONV/DIFCONV_DuDz.c \
../PAREDE/PAREDE_upresc.c \
../PAREDE/PAREDE_f.c \
../PAREDE/PAREDE_kappa.c \
../PAREDE/PAREDE_beta.c \
../PAREDE/PAREDE_sigma.c \
../RAMPA2/RAMPA2_upresc.c \
../RAMPA2/RAMPA2_f.c \
../RAMPA2/RAMPA2_kappa.c \
../RAMPA2/RAMPA2_beta.c \
../RAMPA2/RAMPA2_sigma.c \
../DIFCONVREA/DIFCONVREA_upresc.c \
../DIFCONVREA/DIFCONVREA_f.c \
../DIFCONVREA/DIFCONVREA_ExactSolution.c \
../DIFCONVREA/DIFCONVREA_kappa.c \
../DIFCONVREA/DIFCONVREA_beta.c \
../DIFCONVREA/DIFCONVREA_sigma.c \
../DIFCONVREA/DIFCONVREA_DuDx.c \
../DIFCONVREA/DIFCONVREA_DuDy.c \
../DIFCONVREA/DIFCONVREA_DuDz.c \
../HIGH/HIGH_upresc.c \
../HIGH/HIGH_f.c \
../HIGH/HIGH_kappa.c \
../HIGH/HIGH_beta.c \
../HIGH/HIGH_sigma.c \
../HEMKER/HEMKER_upresc.c \
../HEMKER/HEMKER_f.c \
../HEMKER/HEMKER_kappa.c \
../HEMKER/HEMKER_beta.c \
../HEMKER/HEMKER_sigma.c \
../DIFREA/DIFREA_upresc.c \
../DIFREA/DIFREA_f.c \
../DIFREA/DIFREA_ExactSolution.c \
../DIFREA/DIFREA_kappa.c \
../DIFREA/DIFREA_beta.c \
../DIFREA/DIFREA_sigma.c \
../DIFREA/DIFREA_DuDx.c \
../DIFREA/DIFREA_DuDy.c \
../DIFREA/DIFREA_DuDz.c \
../EXDIFREADOM/EXDIFREADOM_upresc.c \
../EXDIFREADOM/EXDIFREADOM_f.c \
../EXDIFREADOM/EXDIFREADOM_kappa.c \
../EXDIFREADOM/EXDIFREADOM_beta.c \
../EXDIFREADOM/EXDIFREADOM_sigma.c \
../../../00_CommonFiles/Allocation_Operations/mycalloc.c \
../../../00_CommonFiles/CPU_Time_Operations/calculateTime.c \
../../../00_CommonFiles/IO_Operations/myfopen.c \
../../../00_CommonFiles/MatrixVector_Operations/EBE/ebemv3Dndof1.c \
../../../00_CommonFiles/MatrixVector_Operations/CSR/csrmv.c \
../../../00_CommonFiles/BLAS_Operations/ddot.c \
../../../00_CommonFiles/BLAS_Operations/dzero.c \
../../../00_CommonFiles/BLAS_Operations/dcopy.c \
../../../00_CommonFiles/BLAS_Operations/daxpy.c \
../../../00_CommonFiles/BLAS_Operations/dscal.c \
../../../00_CommonFiles/BLAS_Operations/determinant.c \
../../../00_CommonFiles/BLAS_Operations/minor_matrix.c \
../../../00_CommonFiles/BLAS_Operations/inverse_matrix.c \
../../../00_CommonFiles/BLAS_Operations/matrix_multiplication.c \
../../../00_CommonFiles/BLAS_Operations/dmax.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/pgmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/gmres.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/NO_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/LU_precond_EBE_setup_NNOEL4.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/LU_precond_EBE_NNOEL4.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond_setup.c \
../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/ILUp_precond.c \
../../../00_CommonFiles/Reordering/reordering.c \
../../../00_CommonFiles/Reordering/symrcm.c \

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=SSTransportEquation3D

all: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS= -DSSTransportEquation3D -Ddebug -c -lm -lrt -Wall -O0 -g
debug: GCFLAGS= -DSSTransportEquation3D -Ddebug -lm -lrt -Wall -g
debug: all

fast: CFLAGS= -DSSTransportEquation3D -c -lm -lrt -Wall -Ofast -march=native
fast: GCFLAGS= -DSSTransportEquation3D -lm -lrt -Wall -Ofast -march=native
fast: all

gprof: CFLAGS+=-pg -Ddebug
gprof: LDFLAGS+=-pg -Ddebug
gprof: all

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(GCFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
