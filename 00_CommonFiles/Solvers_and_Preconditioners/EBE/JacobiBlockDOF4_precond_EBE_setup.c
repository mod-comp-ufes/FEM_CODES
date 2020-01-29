#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int JacobiBlockDOF4_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int neq = Parameters->neq;
	int nel = Parameters->nel;
	double denominator;
	double **A = MatrixData->A;
	double aux0, aux1, aux2, aux3, aux4, aux5; //auxiliares

	if (tag==1){
		MatrixData->invDeaux = mycalloc("invDeaux of 'Jacobi_precond_EBE_setup'",NDOF*NNOEL*nel,sizeof(double));
		MatrixData->invDe = mycalloc("invDe of 'Jacobi_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++)
			MatrixData->invDe[I] = &(MatrixData->invDeaux[NDOF*NNOEL*I]);
	}

	for (I = 0; I < nel; I++){
		aux0 = A[I][26] * A[I][39] - A[I][27] * A[I][38];
    aux1 = A[I][27] * A[I][37] - A[I][25] * A[I][39];
    aux2 = A[I][25] * A[I][38] - A[I][26] * A[I][37];
    aux3 = A[I][24] * A[I][39] - A[I][27] * A[I][36];
    aux4 = A[I][26] * A[I][36] - A[I][24] * A[I][38];
    aux5 = A[I][24] * A[I][37] - A[I][25] * A[I][36];

    denominator = (A[I][0] * (A[I][13] * (aux0) + A[I][14] * (aux1) + A[I][15] * (aux2)) +
		               A[I][1] * (A[I][12] * (-aux0) + A[I][14] * (aux3) + A[I][15] * (aux4)) +
		               A[I][2] * (A[I][12] * (aux1) + A[I][13] * (-aux3) + A[I][15] * (aux5)) +
		               A[I][3] * (A[I][12] * (-aux2) + A[I][13] * (-aux4) + A[I][14] * (-aux5)));

    MatrixData->invDe[I][0] = 1.0 / denominator;

    aux0 = A[I][78] * A[I][91] - A[I][79] * A[I][90];
    aux1 = A[I][79] * A[I][89] - A[I][77] * A[I][91];
    aux2 = A[I][77] * A[I][90] - A[I][78] * A[I][89];
    aux3 = A[I][76] * A[I][91] - A[I][79] * A[I][88];
    aux4 = A[I][78] * A[I][88] - A[I][76] * A[I][90];
    aux5 = A[I][76] * A[I][89] - A[I][77] * A[I][88];

    denominator = (A[I][52] * (A[I][65] * (aux0) + A[I][66] * (aux1) + A[I][67] * (aux2)) +
		               A[I][53] * (A[I][64] * (-aux0) + A[I][66] * (aux3) + A[I][67] * (aux4)) +
		               A[I][54] * (A[I][64] * (-aux1) + A[I][65] * (-aux3) + A[I][67] * (aux5)) +
		               A[I][55] * (A[I][64] * (-aux2) + A[I][65] * (-aux4) + A[I][66] * (-aux5)));

    MatrixData->invDe[I][1] = 1.0 / denominator;

    aux0 = A[I][130] * A[I][143] - A[I][131] * A[I][142];
    aux1 = A[I][131] * A[I][141] - A[I][129] * A[I][143];
    aux2 = A[I][129] * A[I][142] - A[I][130] * A[I][141];
    aux3 = A[I][128] * A[I][143] - A[I][131] * A[I][140];
    aux4 = A[I][130] * A[I][140] - A[I][128] * A[I][142];
    aux5 = A[I][128] * A[I][141] - A[I][129] * A[I][140];

    denominator = (A[I][104] * (A[I][117] * (aux0) + A[I][118] * (aux1) + A[I][119] * (aux2)) +
		               A[I][105] * (A[I][116] * (-aux0) + A[I][118] * (aux3) + A[I][119] * (aux4)) +
		               A[I][106] * (A[I][116] * (-aux1) + A[I][117] * (-aux3) + A[I][119] * (aux5)) +
		               A[I][107] * (A[I][116] * (-aux2) + A[I][117] * (-aux4) + A[I][118] * (-aux5)));

    MatrixData->invDe[I][2] = 1.0 / denominator;
	}

	/* F preconditioning */
	double *faux = mycalloc("faux of JacobiBlockDOF4_precond_EBE_setup",(neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	JacobiBlockDOF4_precond_EBE (Parameters, MatrixData, FemStructs, faux, F);

	myfree(faux);

	return 0;
}
