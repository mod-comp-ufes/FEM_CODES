#include "../preconditioners.h"
/* p = A z
 * A = D^{-1}
 * D * p = z
 * z out
*/
int JacobiBlockDOF4_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z) {
	int I;
	int lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11; //auxiliar
	double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11; //auxiliar
	double p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double aux0, aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8, aux9; //auxiliar


	for (I = 0; I < neq; I++){
		z[I] = 0;
	}
	z[neq] = 0;

	//OPTION 2: 4X4 INV
	double numerator;
	double **A = MatrixData->A;
	double **denominator = MatrixData->invDe;

	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
		lm1 = lm[I][1];
		lm2 = lm[I][2];
		lm3 = lm[I][3];
		lm4 = lm[I][4];
		lm5 = lm[I][5];
		lm6 = lm[I][6];
		lm7 = lm[I][7];
		lm8 = lm[I][8];
		lm9 = lm[I][9];
		lm10 = lm[I][10];
		lm11 = lm[I][11];

		p0 = p[lm0];
		p1 = p[lm1];
		p2 = p[lm2];
		p3 = p[lm3];
		p4 = p[lm4];
		p5 = p[lm5];
		p6 = p[lm6];
		p7 = p[lm7];
		p8 = p[lm8];
		p9 = p[lm9];
		p10 = p[lm10];
		p11 = p[lm11];

		// BLOCO 1
		// | A0 A1 A2 A3 |
		// | A12 A13 A14 A15 |
		// | A24 A25 A26 A27 |
		// | A36 A37 A38 A39 |
		aux0 = A[I][27] * p3 - A[I][39] * p2;
		aux1 = A[I][38] * p2 - A[I][26] * p3;
		aux2 = A[I][26] * A[I][39] - A[I][27] * A[I][38];
		aux3 = A[I][25] * p3 - A[I][37] * p2;
		aux4 = A[I][27] * A[I][37] - A[I][25] * A[I][39];
		aux5 = A[I][25] * A[I][38] - A[I][26] * A[I][37];
		aux6 = A[I][24] * p3 - A[I][36] * p2;
		aux7 = A[I][27] * A[I][36] - A[I][24] * A[I][39];
		aux8 = A[I][24] * A[I][38] - A[I][26] * A[I][36];
		aux9 = A[I][24] * A[I][37] - A[I][25] * A[I][36];

		numerator = -(A[I][1] * (A[I][14] * (aux0) + A[I][15] * (aux1) + (aux2) * p1) +
		              A[I][2] * (A[I][13] * (-aux0) + A[I][15] * (aux3) + (aux4) * p1) +
		              A[I][3] * (A[I][13] * (-aux1) + A[I][14] * (-aux3) + (aux5) * p1) +
		              (A[I][13] * (-aux2) + A[I][14] * (-aux4) + A[I][15] * (-aux5)) * p0);

		z0 = numerator * denominator[I][0];

		numerator = (A[I][0] * (A[I][14] * (aux0) + A[I][15] * (aux1) + (aux2) * p1) +
		             A[I][2] * (A[I][12] * (-aux0) + A[I][15] * (aux6) + (aux7) * p1) +
		             A[I][3] * (A[I][12] * (-aux1) + A[I][14] * (-aux6) + (aux8) * p1) +
		             (A[I][12] * (-aux2) + A[I][14] * (-aux7) + A[I][15] * (-aux8)) * p0);

		z1 = numerator * denominator[I][0];

		numerator = - (A[I][0] * (A[I][13] * (aux0) + A[I][15] * (-aux3) + (-aux4) * p1) +
		               A[I][1] * (A[I][12] * (-aux0) + A[I][15] * (aux6) + (aux7) * p1) +
		               A[I][3] * (A[I][12] * (aux3) + A[I][13] * (-aux6) + (aux9) * p1) +
		               (A[I][12] * (aux4) + A[I][13] * (-aux7) + A[I][15] * (-aux9)) * p0);

		z2 = numerator * denominator[I][0];

		numerator = (A[I][0] * (A[I][13] * (-aux1) + A[I][14] * (-aux3) + (aux5) * p1) +
		             A[I][1] * (A[I][12] * (aux1) + A[I][14] * (aux6) + (-aux8) * p1) +
		             A[I][2] * (A[I][12] * (aux3) + A[I][13] * (-aux6) + (aux9) * p0) +
		             (A[I][12] * (-aux5) + A[I][13] * (aux8) + A[I][14] * (-aux9)) * p0);

		z3 = numerator * denominator[I][0];

		// BLOCO 2
		// | A52 A53 A54 A55 |
		// | A64 A65 A66 A67 |
		// | A76 A77 A78 A79 |
		// | A88 A89 A90 A91 |
		aux0 = A[I][79] * p7 - A[I][91] * p6;
		aux1 = A[I][90] * p6 - A[I][78] * p7;
		aux2 = A[I][78] * A[I][91] - A[I][79] * A[I][90];
		aux3 = A[I][77] * p7 - A[I][89] * p6;
		aux4 = A[I][79] * A[I][89] - A[I][77] * A[I][91];
		aux5 = A[I][77] * A[I][90] - A[I][78] * A[I][89];
		aux6 = A[I][76] * p7 - A[I][88] * p6;
		aux7 = A[I][79] * A[I][88] - A[I][76] * A[I][91];
		aux8 = A[I][76] * A[I][90] - A[I][78] * A[I][88];
		aux9 = A[I][76] * A[I][89] - A[I][77] * A[I][88];

		numerator = - (A[I][53] * (A[I][66] * (aux0) + A[I][67] * (aux1) + (aux2) * p5) +
		               A[I][54] * (A[I][65] * (-aux0) + A[I][67] * (aux3) + (aux4) * p5) +
		               A[I][55] * (A[I][65] * (-aux1) + A[I][66] * (-aux3) + (aux5) * p5) +
		               (A[I][65] * (-aux2) + A[I][66] * (-aux4) + A[I][67] * (-aux5)) * p4);

		z4 = numerator * denominator[I][1];

		numerator = (A[I][52] * (A[I][66] * (aux0) + A[I][67] * (aux1) + (aux2) * p5) +
		             A[I][54] * (A[I][64] * (-aux0) + A[I][67] * (aux6) + (aux7) * p5) +
		             A[I][55] * (A[I][64] * (-aux1) + A[I][66] * (-aux6) + (aux8) * p5) +
		             (A[I][64] * (-aux2) + A[I][66] * (-aux7) + A[I][67] * (-aux8)) * p4);

		z5 = numerator * denominator[I][1];

		numerator = - (A[I][52] * (A[I][65] * (aux0) + A[I][67] * (-aux3) + (-aux4) * p5) +
		               A[I][53] * (A[I][64] * (-aux0) + A[I][67] * (aux6) + (aux7) * p5) +
		               A[I][55] * (A[I][64] * (aux3) + A[I][65] * (-aux6) + (aux9) * p5) +
		               (A[I][64] * (aux4) + A[I][65] * (-aux7) + A[I][67] * (-aux9)) * p4);

		z6 = numerator * denominator[I][1];

		numerator = (A[I][52] * (A[I][65] * (-aux1) + A[I][66] * (-aux3) + (aux5) * p5) +
		             A[I][53] * (A[I][64] * (aux1) + A[I][66] * (aux6) + (-aux8) * p5) +
		             A[I][54] * (A[I][64] * (aux3) + A[I][65] * (-aux6) + (aux9) * p5) +
		             (A[I][64] * (-aux5) + A[I][65] * (aux8) + A[I][66] * (-aux9)) * p4);

		z7 = numerator * denominator[I][1];

		// BLOCO 3
		// | A104 A105 A106 A107 |
		// | A116 A117 A118 A119 |
		// | A128 A129 A130 A131 |
		// | A140 A141 A142 A143 |
		aux0 = A[I][130] * A[I][143] - A[I][131] * A[I][142];
		aux1 = A[I][131] * p11 - A[I][143] * p10;
		aux2 = A[I][142] * p10 - A[I][130] * p11;
		aux3 = A[I][131] * A[I][141] - A[I][129] * A[I][143];
		aux4 = A[I][129] * p11 - A[I][141] * p10;
		aux5 = A[I][129] * A[I][142] - A[I][130] * A[I][141];
		aux6 = A[I][131] * A[I][140] - A[I][128] * A[I][143];
		aux7 = A[I][128] * p11 - A[I][140] * p10;
		aux8 = A[I][128] * A[I][142] - A[I][130] * A[I][140];
		aux9 = A[I][128] * A[I][141] - A[I][129] * A[I][140];

		numerator = - (A[I][105] * ((aux0) * p6 + A[I][118] * (aux1) + A[I][119] * (aux2)) +
		               A[I][106] * ((aux3) * p9 + A[I][117] * (-aux1) + A[I][119] * (aux4)) +
		               A[I][107] * ((aux5) * p9 + A[I][117] * (-aux2) + A[I][118] * (-aux4)) +
		               (A[I][117] * (-aux0) + A[I][118] * (-aux3) + A[I][119] * (-aux5)) * p8);

		z8 = numerator * denominator[I][2];

		numerator = (A[I][104] * ((aux0) * p9 + A[I][118] * (aux1) + A[I][119] * (aux2)) +
		             A[I][106] * ((aux6) * p9 + A[I][116] * (-aux1) + A[I][119] * (aux7)) +
		             A[I][107] * ((aux8) * p9 + A[I][116] * (-aux2) + A[I][118] * (-aux7)) +
		             (A[I][116] * (-aux0) + A[I][118] * (-aux6) + A[I][119] * (-aux8)) * p8);

		z9 = numerator * denominator[I][2];

		numerator = - (A[I][104] * ((-aux3) * p9 + A[I][117] * (aux1) + A[I][119] * (-aux4)) +
		               A[I][105] * ((aux6) * p9 + A[I][116] * (-aux1) + A[I][119] * (aux7)) +
		               A[I][107] * ((aux9) * p9 + A[I][116] * (aux4) + A[I][117] * (-aux7)) +
		               (A[I][116] * (aux3) + A[I][117] * (aux6) + A[I][119] * (-aux9)) * p8);

		z10 = numerator * denominator[I][2];

		numerator = (A[I][104] * ((aux5) * p9 + A[I][117] * (-aux2) + A[I][118] * (-aux4)) +
		             A[I][105] * ((-aux8) * p9 + A[I][116] * (aux2) + A[I][118] * (aux7)) +
		             A[I][106] * ((aux9) * p9 + A[I][116] * (aux4) + A[I][117] * (-aux7)) +
		             (A[I][116] * (-aux5) + A[I][117] * (aux8) + A[I][118] * (-aux9)) * p8);

		z11 = numerator * denominator[I][2];

		//ASSEMBLY
		z[lm0] += z0;
		z[lm1] += z1;
		z[lm2] += z2;
		z[lm3] += z3;
		z[lm4] += z4;
		z[lm5] += z5;
		z[lm6] += z6;
		z[lm7] += z7;
		z[lm8] += z8;
		z[lm9] += z9;
		z[lm10] += z10;
		z[lm11] += z11;

	}
	z[neq] = 0.0;

  return 0;
}
