#include "../preconditioners.h"
/* p = A z
 * A = (D + U)^{-1} (D + L)^{-1}
 * (L + D) * p = w
 * (D + U) * w = z
 * z out
*/
int SSORBlockDOF4_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z) {

	int I;
	int lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11; //auxiliar
	double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11; //auxiliar
	double w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11; //auxiliar
	double p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data
	double **inv = MatrixData->invDe; //1 / (1 + D)

	double omega;
	omega = atof(&(Parameters->Preconditioner[8])); // In SORw, w is the weight coefficient
	
	for (I = 0; I < neq; I++){
		z[I] = 0;
	}
	z[neq] = 0;

	/* FORMA 1
	 * | x 0 0 0                 |
	 * | 0 x 0 0                 |
	 * | 0 0 x 0                 |
	 * | 0 0 0 x                 |
	 * | x x x x x 0 0 0         |
	 * | x x x x 0 x 0 0         |
	 * | x x x x 0 0 x 0				 |
	 * | x x x x 0 0 0 x				 |
	 * | x x x x x x x x x 0 0 0 |
	 * | x x x x x x x x 0 x 0 0 |
	 * | x x x x x x x x 0 0 x 0 |
	 * | x x x x x x x x 0 0 0 x |
	 */
	 //OPTION 1: (D + wU)(D + wL)
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

	 	// (D + w * L) * w = p
	 	w0 = p0 * inv[I][0];
	 	w1 = p1 * inv[I][1];
	 	w2 = p2 * inv[I][2];
	 	w3 = p3 * inv[I][3];

	 	w4 = (p4 - omega * A[I][48] * w0 - omega * A[I][49] * w1 - omega * A[I][50] * w2 - omega * A[I][51] * w3) * inv[I][4];
	 	w5 = (p5 - omega * A[I][60] * w0 - omega * A[I][61] * w1 - omega * A[I][62] * w2 - omega * A[I][63] * w3) * inv[I][5];
	 	w6 = (p6 - omega * A[I][72] * w0 - omega * A[I][73] * w1 - omega * A[I][74] * w2 - omega * A[I][75] * w3) * inv[I][6];
	 	w7 = (p7 - omega * A[I][84] * w0 - omega * A[I][85] * w1 - omega * A[I][86] * w2 - omega * A[I][87] * w3) * inv[I][7];

	 	w8 = (p8 - omega * A[I][96] * w0 - omega * A[I][97] * w1 - omega * A[I][98] * w2 - omega * A[I][99] * w3 - omega * A[I][100] * w4 - omega * A[I][101] * w5 - omega * A[I][102] * w6 - omega * A[I][103] * w7) * inv[I][8];
	 	w9 = (p9 - omega * A[I][108] * w0 - omega * A[I][109] * w1 - omega * A[I][110] * w2 - omega * A[I][11] * w3 - omega * A[I][112] * w4 - omega * A[I][113] * w5 - omega * A[I][114] * w6 - omega * A[I][115] * w7) * inv[I][9];
	 	w10 = (p10 - omega * A[I][120] * w0 - omega * A[I][121] * w1 - omega * A[I][122] * w2 - omega * A[I][123] * w3 - omega * A[I][124] * w4 - omega * A[I][125] * w5 - omega * A[I][126] * w6 - omega * A[I][127] * w7) * inv[I][10];
	 	w11 = (p11 - omega * A[I][132] * w0 - omega * A[I][133] * w1 - omega * A[I][134] * w2 - omega * A[I][135] * w3 - omega * A[I][136] * w4 - omega * A[I][137] * w5 - omega * A[I][138] * w6 - omega * A[I][139] * w7) * inv[I][11];

	 	//(D + w U) * z = w (
	 	z11 = w11 * inv[I][11];
	 	z10 = w10 * inv[I][10];
	 	z9 = w9 * inv[I][9];
	 	z8 = w8 * inv[I][8];

	 	z7 = (w7 - omega * A[I][92] * z8 - omega * A[I][93] * z9 - omega * A[I][94] * z10 - omega * A[I][95] * z11) * inv[I][7];
	 	z6 = (w6 - omega * A[I][80] * z8 - omega * A[I][81] * z9 - omega * A[I][82] * z10 - omega * A[I][83] * z11) * inv[I][6];
	 	z5 = (w5 - omega * A[I][68] * z8 - omega * A[I][69] * z9 - omega * A[I][70] * z10 - omega * A[I][71] * z11) * inv[I][5];
	 	z4 = (w4 - omega * A[I][56] * z8 - omega * A[I][57] * z9 - omega * A[I][58] * z10 - omega * A[I][59] * z11) * inv[I][4];

	 	z3 = (w3 - omega * A[I][40] * z4 - omega * A[I][41] * z5 - omega * A[I][42] * z6 - omega * A[I][43] * z7 - omega * A[I][44] * z8 - omega * A[I][45] * z9 - omega * A[I][46] * z10 - omega * A[I][47] * z11) * inv[I][3];
	 	z2 = (w2 - omega * A[I][28] * z4 - omega * A[I][29] * z5 - omega * A[I][30] * z6 - omega * A[I][31] * z7 - omega * A[I][32] * z8 - omega * A[I][33] * z9 - omega * A[I][34] * z10 - omega * A[I][35] * z11) * inv[I][2];
	 	z1 = (w1 - omega * A[I][16] * z4 - omega * A[I][17] * z5 - omega * A[I][18] * z6 - omega * A[I][19] * z7 - omega * A[I][20] * z8 - omega * A[I][21] * z9 - omega * A[I][22] * z11 - omega * A[I][23] * z11) * inv[I][1];
	 	z0 = (w0 - omega * A[I][4] * z4 - omega * A[I][5] * z5 - omega * A[I][6] * z6 - omega * A[I][7] * z7 - omega * A[I][8] * z8 - omega * A[I][9] * z9 - omega * A[I][10] * z10 - omega * A[I][11] * z11) * inv[I][0];

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

	 // OPTION 2: (D + wU) D^{-1} (D + wL)
	 /*
	 double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11; //auxiliar
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

	 	// (D + w * L) * w = p
	 	w0 = p0 * inv[I][0];
	 	w1 = p1 * inv[I][1];
	 	w2 = p2 * inv[I][2];
	 	w3 = p3 * inv[I][3];

	 	w4 = (p4 - omega * A[I][48] * w0 - omega * A[I][49] * w1 - omega * A[I][50] * w2 - omega * A[I][51] * w3) * inv[I][4];
	 	w5 = (p5 - omega * A[I][60] * w0 - omega * A[I][61] * w1 - omega * A[I][62] * w2 - omega * A[I][63] * w3) * inv[I][5];
	 	w6 = (p6 - omega * A[I][72] * w0 - omega * A[I][73] * w1 - omega * A[I][74] * w2 - omega * A[I][75] * w3) * inv[I][6];
	 	w7 = (p7 - omega * A[I][84] * w0 - omega * A[I][85] * w1 - omega * A[I][86] * w2 - omega * A[I][87] * w3) * inv[I][7];

	 	w8 = (p8 - omega * A[I][96] * w0 - omega * A[I][97] * w1 - omega * A[I][98] * w2 - omega * A[I][99] * w3 - omega * A[I][100] * w4 - omega * A[I][101] * w5 - omega * A[I][102] * w6 - omega * A[I][103] * w7) * inv[I][8];
	 	w9 = (p9 - omega * A[I][108] * w0 - omega * A[I][109] * w1 - omega * A[I][110] * w2 - omega * A[I][11] * w3 - omega * A[I][112] * w4 - omega * A[I][113] * w5 - omega * A[I][114] * w6 - omega * A[I][115] * w7) * inv[I][9];
	 	w10 = (p10 - omega * A[I][120] * w0 - omega * A[I][121] * w1 - omega * A[I][122] * w2 - omega * A[I][123] * w3 - omega * A[I][124] * w4 - omega * A[I][125] * w5 - omega * A[I][126] * w6 - omega * A[I][127] * w7) * inv[I][10];
	 	w11 = (p11 - omega * A[I][132] * w0 - omega * A[I][133] * w1 - omega * A[I][134] * w2 - omega * A[I][135] * w3 - omega * A[I][136] * w4 - omega * A[I][137] * w5 - omega * A[I][138] * w6 - omega * A[I][139] * w7) * inv[I][11];

	 	//D^{-1} x = w
	 	x0 = A[I][0] * w0;
	 	x1 = A[I][13] * w1;
	 	x2 = A[I][26] * w2;
	 	x3 = A[I][39] * w3;

	 	x4 = A[I][52] * w4;
	 	x5 = A[I][65] * w5;
	 	x6 = A[I][78] * w6;
	 	x7 = A[I][91] * w7;

	 	x8 = A[I][104] * w8;
	 	x9 = A[I][117] * w9;
	 	x10 = A[I][130] * w10;
	 	x11 = A[I][143] * w11;

	 	//(D + w U) * z = x
	 	z11 = x11 * inv[I][11];
	 	z10 = x10 * inv[I][10];
	 	z9 = x9 * inv[I][9];
	 	z8 = x8 * inv[I][8];

	 	z7 = (x7 - omega * A[I][92] * z8 - omega * A[I][93] * z9 - omega * A[I][94] * z10 - omega * A[I][95] * z11) * inv[I][7];
	 	z6 = (x6 - omega * A[I][80] * z8 - omega * A[I][81] * z9 - omega * A[I][82] * z10 - omega * A[I][83] * z11) * inv[I][6];
	 	z5 = (x5 - omega * A[I][68] * z8 - omega * A[I][69] * z9 - omega * A[I][70] * z10 - omega * A[I][71] * z11) * inv[I][5];
	 	z4 = (x4 - omega * A[I][56] * z8 - omega * A[I][57] * z9 - omega * A[I][58] * z10 - omega * A[I][59] * z11) * inv[I][4];

	 	z3 = (x3 - omega * A[I][40] * z4 - omega * A[I][41] * z5 - omega * A[I][42] * z6 - omega * A[I][43] * z7 - omega * A[I][44] * z8 - omega * A[I][45] * z9 - omega * A[I][46] * z10 - omega * A[I][47] * z11) * inv[I][3];
	 	z2 = (x2 - omega * A[I][28] * z4 - omega * A[I][29] * z5 - omega * A[I][30] * z6 - omega * A[I][31] * z7 - omega * A[I][32] * z8 - omega * A[I][33] * z9 - omega * A[I][34] * z10 - omega * A[I][35] * z11) * inv[I][2];
	 	z1 = (x1 - omega * A[I][16] * z4 - omega * A[I][17] * z5 - omega * A[I][18] * z6 - omega * A[I][19] * z7 - omega * A[I][20] * z8 - omega * A[I][21] * z9 - omega * A[I][22] * z11 - omega * A[I][23] * z11) * inv[I][1];
	 	z0 = (x0 - omega * A[I][4] * z4 - omega * A[I][5] * z5 - omega * A[I][6] * z6 - omega * A[I][7] * z7 - omega * A[I][8] * z8 - omega * A[I][9] * z9 - omega * A[I][10] * z10 - omega * A[I][11] * z11) * inv[I][0];

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
	 */

	return 0;
}
