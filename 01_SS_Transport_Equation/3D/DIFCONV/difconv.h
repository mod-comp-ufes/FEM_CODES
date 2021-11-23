# include "../01_CommonFiles/SSTransportEquation3D.h"

double DIFCONV_upresc(double, double, double, double, double);
double DIFCONV_f(double, double, double, double);
double DIFCONV_ExactSolution(double, double, double, double);
void DIFCONV_ExactSolutionAllPoints(NodeType *, double, double *, double);
void DIFCONV_kappa(double *, double *, double *, double, double, double, double);
void DIFCONV_beta(double *, double *, double *, double, double, double, double, double, double, double);
double DIFCONV_sigma(double, double, double);
double DIFCONV_DuDx(double, double, double, double);
double DIFCONV_DuDy(double, double, double, double);
double DIFCONV_DuDz(double, double, double, double);
