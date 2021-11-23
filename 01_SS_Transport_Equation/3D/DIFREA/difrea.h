# include "../01_CommonFiles/SSTransportEquation3D.h"

double DIFREA_upresc(double, double, double, double, double);
double DIFREA_f(double, double, double, double);
double DIFREA_ExactSolution(double, double, double, double Cst);
void DIFREA_ExactSolutionAllPoints(NodeType *, double, double *, double);
void DIFREA_kappa(double *, double *, double *, double, double, double, double);
void DIFREA_beta(double *, double *, double *, double, double, double, double, double, double, double);
double DIFREA_sigma(double, double, double);
double DIFREA_DuDx(double, double, double, double Cst);
double DIFREA_DuDy(double, double, double, double Cst);
double DIFREA_DuDz(double, double, double, double Cst);
