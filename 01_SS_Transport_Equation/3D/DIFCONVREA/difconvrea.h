
#include "../01_CommonFiles/SSTransportEquation3D.h"
#include <math.h>

double DIFCONVREA_upresc(double, double, double, double, double);
double DIFCONVREA_f(double, double, double, double);
double DIFCONVREA_ExactSolution(double, double, double, double Cst);
void DIFCONVREA_ExactSolutionAllPoints(NodeType *, double, double *, double);
void DIFCONVREA_kappa(double *, double *, double *, double, double, double, double);
void DIFCONVREA_beta(double *, double *, double *, double, double, double, double, double, double, double);
double DIFCONVREA_sigma(double, double, double);
double DIFCONVREA_DuDx(double, double, double, double Cst);
double DIFCONVREA_DuDy(double, double, double, double Cst);
double DIFCONVREA_DuDz(double, double, double, double Cst);
