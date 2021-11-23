
# include "../01_CommonFiles/SSTransportEquation3D.h"

double APLI4_upresc(double, double, double, double, double);
double APLI4_f(double, double, double, double);
double APLI4_ExactSolution(double x, double y, double z, double A);
void APLI4_ExactSolutionAllPoints(NodeType *Node, double nnodes, double *u, double A);
void APLI4_kappa(double *, double *, double *, double, double, double, double);
void APLI4_beta(double *, double *, double *, double, double, double, double, double, double, double);
double APLI4_sigma(double, double, double);
double APLI4_DuDx(double x, double y, double z, double A);
double APLI4_DuDy(double x, double y, double z, double A);
double APLI4_DuDz(double x, double y, double z, double A);
