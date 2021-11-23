// Retirada do artigo 'Galerkin and Least Squares Methods to solve a 3D convection-diffusion-reaction equation with variable coefficients' Romão e Moura, 2012

# include "../01_CommonFiles/SSTransportEquation3D.h"

double APLI6_upresc(double, double, double, double, double);
double APLI6_f(double, double, double, double);
double APLI6_ExactSolution(double x, double y, double z, double A);
void APLI6_ExactSolutionAllPoints(NodeType *Node, double nnodes, double *u, double A);
void APLI6_kappa(double *, double *, double *, double, double, double, double);
void APLI6_beta(double *, double *, double *, double, double, double, double, double, double, double);
double APLI6_sigma(double, double, double);
double APLI6_DuDx(double x, double y, double z, double A);
double APLI6_DuDy(double x, double y, double z, double A);
double APLI6_DuDz(double x, double y, double z, double A);
