preconditioner="ILU1"
precond="ILU1"

explosion() {
echo "${experiment}_$precond	: Experiments
${experiment}	: Problem Title
1.2	: Reynolds Number
1.2e-2	: Vel Max
1e-8	: Solver Tolerance
0	: SolverToleranceCase
5e-4	: SolverToleranceNonLin
1000	: LinearMaxIter
45	: KrylovBasisVectorsQuantity
GMRES	: Solver
$preconditioner	: Preconditioner
NOT	: reordering
CSR	: MatrixVectorProductScheme
SUPG	: StabilizationForm
$nnodes	: nnodes
$nel	: nel"
}

experiment="EXATA"
mkdir -p 0_Parameters/$precond
for tam in 1; do
	case $tam in
		1) nnodes=12
		   nel=12 ;;
	esac
	exe="SSNavierStokesEquations2D"
	dat="0_Parameters/$precond/${experiment}_${precond}_${nnodes}_${nel}.dat"
	output="../03_output/"
	newoutput=${output}/${experiment}_${precond}_${nnodes}_${nel}
	explosion > $dat
	./$exe $dat > $output/tempo.txt
	mkdir -p $newoutput
	find $output -maxdepth 1 -type f -exec mv {} $newoutput \;
done