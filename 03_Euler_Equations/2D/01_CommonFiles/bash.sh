explosion() {
echo "Exp_$precond	:Experiments
EXPLOSION	:Problem title
1e-6	:Solver tolerance
1e-2	:Correction tolerance
1e-1	:Time integration tolerance
1e-8	:Stabilization coefficient tolerance
1000	:Maximum number of solver iterations
30	:Number of vectors in krylov basis
GMRES	:Solver
$precond	:Preconditioner
NOT	:Scaling (NOT or BlockScaling)	
NOT	:Reordering
Predictor1	:Time integration method
NO	:Stop at steady state
0.5	:Alpha
0.001	:Step time
0.25	:Final time
NO	:Dimensionless
1.0	:Mach of reference
10	:Number of multicorrection
ITERATION	:Stop multicorrection
CSR	:Matrix vector product scheme
SUPG	:Stabilization form
YZBeta	:Operator captures
1.0	:invY[0]=1/U1
0.0	:invY[1]=1/U2
0.0	:invY[2]=1/U3
0.4	:invY[3]=1/U4
$nnodes	:Number of nodes
$nel	:Number of elements"
}

#EXPLOSION
precond="ILU1"
mkdir -p 0_Parameters/$precond
for tam in 1 2 3; do
	case $tam in
		1) nnodes=3434
		   nel=6666 ;;
		2) nnodes=8572
		   nel=16822 ;;
		3) nnodes=14954
		   nel=29482 ;;
	esac
	exe="EulerEquations2D"
	dat="0_Parameters/$precond/Exp_${precond}_${nnodes}_${nel}.dat"
	output="../../../../OUTPUT_DATA/"
	newoutput=${output}/Exp_${precond}_${nnodes}_${nel}
	explosion > $dat
	./$exe $dat > $output/tempo.txt
	mkdir -p $newoutput
	find $output -maxdepth 1 -type f -exec mv {} $newoutput \;
done
