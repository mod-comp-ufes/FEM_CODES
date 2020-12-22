preconditioner="AMG 1 1 0.2 1.5 2 0 0 0 0"
precond="AMG_1_0.2"


explosion() {
echo "${experiment}_$precond	:Experiments
${experiment}	:Problem title
1e-6	:Solver tolerance
1e-2	:Correction tolerance
1e-1	:Time integration tolerance
1e-8	:Stabilization coefficient tolerance
1000	:Maximum number of solver iterations
30	:Number of vectors in krylov basis
GMRES	:Solver
$preconditioner	:Preconditioner
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

experiment="EXPLOSION"
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
	dat="0_Parameters/$precond/${experiment}_${precond}_${nnodes}_${nel}.dat"
	output="../../../../OUTPUT_DATA/"
	newoutput=${output}/${experiment}_${precond}_${nnodes}_${nel}
	explosion > $dat
	./$exe $dat > $output/tempo.txt
	mkdir -p $newoutput
	find $output -maxdepth 1 -type f -exec mv {} $newoutput \;
done



OBLIQUO() {
echo "${experiment}_$precond	:Experiments
${experiment}	:Problem title
1e-6	:Solver tolerance
1e-2	:Correction tolerance
1e-1	:Time integration tolerance
1e-8	:Stabilization coefficient tolerance
1000	:Maximum number of solver iterations
30	:Number of vectors in krylov basis
GMRES	:Solver
$preconditioner	:Preconditioner
NOT	:Scaling (NOT or BlockScaling)	
NOT	:Reordering
Predictor1	:Time integration method
NO	:Stop at steady state
0.5	:Alpha
0.001	:Step time
3.0	:Final time
NO	:Dimensionless
2.0	:Mach of reference
10	:Number of multicorrection
ITERATION	:Stop multicorrection
CSR	:Matrix vector product scheme
SUPG	:Stabilization form
YZBeta	:Operator captures
1.0	:invY[0]=1/U1
1.015426612	:invY[1]=1/U2
-5.758770505	:invY[2]=1/U3
1.056607762	:invY[3]=1/U4
$nnodes	:Number of nodes
$nel	:Number of elements"
}

experiment="OBLIQUO"
mkdir -p 0_Parameters/$precond
for tam in 1 2 3; do
	case $tam in
		1) nnodes=3436
		   nel=6670 ;;
		2) nnodes=8558
		   nel=16794 ;;
		3) nnodes=15304
		   nel=30178 ;;
	esac
	exe="EulerEquations2D"
	dat="0_Parameters/$precond/${experiment}_${precond}_${nnodes}_${nel}.dat"
	output="../../../../OUTPUT_DATA/"
	newoutput=${output}/${experiment}_${precond}_${nnodes}_${nel}
	OBLIQUO > $dat
	./$exe $dat > $output/tempo.txt
	mkdir -p $newoutput
	find $output -maxdepth 1 -type f -exec mv {} $newoutput \;
done



REFLETIDO() {
echo "${experiment}_$precond	:Experiments
${experiment}	:Problem title
1e-6	:Solver tolerance
1e-2	:Correction tolerance
1e-1	:Time integration tolerance
1e-8	:Stabilization coefficient tolerance
1000	:Maximum number of solver iterations
30	:Number of vectors in krylov basis
GMRES	:Solver
$preconditioner	:Preconditioner
NOT	:Scaling (NOT or BlockScaling)	
NOT	:Reordering
Predictor1	:Time integration method
NO	:Stop at steady state
0.5	:Alpha
0.001	:Step time
3.0	:Final time
NO	:Dimensionless
2.9	:Mach of reference
10	:Number of multicorrection
ITERATION	:Stop multicorrection
CSR	:Matrix vector product scheme
SUPG	:Stabilization form
YZBeta	:Operator captures
1.0		:invY[0]=1/U1
0.344827586	:invY[1]=1/U2
0.0		:invY[2]=1/U3
0.166924983	:invY[3]=1/U4
$nnodes	:Number of nodes
$nel	:Number of elements"
}

experiment="REFLETIDO"
mkdir -p 0_Parameters/$precond
for tam in 1 2 3; do
	case $tam in
		1) nnodes=3523
		   nel=6788 ;;
		2) nnodes=8895
		   nel=17380 ;;
		3) nnodes=15263
		   nel=29986 ;;
	esac
	exe="EulerEquations2D"
	dat="0_Parameters/$precond/${experiment}_${precond}_${nnodes}_${nel}.dat"
	output="../../../../OUTPUT_DATA/"
	newoutput=${output}/${experiment}_${precond}_${nnodes}_${nel}
	REFLETIDO > $dat
	./$exe $dat > $output/tempo.txt
	mkdir -p $newoutput
	find $output -maxdepth 1 -type f -exec mv {} $newoutput \;
done



TUNEL() {
echo "${experiment}_$precond	:Experiments
${experiment}	:Problem title
1e-6	:Solver tolerance
1e-2	:Correction tolerance
1e-1	:Time integration tolerance
1e-8	:Stabilization coefficient tolerance
1000	:Maximum number of solver iterations
30	:Number of vectors in krylov basis
GMRES	:Solver
$preconditioner	:Preconditioner
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
0.714285714	:invY[0]=1/U1
0.238095238	:invY[1]=1/U2
0.0		:invY[2]=1/U3
0.113636363	:invY[3]=1/U4
$nnodes	:Number of nodes
$nel	:Number of elements"
}

experiment="TUNEL"
mkdir -p 0_Parameters/$precond
for tam in 1 2 3; do
	case $tam in
		1) nnodes=3806
		   nel=7342 ;;
		2) nnodes=8549
		   nel=16696 ;;
		3) nnodes=15112
		   nel=29687 ;;
	esac
	exe="EulerEquations2D"
	dat="0_Parameters/$precond/${experiment}_${precond}_${nnodes}_${nel}.dat"
	output="../../../../OUTPUT_DATA/"
	newoutput=${output}/${experiment}_${precond}_${nnodes}_${nel}
	TUNEL > $dat
	./$exe $dat > $output/tempo.txt
	mkdir -p $newoutput
	find $output -maxdepth 1 -type f -exec mv {} $newoutput \;
done
