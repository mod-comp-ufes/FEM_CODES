# Common files

Implementação de funções comuns aos métodos numéricos utilizados.

#### Allocation
Modularização das operações de alocação de memória `calloc()` e `free()`.

#### BLAS
Operações matrix-vetor-escalar.

#### C
Operações de cópia e definir zero.

#### CPU Time
Cálculo de tempo de simulações, de acordo com o problema definido. Problemas definidos: `SSTranspEquation2D`, `SSTransportEquation3D`, `TranspEquation2D`, `EulerEquations2D`, `SSNavierStokesEquations2D`, `SSNavierStokesEquations3D`, `NavierStokesEquations2D`, `NavierStokesEquations3D`, `SS_StokesEquations3D`, `SS7_StokesEquations3D`, `PoissonEquation3D`.

#### IO
Modularização da função de abertura de arquivo `fopen()`.

#### Matrix Vector
Operações de matriz vetor de acordo com a respresentação utilizada (CSR, EBE, EDE).

#### Reordering
Algoritmos de reordenamento de matrizes: Cuthill-McKee (SYMRCM) e Spectral Ordering.

#### Solvers and Precondicioners
Implementação de solvers e precondicionadores de acordo com a representação utilizada (CSR, EBE, EDE).
Solvers: GMRES com e sem precondicionamento, PCG (non-linear)
Precondicionadores: Diag, ILUp.

#### Time integration
