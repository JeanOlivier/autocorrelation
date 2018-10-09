# autocorrelation.{c, py}
C and Python interfaces to compute the autocorrelation vs. lag of a signal. 

## Description
`autocorrelation.c` provides functions to compute the autocorrelation of a signal for the *k* first lags/delays.
`autocorrelation.py` is a python interface for the C code.

## Algorithm and optimisation
The algorithm used is the following:

![](https://latex.codecogs.com/gif.latex?a_k%20%3D%20%5Cfrac%7B1%7D%7B%28N-k%29%5Csigma%5E2%7D%5Csum_%7Bi%3D1%7D%5E%7BN-k%7D%28x_i-%5Clangle%20x%5Crangle%29%28x_%7Bi&plus;k%7D-%5Clangle%20x%5Crangle%29),

where ![](https://latex.codecogs.com/gif.latex?x_k) is the ![](https://latex.codecogs.com/gif.latex?k^\text{th}) point of the signal, ![](https://latex.codecogs.com/gif.latex?N) is its length, ![](https://latex.codecogs.com/gif.latex?\langle%20x\rangle) is its mean, and ![](https://latex.codecogs.com/gif.latex?\sigma^2) is its variance.

We re-express it as :

![](https://latex.codecogs.com/gif.latex?a_k%20%3D%20%5CBigg%28%20r_k%20-%20%5Cfrac%7BM%5E2%7D%7BN%7D%20&plus;%20%5Cfrac%7BM%7D%7BN%7D%5Cleft%28%5Cbeta_k&plus;%5Cgamma_k%5Cright%29%20-%20%5Cfrac%7BM%5E2k%7D%7BN%5E2%7D%20%5CBigg%29%5CBigg/%5CBigg%28%20%5Cleft%28N-k%5Cright%29%5Csigma%5E2%20%5CBigg%29),

with 

![](https://latex.codecogs.com/gif.latex?r_k%20%3D%20%5Csum%5E%7BN-k%7D_%7Bi%3D1%7D%20x_ix_%7Bi&plus;k%7D%20%5Cquad%3B%5Cquad%20%5Cbeta_k%20%3D%20%5Csum_%7Bi%3DN-k&plus;1%7D%5E%7BN%7D%20x_i%20%5Cquad%3B%5Cquad%20%5Cgamma_k%20%3D%20%5Csum_%7Bi%3D1%7D%5E%7Bk%7D%20x_i%20%5Cquad%5Ctext%7Band%7D%5Cquad%20M%3D%5Csum_%7Bi%3D1%7D%5EN%20x_i),

in order to speed up calculations with simple parallelizable forms for ![](https://latex.codecogs.com/gif.latex?r_k) and ![](https://latex.codecogs.com/gif.latex?M) and correcting for the induced errors using the relatively simple ![](https://latex.codecogs.com/gif.latex?\gamma_k) and ![](https://latex.codecogs.com/gif.latex?\beta_k).

In the ![](https://latex.codecogs.com/gif.latex?k=1) case, it falls back to the variance normalized by itself (a.k.a. `1`) and the corrections vanish.

## Other details

- It can comput the autocorrelations using the whole bytes or for each bits of every byte.
- The whole input file is loaded into memory, but it's otherwise trying to be memory efficient.
- `MPFR` is used for precision when working with big numbers.
- `openMP` is used for parallelisation; some of the code might be suboptimal for single threaded execution.

## Compiling
- To compile a binary: `gcc -O3 autocorrelation.c -o autocorrelation.out -Wall -lmpfr -lgmp -fopenmp`
- To compile a shared library: `gcc -O3 -fPIC -shared autocorrelation.c -o autocorrelation.so -Wall -lmpfr -lgmp -fopenmp`

For linux, add `-DUNIX`.
