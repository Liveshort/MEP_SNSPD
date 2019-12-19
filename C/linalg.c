#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

// function that solves a tridiagonal matrix
//     input arguments are the target vector, the dimension of the problem (everthing is square
//     or a vector), the diagonal of the matrix, the off diagonal of the matrix and the right
//     hand side vector that we want to solve against. This is an implementation of the Thomas
//     algorithm.
int TDM_solve(double * trg, size_t J, double * lhsDiag_n, double * lhsOffDiag_n, double * rhs) {
    // implementation of the Thomas algorithm
    // determine gam(ma) from the implementation on [http://www.industrial-maths.com/ms6021_thomas.pdf]
    // note that the gamma from the theory is put in the off diagonal, which we won't need later on
    // this saves memory allocation

    // determine the first round (which is different from the rest)
    rhs[0] = rhs[0]/lhsDiag_n[0];
    lhsOffDiag_n[0] = lhsOffDiag_n[0]/lhsDiag_n[0];
    // perform the downward sweep
    for (unsigned j=1; j<J-1; ++j) {
        rhs[j] = (rhs[j] - lhsOffDiag_n[j]*rhs[j - 1])/(lhsDiag_n[j] - lhsOffDiag_n[j]*lhsOffDiag_n[j-1]);
        lhsOffDiag_n[j] = lhsOffDiag_n[j]/(lhsDiag_n[j] - lhsOffDiag_n[j]*lhsOffDiag_n[j-1]);
    }
    // calculate the last row seperately again
    rhs[J - 1] = (rhs[J - 1] - lhsOffDiag_n[J - 1]*rhs[J - 2])/(lhsDiag_n[J - 1] - lhsOffDiag_n[J - 1]*lhsOffDiag_n[J - 2]);

    // the last equation is already solved, so input that answer directly into the target array
    trg[J - 1] = rhs[J - 1];
    // perform the upward sweep
    for (size_t j=J-2; j != (size_t) -1; --j) {
        trg[j] = rhs[j] - lhsOffDiag_n[j]*trg[j + 1];
    }

    return 0;
}
