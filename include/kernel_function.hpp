#ifndef __kernel_function__
#define __kernel_function__

#include "myHeaders.hpp"
#include "points_dt.hpp"

template <class Operator>
class kernel_function
{
    Operator *kernel;
public :
    kernel_function(Operator *A){
        this->kernel = A;
    }
    // Returns individual entries of the matrix:
    // this function can be complemented through a user kernel class
    // try to use the same function, expected the kernel evaluation at ith and jth gridPoints
    // Essentially it queries into K_{ij}
    Vec getRow(int j, int n_col_start, int n_cols)
    {
        Vec row(n_cols);
        // #pragma omp parallel for
        for (int k = 0; k < n_cols; k++)
            row(k) = kernel->getMatrixEntry(j, k + n_col_start);
        return row;
    }

Vec getCol(int k, int n_row_start, int n_rows){
    Vec col(n_rows);
    // #pragma omp parallel for
    for (int j = 0; j < n_rows; ++j)
            col(j) = kernel->getMatrixEntry(j + n_row_start, k);
    return col;
}

Vec getDiag1(const int n_row_start, const int n_col_start, const int n_rows, const int n_cols){
    int N = std::max(n_rows, n_cols);
    Vec diag(N);
    int row_ind, col_ind;
    // #pragma omp parallel for
    for (int j = 0; j < N; ++j){
        if (n_cols > n_rows){
            row_ind = mod(n_row_start - n_col_start + j, n_rows);
            col_ind = j;
        }
        else{
            row_ind = j;
            col_ind = mod(n_col_start - n_row_start + j, n_cols);
        }
        diag(j) = kernel->getMatrixEntry(row_ind, col_ind);
    }
    return diag;
}

Vec getDiag2(const int n_row_start, const int n_col_start, const int n_rows, const int n_cols){
    int N = std::max(n_rows, n_cols);
    Vec diag(N);
    int row_ind, col_ind;
    // #pragma omp parallel for
    for (int j = 0; j < N; ++j){
        if (n_cols > n_rows){
            row_ind = mod(n_row_start + n_col_start - j, n_rows);
            col_ind = j;
        }
        else{
            row_ind = j;
            col_ind = mod(n_col_start + n_row_start - j, n_cols);
        }
        diag(j) = kernel->getMatrixEntry(row_ind, col_ind);
    }
    return diag;
}

Mat getMatrix(const int n_row_start, const int n_col_start, const int n_rows, const int n_cols)
{
    Mat mat(n_rows, n_cols);
    // #pragma omp parallel for
    for (int j = 0; j < n_rows; ++j)
        // #pragma omp parallel for
        for (int k = 0; k < n_cols; ++k)
            mat(j, k) = kernel->getMatrixEntry(j + n_row_start, k + n_col_start);
    return mat;
}
// Like HODLR library it has the assemble feature that either constructs all 
// the matrix operators or just an memory efficient way 
 ~kernel_function(){
    // Destructor
}
};

#endif