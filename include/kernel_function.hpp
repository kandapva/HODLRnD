#ifndef __kernel_function__
#define __kernel_function__

#include "myHeaders.hpp"
#include "points_dt.hpp"


template <class Operator>
class kernel_function
{
    Operator *kernel;
public :
    kernel_function(Operator*& A){
        this->kernel = A;
    }
    // Returns individual entries of the matrix:
    // this function can be complemented through a user kernel class
    // try to use the same function, expected the kernel evaluation at ith and jth gridPoints
    // Essentially it queries into K_{ij}
    Vec getRow(int j, std::vector<size_t> &targets)
    {
        int n_cols = targets.size();
        Vec row(n_cols);
        // #pragma omp parallel for
        for (int k = 0; k < n_cols; k++){
            row(k) = kernel->getMatrixEntry(j, targets[k]);
        }
            
        return row;
    }

    Vec getCol(int k, std::vector<size_t> &sources)
    {
        int n_rows = sources.size();
        Vec col(n_rows);
        // #pragma omp parallel for
        for (int j = 0; j < n_rows; ++j)
            col(j) = kernel->getMatrixEntry(sources[j], k);
    return col;
    }

// const int n_row_start, const int n_col_start, const int n_rows, const int n_cols
Mat getMatrix(std::vector<size_t> &sources, std::vector<size_t> &targets){
    int n_rows = sources.size();
    int n_cols = targets.size();
    Mat mat(n_rows, n_cols);
    // #pragma omp parallel for
    for (int j = 0; j < n_rows; ++j)
        // #pragma omp parallel for
        for (int k = 0; k < n_cols; ++k)
            mat(j, k) = kernel->getMatrixEntry(sources[j], targets[k]);
    return mat;
}

void maxAbsVector(const Vec &v, const std::set<size_t> &allowed_indices, dtype &max, size_t &index){
    std::set<size_t>::iterator it;
    index = *allowed_indices.begin();
    max = v(index);

    for (it = allowed_indices.begin(); it != allowed_indices.end(); it++)
    {
        if (fabs(v(*it)) > fabs(max))
        {
            index = *it;
            max = v(index);
        }
    }
}

void ACA_FAST(Mat &L, Mat &R, 
                double tolerance_or_rank, std::vector<size_t>& sources, std::vector<size_t>& targets)
{
    int n_rows = sources.size();
    int n_cols = targets.size();

    // Indices which have been used:
    std::vector<size_t> row_ind;
    std::vector<size_t> col_ind;

    // Indices that are remaining:
    std::set<size_t> remaining_row_ind;
    std::set<size_t> remaining_col_ind;

    // Bases:
    std::vector<Vec> u;
    std::vector<Vec> v;

    for (int k = 0; k < n_rows; k++)
    {
        remaining_row_ind.insert(k);
    }

    for (int k = 0; k < n_cols; k++)
    {
        remaining_col_ind.insert(k);
    }

    dtype max, gamma;

    // Initialize the matrix norm and the the first row index
    dtype_base matrix_norm = 0;
    row_ind.push_back(0);
    remaining_row_ind.erase(0);

    // Stores the pivot entry of the considered row / col:
    size_t pivot;

    int target_rank = 0;
    // This would get updated:
    int computed_rank = 0;
    Vec row, col;

    double tolerance = 0;
    // These quantities in finding the stopping criteria:
    dtype_base row_squared_norm, row_norm, col_squared_norm, col_norm;

    if (tolerance_or_rank < 1)
        tolerance = tolerance_or_rank;
    else
        target_rank = tolerance_or_rank;

    // So these would be particularly useful for poorly conditioned matrices:
    int max_tries = 10;
    int count;

    // Repeat till the desired tolerance / rank is obtained
    do
    {
        // Generation of the row
        // Row of the residuum and the pivot column
        // By calling row_ind.back(), we are getting the last pushed number
        row = this->getRow(sources[row_ind.back()], targets);

        for (int i = 0; i < computed_rank; i++)
        {
            row = row - u[i](row_ind.back()) * v[i];
        }

        this->maxAbsVector(row, remaining_col_ind, max, pivot);
        count = 0;

        // Alternating upon each call:
        bool eval_at_end = false;
        // Toggling randomness
        bool use_randomization = true;

        // This while loop is needed if in the middle of the algorithm the
        // row happens to be exactly the linear combination of the previous rows
        // upto some tolerance. i.e. prevents from ACA throwing false positives
        while (fabs(max) < tolerance &&
               count < max_tries &&
               remaining_col_ind.size() > 0 &&
               remaining_row_ind.size() > 0)
        {
            row_ind.pop_back();
            int new_row_ind;

            // When rank < 3, we will just choose entries from the ends of the matrix:
            if (computed_rank < 3)
            {
                if (eval_at_end == true)
                {
                    new_row_ind = *--remaining_row_ind.end();
                }

                else
                {
                    new_row_ind = *remaining_row_ind.begin();
                }

                eval_at_end = !(eval_at_end);
            }

            // However, when we have rank >=3, we will choose the entries such that
            // the newly picked entry is at the mid-point of the already chosen ones:
            else
            {
                if (use_randomization == true)
                {
                    std::set<size_t>::const_iterator it(remaining_row_ind.begin());
                    std::advance(it, rand() % remaining_row_ind.size());
                    new_row_ind = *it;
                }

                else
                {
                    std::vector<size_t> row_ind_sort(row_ind);
                    std::sort(row_ind_sort.begin(), row_ind_sort.end());
                    std::vector<size_t> row_ind_diff(row_ind_sort.size() - 1);

                    int max = 0;
                    int idx = 0;

                    for (size_t i = 0; i < row_ind_sort.size() - 1; i++)
                    {
                        row_ind_diff[i] = row_ind_sort[i + 1] - row_ind_sort[i];
                        if (row_ind_diff[i] > (unsigned) max)
                        {
                            idx = i;
                            max = row_ind_diff[i];
                        }
                    }

                    new_row_ind = row_ind_sort[idx] + max / 2;
                }

                use_randomization = !(use_randomization);
            }

            row_ind.push_back(new_row_ind);
            remaining_row_ind.erase(new_row_ind);
            // Generation of the row
            // Row of the residuum and the pivot column
            row = this->getRow(sources[new_row_ind], targets);
            for (int i = 0; i < computed_rank; i++)
            {
                row = row - u[i](row_ind.back()) * v[i];
            }

            this->maxAbsVector(row, remaining_col_ind, max, pivot);
            count++;
        }

        // In case it failed to resolve in the previous step,
        // we break out of the dowhile loop:
        if (count == max_tries ||
            remaining_col_ind.size() == 0 ||
            remaining_row_ind.size() == 0)
        {
            break;
        }

        // Resetting count back to zero for columns:
        count = 0;

        col_ind.push_back(pivot);
        remaining_col_ind.erase(pivot);
        // Normalizing constant
        gamma = dtype_base(1.0) / max;

        // Generation of the column
        // Column of the residuum and the pivot row
        col = this->getCol(targets[col_ind.back()], sources);
        for (int i = 0; i < computed_rank; i++)
        {
            col = col - v[i](col_ind.back()) * u[i];
        }

        this->maxAbsVector(col, remaining_row_ind, max, pivot);
        // Repeating the same randomization we carried out for the rows, now for the columns:
        while (fabs(max) < tolerance &&
               count < max_tries &&
               remaining_col_ind.size() > 0 &&
               remaining_row_ind.size() > 0)
        {
            col_ind.pop_back();

            int new_col_ind;

            if (col_ind.size() < 3)
            {
                if (eval_at_end)
                {
                    new_col_ind = *remaining_col_ind.end();
                }

                else
                {
                    new_col_ind = *remaining_col_ind.begin();
                }

                eval_at_end = !eval_at_end;
            }

            else
            {
                if (use_randomization == true)
                {
                    std::set<size_t>::const_iterator it(remaining_col_ind.begin());
                    std::advance(it, rand() % remaining_col_ind.size());
                    new_col_ind = *it;
                }

                else
                {
                    std::vector<size_t> col_ind_sort(col_ind);
                    std::sort(col_ind_sort.begin(), col_ind_sort.end());
                    std::vector<size_t> col_ind_diff(col_ind_sort.size() - 1);

                    int max = 0;
                    int idx = 0;

                    for (size_t i = 0; i < col_ind_sort.size() - 1; i++)
                    {
                        col_ind_diff[i] = col_ind_sort[i + 1] - col_ind_sort[i];
                        if (col_ind_diff[i] > (unsigned) max)
                        {
                            idx = i;
                            max = col_ind_diff[i];
                        }
                    }

                    new_col_ind = col_ind_sort[idx] + max / 2;
                }

                use_randomization = !(use_randomization);
            }

            col_ind.push_back(new_col_ind);
            remaining_col_ind.erase(new_col_ind);

            // Generation of the column
            // Column of the residuum and the pivot row:
            col = this->getCol(targets[new_col_ind], sources);
            for (int i = 0; i < computed_rank; i++)
            {
                col = col - v[i](col_ind.back()) * u[i];
            }

            this->maxAbsVector(col, remaining_row_ind, max, pivot);
            count++;
        }

        row_ind.push_back(pivot);
        remaining_row_ind.erase(pivot);

        // New vectors
        u.push_back(gamma * col);
        v.push_back(row);

        // New approximation of matrix norm
        row_squared_norm = row.squaredNorm();
        row_norm = sqrt(row_squared_norm);

        col_squared_norm = col.squaredNorm();
        col_norm = sqrt(col_squared_norm);

        // Updating the matrix norm:
        matrix_norm += std::abs(gamma * gamma * row_squared_norm * col_squared_norm);

        for (int j = 0; j < computed_rank; j++)
        {
            matrix_norm += 2.0 * std::abs(u[j].dot(u.back())) * std::abs(v[j].dot(v.back()));
        }

        computed_rank++;
    } while (((tolerance_or_rank < 1) ? computed_rank * (n_rows + n_cols) * row_norm * col_norm >
                                            fabs(max) * tolerance * matrix_norm
                                      : computed_rank < target_rank) &&
             computed_rank < fmin(n_rows, n_cols));

    // If the computed_rank is >= to full-rank
    // then return the trivial full-rank decomposition
    if (computed_rank >= fmin(n_rows, n_cols) - 1)
    {
        if (n_rows < n_cols)
        {
            L = Mat::Identity(n_rows, n_rows);
            R = this->getMatrix(sources,targets).transpose();
            computed_rank = n_rows;
        }

        else
        {
            L = this->getMatrix(sources, targets);
            R = Mat::Identity(n_cols, n_cols);
            computed_rank = n_cols;
        }
    }

    // This is when ACA has succeeded:
    else
    {
        L = Mat(n_rows, computed_rank);
        R = Mat(n_cols, computed_rank);

        for (int j = 0; j < computed_rank; j++)
        {
            L.col(j) = u[j];
            R.col(j) = v[j];
        }
    }
}
// Like HODLR library it has the assemble feature that either constructs all 
// the matrix operators or just an memory efficient way
~kernel_function()
{
    // Destructor
}
};

#endif