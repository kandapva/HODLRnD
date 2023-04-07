#ifndef __LowRank_matrix__
#define __LowRank_matrix__

#include "myHeaders.hpp"
#include "kernel_function.hpp"

template <class Kernel>
class LowRankMat
{
    Mat L,R;
    bool is_mem_efficient = false;
    std::vector<size_t> row_id;
    std::vector<size_t> col_id;
    std::vector<size_t> row_basis;
    std::vector<size_t> col_basis;
    kernel_function<Kernel> *userkernel;
    ColPivHouseholderQR<Mat> *K;

public:
    LowRankMat(){
    }
    LowRankMat(kernel_function<Kernel> *&usr_, std::vector<size_t> &sources, std::vector<size_t> &targets,bool mem=false)
    {
        this->userkernel = usr_;
        row_id.assign(sources.begin(), sources.end());
        col_id.assign(targets.begin(), targets.end());
        if(SYS_SIZE < col_id.size() || mem){
            //std::cout << "MEM Eff" << std::endl;
            ACA_MEM_EFF(sources, targets);
            is_mem_efficient = true;
            //std::cout << "MEM done" << std::endl;
        }
        else{
            ACA_FAST(sources, targets);
        }
    }
    int rank(){
        return L.cols();
    }
    Vec operator * (Vec x)
    {
        Vec b;
        if (is_mem_efficient)
        {
            //std::cout << "(" << b.size() << ")" << std::endl;
            //std::cout << "(" << Ac.rows() << "," << Ac.cols() << ")" << std::endl;
            //std::cout << "(" << Ar.rows() << "," << Ar.cols() << ")" << std::endl;
            // Vec t0 = Ar * x;
            // std::cout << "(" << L.rows() << "," << L.cols() << ")" << std::endl;
            // std::cout << "(" << R.rows() << "," << R.cols() << ")" << std::endl;
            // Vec t1 = L.triangularView<Eigen::Lower>().solve(t0);
            // //std::cout << "(" << t1.size() << ")" << std::endl;
            // Vec t2 = R.triangularView<Eigen::Upper>().solve(t1);
            // //std::cout << "(" << t2.size() << ")" << std::endl;
            Mat Ac = userkernel->getMatrix(row_id, col_basis);
            Mat Ar = userkernel->getMatrix(row_basis, col_id);
            // Vec t2 = K->solve(Ar * x);
            // b = Ac * t2;
            Vec t0 = Ar * x;
            Vec t1 = L.triangularView<Eigen::Lower>().solve(t0);
            Vec t2 = R.triangularView<Eigen::Upper>().solve(t1);
            b = Ac * t2;
        }
        else
            b = L * (R.transpose() * x);
        return b;
    }
    void ACA_FAST(std::vector<size_t> &sources, std::vector<size_t> &targets)
    {
        double tolerance_or_rank = eps_ACA;
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
            row = userkernel->getRow(sources[row_ind.back()], targets);

            for (int i = 0; i < computed_rank; i++)
            {
                row = row - u[i](row_ind.back()) * v[i];
            }

            userkernel->maxAbsVector(row, remaining_col_ind, max, pivot);
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
                row = userkernel->getRow(sources[new_row_ind], targets);
                for (int i = 0; i < computed_rank; i++)
                {
                    row = row - u[i](row_ind.back()) * v[i];
                }

                userkernel->maxAbsVector(row, remaining_col_ind, max, pivot);
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
            col = userkernel->getCol(targets[col_ind.back()], sources);
            for (int i = 0; i < computed_rank; i++)
            {
                col = col - v[i](col_ind.back()) * u[i];
            }

            userkernel->maxAbsVector(col, remaining_row_ind, max, pivot);
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
                col = userkernel->getCol(targets[new_col_ind], sources);
                for (int i = 0; i < computed_rank; i++)
                {
                    col = col - v[i](col_ind.back()) * u[i];
                }

                userkernel->maxAbsVector(col, remaining_row_ind, max, pivot);
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
                R = userkernel->getMatrix(sources, targets).transpose();
                computed_rank = n_rows;
            }

            else
            {
                L = userkernel->getMatrix(sources, targets);
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

    void ACA_MEM_EFF(std::vector<size_t> &sources, std::vector<size_t> &targets)
    {
        is_mem_efficient = true;
        double tolerance_or_rank = eps_ACA;
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
            row = userkernel->getRow(sources[row_ind.back()], targets);

            for (int i = 0; i < computed_rank; i++)
            {
                row = row - u[i](row_ind.back()) * v[i];
            }

            userkernel->maxAbsVector(row, remaining_col_ind, max, pivot);
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
                row = userkernel->getRow(sources[new_row_ind], targets);
                for (int i = 0; i < computed_rank; i++)
                {
                    row = row - u[i](row_ind.back()) * v[i];
                }

                userkernel->maxAbsVector(row, remaining_col_ind, max, pivot);
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
            col = userkernel->getCol(targets[col_ind.back()], sources);
            for (int i = 0; i < computed_rank; i++)
            {
                col = col - v[i](col_ind.back()) * u[i];
            }

            userkernel->maxAbsVector(col, remaining_row_ind, max, pivot);
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
                col = userkernel->getCol(targets[new_col_ind], sources);
                for (int i = 0; i < computed_rank; i++)
                {
                    col = col - v[i](col_ind.back()) * u[i];
                }

                userkernel->maxAbsVector(col, remaining_row_ind, max, pivot);
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
                R = userkernel->getMatrix(sources, targets).transpose();
                computed_rank = n_rows;
            }

            else
            {
                L = userkernel->getMatrix(sources, targets);
                R = Mat::Identity(n_cols, n_cols);
                computed_rank = n_cols;
            }
            is_mem_efficient = false;
        }

        // This is when ACA has succeeded:
        else
        {
            for(int i=0; i<computed_rank; i++){
                row_basis.push_back(row_id[row_ind[i]]);
                
                col_basis.push_back(col_id[col_ind[i]]);
            }
            // row_basis.assign(row_ind.begin(),row_ind.end());
            // col_basis.assign(col_ind.begin(), col_ind.end());
            L = Mat::Zero(computed_rank, computed_rank);
            R = Mat::Zero(computed_rank, computed_rank);
            if (computed_rank > 0)
            {
                // Considered this QR mainly because of its comprise with speed and stability
                // Mat A = userkernel->getMatrix(row_basis,col_basis);
                // K = new ColPivHouseholderQR<Mat>(A); 
                for (int i = 0; i < computed_rank; i++)
                {
                    L(i, i) = 1.0;
                    if (i >= 1)
                    {
                        for (int j = 0; j <= i - 1; j++)
                        {
                            L(i, j) = u[j](row_ind[i]);
                        }
                    }
                }
                for (int i = 0; i < computed_rank; i++)
                {
                    R(i, i) = v[i](col_ind[i]);
                    if (i >= 1)
                    {
                        for (int j = 0; j <= i - 1; j++)
                        {
                            R(j, i) = v[j](col_ind[i]);
                            //R(j, i) = v[i](col_ind[j]);
                        }
                    }
                }
            }
        }
    }

    ~LowRankMat(){
    }
};

#endif