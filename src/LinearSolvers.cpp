#include <iostream>

#include "LinearSolvers.h"


void successiveOverRelaxation(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::Ref<Eigen::VectorXd> x, double omega, double tol, int max_iter) {
    // Jacobi Iterative System Solver
    int n = A.rows();
    double e = tol + 1.0;
    double new_x = 0.0;

    // Iteration Loop
    for (int it = 0; it < max_iter; ++it) {
        e = 0.0;

        // Row-wise loop
        for (int i = 0; i < n; ++i) {
            double sigma = 0.0;

            // Col-wise loop
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sigma += A(i, j) * x[j];
                }
            }
            // Compute Error and Update
            new_x = (b[i] - sigma) / A(i, i);  // calculate the new x[i] based on the other current x[j] values
            new_x = x[i] + omega * (new_x - x[i]);  // adjust x[i] towards the new value

            e = std::max(e, std::abs(new_x - x[i]));  // track the maximum change for convergence check
            x[i] = new_x;
        }

        // Convergence Check
        if (e < tol) {
            std::cout << "Converged after " << it + 1 << " iterations." << std::endl;;
            return;
        }
    }

    std::cout << "Solver failed to converged after " << max_iter << " iterations." << std::endl;
}

void sparseSuccessiveOverRelaxation(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::VectorXd& b, Eigen::Ref<Eigen::VectorXd> x, double omega, double tol, int max_iter) {
    double bNorm = b.norm();        // Norm of the right-hand side vector


    // Iteration Loop
    for (int iteration = 0; iteration < max_iter; ++iteration) {
        double r_i = 0.0; // Row-wise Residual
        double r_sum = 0.0; // Running sum of squared residuals

        // Row-wise loop
        for (int i = 0; i < A.outerSize(); ++i) {
            double sigma = 0.0;
            double diag = 0.0;

            // Col-wise loop
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator iit(A, i); iit; ++iit) {
                int j = iit.col();
                double aij = iit.value();

                if (i == j) {
                    diag = aij;
                }
                else {
                    sigma += aij * x[j];
                }

            }
            // Update Solution
            x[i] = (1 - omega) * x[i] + omega * (b[i] - sigma) / diag;

            // Compute Error
            r_i = b[i] - (sigma + diag * x[i]);
            r_sum += r_i * r_i;
        }

        // Convergence Check
        double rNorm = std::sqrt(r_sum);
        if (rNorm/bNorm < tol) {
            std::cout << "Converged after " << iteration + 1 << " iterations." << std::endl;;
            return;
        }
    }

    std::cout << "Solver failed to converged after " << max_iter << " iterations." << std::endl;
}
