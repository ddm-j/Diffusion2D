#pragma once

// LinearSolvers.h

#ifndef LINEAR_SOLVERS_H  // Include guard to prevent multiple inclusion
#define LINEAR_SOLVERS_H

#include <Eigen/Sparse>
#include <Eigen/Dense>

// Successive Over Relaxation Solver
void successiveOverRelaxation(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::Ref<Eigen::VectorXd> x, double omega = 1.0, double tol = 0.001, int max_iter = 100000);

// Sparse SOR
void sparseSuccessiveOverRelaxation(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::VectorXd& b, Eigen::Ref<Eigen::VectorXd> x, double omega = 1.0, double tol = 0.001, int max_iter = 100000);


#endif // LINEAR_SOLVER_H
