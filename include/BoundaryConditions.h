#pragma once

#include <Eigen/Dense>

class BoundaryConditions {
public:
    Eigen::Vector4i bcs; // Boundary Condition Types: 0 - Dirichlet, 1 - Neumann
    Eigen::Vector4d bcvs; // Boundary Condition Values

    // Class Methods
    BoundaryConditions(int N, double Tn, int S, double Ts, int E, double Te, int W, double Tw);
    BoundaryConditions(const Eigen::VectorXi& types, const Eigen::VectorXd& values);
};