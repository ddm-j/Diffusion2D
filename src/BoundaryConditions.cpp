#include <BoundaryConditions.h>

// Boundary Conditions
BoundaryConditions::BoundaryConditions(int N, double Tn, int S, double Ts, int E, double Te, int W, double Tw) {
    bcs << N, S, E, W;
    bcvs << Tn, Ts, Te, Tw;
}

BoundaryConditions::BoundaryConditions(const Eigen::VectorXi& types, const Eigen::VectorXd& values)
    : bcs(types), bcvs(values) {
    if (bcs.size() != bcvs.size()) {
        throw std::runtime_error("Boundary condition types and values must have the same size.");
    }
}