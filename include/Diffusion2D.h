#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "LinearSolvers.h"
#include "VTRWriter.h"
#include "BoundaryConditions.h"
#include "Mesh.h"


class Diffusion2D {
public:
    Mesh mesh;
    BoundaryConditions bounds;
    VTRWriter* writer;

    // Solution Vector and Field
    Eigen::VectorXd solVec;
    Eigen::MatrixXd solField;

    // Mesh Setting
    int Nx;
    int Ny;

    // Constructor
    Diffusion2D(Mesh& mesh, const BoundaryConditions& bounds, double source, std::string name);

    // Class Methods
    void setDirichletBCs();

    void setNeumannBCs();

    void fillSolution();

    void solveSteady(double tol = 1e-6);

    double computeAutoTimestep(double sf = 0.8);

    void solveUnsteady(Eigen::Ref<Eigen::VectorXd> x, double finalTime, double dt, double write_freq, double tol = 1e-6);

    void solveUnsteady(double finalTime, double dt, double write_freq, double tol = 1e-6);
};
