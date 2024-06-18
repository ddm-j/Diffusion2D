#pragma once

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "BoundaryConditions.h"

class Mesh {
private:
    BoundaryConditions* bounds = nullptr;
    double ax, ay;
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>* b_matrix;


    void updateBoundary(const double source, const int face_index);

public:
    // Mesh Dimensions
    double lx, ly; // Dimensions of the Domain
    int nx, ny; // Node Counts
    double dx, dy; // Mesh Spacing
    int Nunk;

    // Linear Alegbra Things
    Eigen::SparseMatrix<double, Eigen::RowMajor> A; // The A Matrix
    Eigen::VectorXd b; // The B Vector

    // Constructors
    // Full Mesh Spec.
    Mesh(double lx, double ly, int nx, int ny);
    // Square Mesh Spec.
    Mesh(double l, int n);


    // Class Methods
    void initializeMatrix(BoundaryConditions* bounds, double f);
    void updateBoundaries(const double source);
};
