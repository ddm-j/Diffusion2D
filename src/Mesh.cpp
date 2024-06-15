#include "Mesh.h"

// Mesh
Mesh::Mesh(double lx, double ly, int nx, int ny) {
    // Set Variables
    this->lx = lx;
    this->ly = ly;
    this->nx = nx;
    this->ny = ny;

    // Compute Mesh Spacing
    this->dx = lx / (nx - 1);
    this->dy = ly / (ny - 1);

    // Number of Unknowns
    this->Nunk = (nx - 2) * (ny - 2);
}

Mesh::Mesh(double l, int n) : Mesh(l, l, n, n) {}

void Mesh::initializeMatrix(const BoundaryConditions bounds, double f) {
    // Size the Matrix Accordingly
    int N = (this->nx - 2) * (this->ny - 2);
    this->A.resize(N, N);
    this->b.setConstant(N, -f);

    // Useful Variables
    int xe = nx - 1;
    int ye = ny - 1;
    double diag_coeff = -2.0 / std::pow(dx, 2) - 2.0 / std::pow(dy, 2);
    double ay = 1 / std::pow(dy, 2);
    double ax = 1 / std::pow(dx, 2);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(5 * N); // Approximate non-zero entries

    for (int j = 1; j < ye; j++) {
        for (int i = 1; i < xe; i++) {
            int index = (i - 1) * (nx - 2) + (j - 1); // Row index in the matrix

            // Reset Diagonal Value
            double diag = diag_coeff;

            // North Boundary
            if (i == 1) {
                if (bounds.bcs[0] == 0) { // Dirichlet
                    b(index) -= ay * bounds.bcvs[0];
                }
                else { // Neumann
                    diag += ay;
                    b(index) -= ay * bounds.bcvs[0] * dy;
                }
            }

            // South Boundary
            if (i == ny - 2) {
                if (bounds.bcs[1] == 0) { // Dirichlet
                    b(index) -= ay * bounds.bcvs[1];
                }
                else { // Neumann
                    diag += ay;
                    b(index) += ay * bounds.bcvs[1] * dy;
                }
            }

            // West Boundary
            if (j == 1) {
                if (bounds.bcs[3] == 0) { // Dirichlet
                    b(index) -= ax * bounds.bcvs[3];
                }
                else { // Neumann
                    diag += ax;
                    b(index) += ay * bounds.bcvs[3] * dx;
                }
            }

            // East Boundary
            if (j == nx - 2) {
                if (bounds.bcs[2] == 0) { // Dirichlet
                    b(index) -= ax * bounds.bcvs[2];
                }
                else { // Neumann
                    diag += ax;
                    b(index) -= ay * bounds.bcvs[2] * dx;
                }
            }

            // Diagonal element
            tripletList.push_back(Eigen::Triplet<double>(index, index, diag));
            // Other neighbors
            if (i > 1) tripletList.push_back(Eigen::Triplet<double>(index, index - (nx - 2), ay)); // North
            if (i < ny - 2) tripletList.push_back(Eigen::Triplet<double>(index, index + (nx - 2), ay)); // South
            if (j > 1) tripletList.push_back(Eigen::Triplet<double>(index, index - 1, ax)); // West
            if (j < nx - 2) tripletList.push_back(Eigen::Triplet<double>(index, index + 1, ax)); // East
        }

        // Set the Sparse Matrix
        A.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    // Build sparse matrix from triplets
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}
