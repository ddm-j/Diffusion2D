// 2DPoissonSolver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "LinearSolvers.h"
#include "vtkWriter.h"

class BoundaryConditions {
public:
    Eigen::Vector4i bcs; // Boundary Condition Types: 0 - Dirichlet, 1 - Neumann
    Eigen::Vector4d bcvs; // Boundary Condition Values

    // Useful bools
    bool hasDirichlet = false;
    bool hasNeumann = false;

    BoundaryConditions(int N, double Tn, int S, double Ts, int E, double Te, int W, double Tw) {
        bcs << N, S, E, W;
        bcvs << Tn, Ts, Te, Tw;
    }

    BoundaryConditions(const Eigen::VectorXi& types, const Eigen::VectorXd& values)
        : bcs(types), bcvs(values) {
        if (bcs.size() != bcvs.size()) {
            throw std::runtime_error("Boundary condition types and values must have the same size.");
        }
    }
};

class Mesh {
public:
    // Mesh Dimensions
    double lx, ly; // Dimensions of the Domain
    int nx, ny; // Node Counts
    double dx, dy; // Mesh Spacing
    int Nunk;

    // Linear Alegbra Things
    Eigen::SparseMatrix<double, Eigen::RowMajor> A; // The A Matrix
    Eigen::VectorXd b; // The B Vector
    
    // Full Mesh Spec.
    Mesh(double lx, double ly, int nx, int ny) {
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

    // Square Mesh Spec.
    Mesh(double l, int n) : Mesh(l, l, n, n) {}

    // 

    // Sparse Matrix Creation
    void initializeMatrix(const BoundaryConditions bounds, double f) {
        // Size the Matrix Accordingly
        int N = (this->nx - 2) * (this->ny - 2);
        this->A.resize(N, N);
        this->b.setConstant(N, f);

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

    void printFullMatrix() {
        // Convert sparse matrix A to a dense matrix
        Eigen::MatrixXd denseA = Eigen::MatrixXd(A);

        // Print the dense matrix
        std::cout << "Full matrix representation:\n" << denseA << std::endl;
    }

};

class Poisson2D {
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

    Poisson2D(Mesh mesh, const BoundaryConditions& bounds, double source) : mesh(mesh), bounds(bounds) {
        // Set the Mesh Settings for Access
        this->Nx = this->mesh.nx;
        this->Ny = this->mesh.ny;

        // Initialize the A & b
        this->mesh.initializeMatrix(this->bounds, source);

        // Resize the Solution Vector/Field
        solVec.resize(this->mesh.A.rows());
        solField.resize(this->mesh.nx, this->mesh.ny);
        solField.setConstant(0.0);

        // Set the Dirichlet BC
        this->setDirichletBCs();

        // Initialize the VTR Writer
        writer = new VTRWriter(Nx, Ny, this->mesh.dx, this->mesh.dy, "poisson");
    }

    void setDirichletBCs() {
        // Boundary parameters [startRow, startCol, numRows, numCols] for North, South, East, West
        std::vector<std::vector<int>> params = {
            {0, 0, 1, Nx},       // North
            {Ny - 1, 0, 1, Nx},  // South
            {0, Nx - 1, Ny, 1},  // East
            {0, 0, Ny, 1}        // West
        };

        // Set the Values of the Dirichlet BCs
        for (int i = 0; i < 4; ++i) {
            if (this->bounds.bcs[i] == 0) {
                solField.block(params[i][0], params[i][1], params[i][2], params[i][3]).setConstant(this->bounds.bcvs[i]);
            }
        }
    }

    void setNeumannBCs() {
        // Mesh Spacing
        std::vector<double> d = { this->mesh.dy, this->mesh.dy, this->mesh.dx, this->mesh.dx };

        // Set the Values of the Neumann BCs
        for (int i = 0; i < 4; ++i) {
            if (this->bounds.bcs[i] == 1) {
                // Gradient Value
                double g = this->bounds.bcvs[i] * d[i];

                switch (i) {
                case 0: // North
                    solField.block(0, 0, 1, Nx) = solField.block(1, 0, 1, Nx).array() + g;
                
                case 1: // South
                    solField.block(Ny - 1, 0, 1, Nx) = solField.block(Ny - 2, 0, 1, Nx).array() - g;

                case 2: // East
                    solField.block(0, Nx - 1, Ny, 1) = solField.block(0, Nx - 2, Ny, 1).array() + g;

                case 3: // West
                    solField.block(0, 0, Ny, 1) = solField.block(0, 1, Ny, 1).array() - g;
                }
            }
        }
    }

    void fillSolution() {
        // Map solVec into the internal nodes of solField
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> internalNodes(solVec.data(), Ny - 2, Nx - 2);
        solField.block(1, 1, Ny - 2, Nx - 2) = internalNodes;
    }

    void solveSteady(double tol=1e-6) {
        // Setup Solver
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
        cg.setTolerance(tol);
        cg.compute(mesh.A);

        // Solve
        solVec = cg.solve(mesh.b);
        fillSolution();

        // Project Neumann BC
        setNeumannBCs();

        // Save Output
        writer->writeVTRFile(solField, 0);
    }

    double computeAutoTimestep(double sf = 0.8) {
        double dt = 0.5 * sf * (mesh.dx * mesh.dx * mesh.dy * mesh.dy) / (mesh.dx * mesh.dx + mesh.dy * mesh.dy);
        return dt;
    }

    void solveUnsteady(Eigen::VectorXd x, double finalTime, double dt, double write_freq, double tol=1e-6) {
        // Adjust timestep if necessary
        dt = (dt == 0) ? computeAutoTimestep() : dt;

        // Convergence
        Eigen::VectorXd r;
        r.setConstant(0.0);
        double bNorm = mesh.b.norm();

        // Time stepping
        int out_cnt = 0;
        double lastWrite = 0.0;
        for (double t = 0.0; t <= finalTime; t += dt) {
            std::cout << "Timestepping: " << t << std::endl;

            // Compute New Solution
            solVec = x + dt * (mesh.A * x + mesh.b);

            // Check Convergence
            r = solVec - x;
            double rNorm = r.norm();
            if (rNorm / bNorm < tol) {
                std::cout << "Problem Converged." << std::endl;
                break;
            }
            else {
                std::cout << "Current Residual: " << rNorm << std::endl;
                x = solVec;
            }

            // Output 
            if (t - lastWrite >= write_freq || t == 0) {
                std::cout << "Writing output." << std::endl;

                // Fill and Project
                fillSolution();
                setNeumannBCs();

                // Write out
                writer->writeVTRFile(solField, out_cnt);

                // Update Counters
                lastWrite = t;
                out_cnt += 1;
            }
        }

    }
};

int main()
{

    // Create Some BCs
    //                     N     S     E     W
    BoundaryConditions bcs(1, 0, 1, 0, 0, 0, 0, 0);

    // Create a Mesh
    Mesh mytestmesh(10.0, 10);
    
    // Define the Problem
    Poisson2D problem(mytestmesh, bcs, -1.0);

    // Get the steady solution
    // std::cout << "Eigen's CG Solver: " << std::endl;
    problem.solveSteady();

    // Solve the Unsteady Problem
    //double tFinal = 1000.0;
    //Eigen::VectorXd x(mytestmesh.Nunk);
    //x.setConstant(0.0);
    //problem.solveUnsteady(x, tFinal, 0.0, 0.1);

    // Test the SOR Solver
    //Eigen::VectorXd x = Eigen::VectorXd::Zero(problem.mesh.A.rows());
    //sparseSuccessiveOverRelaxation(problem.mesh.A, problem.mesh.b, x, 1.7, 1e-6);
    //std::cout << "My SOR Solver: ";
    //std::cout << x.transpose() << std::endl;

    //// Create a 5x5 sparse matrix
    //Eigen::SparseMatrix<double, Eigen::RowMajor> mat(5, 5);

    //// Fill the matrix with some values for demonstration
    //mat.insert(0, 1) = 10;  // mat(0,1) = 10
    //mat.insert(1, 2) = 20;  // mat(1,2) = 20
    //mat.insert(2, 3) = 30;  // mat(2,3) = 30
    //mat.insert(3, 4) = 40;  // mat(3,4) = 40
    //mat.insert(4, 0) = 50;  // mat(4,0) = 50

    //// Finalize matrix assembly
    ////mat.makeCompressed();

    //Eigen::MatrixXd denseMat = Eigen::MatrixXd(mat);
    //std::cout << denseMat << std::endl;

    //// Iterate over all elements in the matrix
    //int cnt = 0;
    //for (int k = 0; k < mat.outerSize(); ++k) {
    //    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, k); it; ++it) {
    //        // it.row()       -- row index
    //        // it.col()       -- column index
    //        // it.value()     -- value of the element
    //        std::cout << "Element (" << it.row() << "," << it.col() << ") = " << it.value() << std::endl;
    //        ++cnt;
    //    }
    //}
    //std::cout << "Total iterations: " << cnt << ". Out of " << mat.rows() * mat.cols() << " possible iterations." << std::endl;

    //return 0;

}

