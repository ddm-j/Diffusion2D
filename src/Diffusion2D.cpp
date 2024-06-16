// 2DPoissonSolver.cpp : A Rectilinear Solver for 2D Diffusion Equations
//

#include "Diffusion2D.h"

// Pybind11 isht
#ifdef BUILD_PYTHON_BINDINGS
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#endif

Diffusion2D::Diffusion2D(Mesh& mesh, const BoundaryConditions& bounds, double source, std::string name) : mesh(mesh), bounds(bounds) {
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
    writer = new VTRWriter(Nx, Ny, this->mesh.dx, this->mesh.dy, name);
}

void Diffusion2D::setDirichletBCs() {
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

void Diffusion2D::setNeumannBCs() {
    // Mesh Spacing
    std::vector<double> d = { this->mesh.dy, this->mesh.dy, this->mesh.dx, this->mesh.dx };

    // Set the Values of the Neumann BCs
    for (int i = 0; i < 4; ++i) {
        if (this->bounds.bcs[i] == 1) {
            // Gradient Value
            double g = this->bounds.bcvs[i] * d[i];

            switch (i) {
            case 0: // North
                //std::cout << "Setting North Neumann" << std::endl;
                solField.block(0, 0, 1, Nx) = solField.block(1, 0, 1, Nx).array() + g;
                break;
            case 1: // South
                //std::cout << "Setting South Neumann" << std::endl;
                solField.block(Ny - 1, 0, 1, Nx) = solField.block(Ny - 2, 0, 1, Nx).array() - g;
                break;
            case 2: // East
                //std::cout << "Setting East Neumann" << std::endl;
                solField.block(0, Nx - 1, Ny, 1) = solField.block(0, Nx - 2, Ny, 1).array() + g;
                break;
            case 3: // West
                //std::cout << "Setting West Neumann" << std::endl;
                solField.block(0, 0, Ny, 1) = solField.block(0, 1, Ny, 1).array() - g;
                break;
            }
        }
    }
}

void Diffusion2D::fillSolution() {
    // Map solVec into the internal nodes of solField
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> internalNodes(solVec.data(), Ny - 2, Nx - 2);
    solField.block(1, 1, Ny - 2, Nx - 2) = internalNodes;
}

void Diffusion2D::solveSteady(double tol) {
    // Setup Solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
    cg.setTolerance(tol);
    cg.compute(mesh.A);

    // Solve
    solVec = cg.solve(mesh.b);
    fillSolution();
    std::cout << solField << std::endl;
    // Project Neumann BC
    setNeumannBCs();

    // Save Output
    writer->writeVTRFile(solField, 0);
}

double Diffusion2D::computeAutoTimestep(double sf) {
    double dt = 0.5 * sf * (mesh.dx * mesh.dx * mesh.dy * mesh.dy) / (mesh.dx * mesh.dx + mesh.dy * mesh.dy);
    return dt;
}

void Diffusion2D::solveUnsteady(Eigen::Ref<Eigen::VectorXd> x, double finalTime, double dt, double write_freq, double tol) {
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
        solVec = x + dt * (mesh.A * x - mesh.b);

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

void Diffusion2D::solveUnsteady(double finalTime, double dt, double write_freq, double tol) {
    // Initialize the solution vector to zero
    Eigen::VectorXd x(mesh.Nunk);
    x.setConstant(0.0);

    solveUnsteady(x, finalTime, dt, write_freq, tol);
}

int main()
{

    // Create Some BCs
    //                     N     S     E     W
    //BoundaryConditions bcs(0, 1, 0, 0, 0, 0, 0, 1);

    // Create a Mesh
    //Mesh mytestmesh(10.0, 100);
    
    // Define the Problem
    //Diffusion2D problem(mytestmesh, bcs, 0.0, "unsteady");

    // Get the steady solution
    // problem.solveSteady();

    // Solve the Unsteady Problem
    //double tFinal = 1000.0;
    //Eigen::VectorXd x(mytestmesh.Nunk);
    //x.setConstant(0.0);
    //problem.solveUnsteady(tFinal, 0.0, 0.1);

    // Testing Lua!
    // An example function
    float a = 4.0;
    float b = 2.0;

    std::string script_path = "myscript.lua";

    lua_State* L = luaL_newstate();
    luaL_dofile(L, script_path.c_str());
    lua_getglobal(L, "Pythagoras");
    if (lua_isfunction(L, -1)) {
        std::cout << "LUA File is a Function" << std::endl;
        // Push Function arguments onto the stack left to right (on the Lua function signature)
        lua_pushnumber(L, a);
        lua_pushnumber(L, b);

        // Lua Function call requirements
        constexpr int NUM_ARGS = 2;
        constexpr int NUM_RETURNS = 1;

        // Call the function
        lua_pcall(L, NUM_ARGS, NUM_RETURNS, 0);
        lua_Number lua_ret = lua_tonumber(L, -1);
        float c_sq = (float)lua_ret;

        std::cout << "Result of Lua Function: " << c_sq << std::endl;
    }
    lua_close(L);


}

#ifdef BUILD_PYTHON_BINDINGS
namespace py = pybind11;

PYBIND11_MODULE(Diffusion2D, m) {
    // Boundary Condition Class
    py::class_<BoundaryConditions>(m, "BoundaryConditions")
        .def(py::init<int, double, int, double, int, double, int, double>());

    // Mesh Class
    py::class_<Mesh>(m, "Mesh")
        .def(py::init<double, double, int, int>())
        .def_readonly("Nunk", &Mesh::Nunk)
        .def_readonly("nx", &Mesh::nx)
        .def_readonly("ny", &Mesh::ny);

    // Problem Class
    py::class_<Diffusion2D>(m, "Diffusion2D")
        .def(py::init<Mesh&, const BoundaryConditions&, double, std::string>())
        .def("solveSteady", &Diffusion2D::solveSteady,
            py::arg("tol") = 1e-6)
        .def("solveUnsteady",
            py::overload_cast<Eigen::Ref<Eigen::VectorXd>, double, double, double, double>(&Diffusion2D::solveUnsteady),
            py::arg("x"), py::arg("finalTime"), py::arg("dt"), py::arg("write_freq"), py::arg("tol") = 1e-6)
        .def_readonly("solField", &Diffusion2D::solField)
        .def_readonly("solVec", &Diffusion2D::solVec);
};

#endif
