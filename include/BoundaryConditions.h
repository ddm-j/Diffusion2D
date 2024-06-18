#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <lua.hpp>

class BoundaryConditions {
private:
    lua_State* L = nullptr; // Lua Virtual Machine
    std::map<char, int> bcMap = {
        {'N', 0},
        {'S', 1},
        {'E', 2},
        {'W', 3}
    };

    std::vector<int> luaRefs;

    
    void initializeLua();
    void loadLuaFunction(const std::string funcName, const std::string UDF, const bool is_file, int& ref);
    double callLuaFunction(int& ref, double arg);

public:
    Eigen::Vector4i bcs; // Boundary Condition Types: 0 - Dirichlet, 1 - Neumann, 2 - UDF Dirichlet, 3 - UDF Neumann
    Eigen::Vector4d bcvs; // Boundary Condition Values
    std::vector<int> udfBounds;
    bool has_UDFs = false;


    // Class Methods
    /*BoundaryConditions(int N, double Tn, int S, double Ts, int E, double Te, int W, double Tw);
    BoundaryConditions(const Eigen::VectorXi& types, const Eigen::VectorXd& values);*/
    BoundaryConditions() : luaRefs(4, -1) {} // Initialize luaRefs to -1
    
    void addBC(char face, const int type, const double value);
    void addBC(char face, const int type, const std::string funcName, const std::string UDF);
    void updateBCs(const double arg);
};