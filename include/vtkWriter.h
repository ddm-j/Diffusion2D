#pragma once

#include "tinyxml2.h"
#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <fstream>

class VTRWriter {
public:
    VTRWriter(int Nx, int Ny, double dx, double dy, const std::string& baseFilename);
    void writeVTRFile(const Eigen::MatrixXd& solField, int timestep);

private:
    int Nx, Ny;
    double dx, dy;
    std::string baseFilename;
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLElement* pointData;
    tinyxml2::XMLElement* dataArray;

    void createMetadata();

    std::string generateFilename(int timestep);
};