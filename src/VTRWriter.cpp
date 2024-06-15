#include "VTRWriter.h"

VTRWriter::VTRWriter(int Nx, int Ny, double dx, double dy, const std::string& baseFilename)
    : Nx(Nx), Ny(Ny), dx(dx), dy(dy), baseFilename(baseFilename) {
    createMetadata();
    createSolutionDir();
}

void VTRWriter::createMetadata() {
    // Create XML File Structure

    // Header
    tinyxml2::XMLElement* vtkFile = doc.NewElement("VTKFile");
    vtkFile->SetAttribute("type", "RectilinearGrid");
    vtkFile->SetAttribute("version", "0.1");
    vtkFile->SetAttribute("byte_order", "LittleEndian");
    vtkFile->SetAttribute("header_type", "UInt64");
    doc.InsertFirstChild(vtkFile);

    // RectilinearGrid element
    tinyxml2::XMLElement* rectilinearGrid = doc.NewElement("RectilinearGrid");
    rectilinearGrid->SetAttribute("WholeExtent", ("0 " + std::to_string(Nx-1) + " 0 " + std::to_string(Ny-1) + " 0 0").c_str());
    vtkFile->InsertEndChild(rectilinearGrid);

    // Piece Element
    tinyxml2::XMLElement* piece = doc.NewElement("Piece");
    piece->SetAttribute("Extent", ("0 " + std::to_string(Nx-1) + " 0 " + std::to_string(Ny-1) + " 0 0").c_str());
    rectilinearGrid->InsertEndChild(piece);

    // Coordinates
    tinyxml2::XMLElement* coordinates = doc.NewElement("Coordinates");
    piece->InsertEndChild(coordinates);

    // X Coordinates
    tinyxml2::XMLElement* xCoords = doc.NewElement("DataArray");
    xCoords->SetAttribute("type", "Float32");
    xCoords->SetAttribute("Name", "X-Axis");
    xCoords->SetAttribute("format", "ascii");
    xCoords->SetAttribute("NumberOfComponents", "1");
    std::string xValues = "\n";
    for (int i = 0; i < Ny; ++i) {
        xValues += data_tabs + std::to_string(i * dy) + "\n";
    }
    xCoords->SetText(xValues.c_str());
    coordinates->InsertEndChild(xCoords);

    // Y Coordinates
    tinyxml2::XMLElement* yCoords = doc.NewElement("DataArray");
    yCoords->SetAttribute("type", "Float32");
    yCoords->SetAttribute("Name", "Y-Axis");
    yCoords->SetAttribute("format", "ascii");
    yCoords->SetAttribute("NumberOfComponents", "1");
    std::string yValues = "\n";
    for (int j = Nx-1; j >= 0; --j) {
        yValues += data_tabs + std::to_string(j * dx) + "\n";
    }
    yCoords->SetText(yValues.c_str());
    coordinates->InsertEndChild(yCoords);

    // Z Coordinates
    tinyxml2::XMLElement* zCoords = doc.NewElement("DataArray");
    zCoords->SetAttribute("type", "Float32");
    zCoords->SetAttribute("Name", "Z-Axis");
    zCoords->SetAttribute("format", "ascii");
    zCoords->SetAttribute("NumberOfComponents", "1");
    std::string zValues = "\n";
    zValues += data_tabs + "0.0" + "\n";
    zCoords->SetText(zValues.c_str());
    coordinates->InsertEndChild(zCoords);

    // Point Data
    pointData = doc.NewElement("PointData");
    pointData->SetAttribute("Scalars", "Temp");
    piece->InsertEndChild(pointData);

    // Data Array (Actual Data)
    dataArray = doc.NewElement("DataArray");
    dataArray->SetAttribute("type", "Float32");
    dataArray->SetAttribute("Name", "Temp");
    dataArray->SetAttribute("format", "ascii");
    dataArray->SetAttribute("NumberOfComponents", "1");
    pointData->InsertEndChild(dataArray);
}

void VTRWriter::writeVTRFile(const Eigen::MatrixXd& solField, int timestep) {
    // Create String of Scalar Values
    std::string scalarValues = "\n";
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            scalarValues += data_tabs + std::to_string(solField(i, j)) + "\n";
        }
    }

    // Set Data to XML Object
    dataArray->SetText(scalarValues.c_str());

    // Save out File
    std::string filename = generateFilename(timestep);
    tinyxml2::XMLPrinter printer;
    doc.Print(&printer); // Use the printer to format the XML
    std::ofstream fileStream(baseFilename + "/" + filename);
    fileStream << printer.CStr(); // Write the formatted XML to the file
    fileStream.close();
}

std::string VTRWriter::generateFilename(int timestep) {
    return baseFilename + "_" + std::to_string(timestep) + ".vtr";
}

void VTRWriter::createSolutionDir() {
    // Filesystem Namespace
    namespace fs = std::filesystem;
    std::regex pattern(baseFilename + "_.*\\.vtr");

    // Create Directory if Necessary
    if (!fs::exists(baseFilename)) {
        fs::create_directory(baseFilename);
    }
    else {
        // We should remove files that match the pattern
        for (const auto& entry : fs::directory_iterator(baseFilename)) {
            //std::cout << entry.path().filename().string() << std::endl;
            if (fs::is_regular_file(entry)) {
                std::cout << entry.path().filename().string() << std::endl;
                if (std::regex_match(entry.path().filename().string(), pattern)) {
                    fs::remove(entry);
                }
            }
        }
    }
}