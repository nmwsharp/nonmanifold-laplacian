#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/halfedge_factories.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

// Polyscope visualization handle, to quickly add data to the surface
bool withGUI = true;
polyscope::SurfaceMesh* psMesh;

// Parameters
float mollifyFactor = 0.;


template <typename T>
void saveMatrix(std::string filename, SparseMatrix<T>& matrix) {

  // WARNING: this follows matlab convention and thus is 1-indexed

  std::cout << "Writing sparse matrix to: " << filename << std::endl;

  std::ofstream outFile(filename);
  if (!outFile) {
    throw std::runtime_error("failed to open output file " + filename);
  }

  // Write a comment on the first line giving the dimensions
  outFile << "# sparse " << matrix.rows() << " " << matrix.cols() << std::endl;

  outFile << std::setprecision(16);

  for (int k = 0; k < matrix.outerSize(); ++k) {
    for (typename SparseMatrix<T>::InnerIterator it(matrix, k); it; ++it) {
      T val = it.value();
      size_t iRow = it.row();
      size_t iCol = it.col();

      outFile << (iRow + 1) << " " << (iCol + 1) << " " << val << std::endl;
    }
  }

  outFile.close(); 
}

template <typename T>
void saveMatrix(std::string filename, DenseMatrix<T>& matrix) {

  std::cout << "Writing dense matrix to: " << filename << std::endl;

  std::ofstream outFile(filename);
  if (!outFile) {
    throw std::runtime_error("failed to open output file " + filename);
  }

  // Write a comment on the first line giving the dimensions
  outFile << "# dense " << matrix.rows() << " " << matrix.cols() << std::endl;

  outFile << std::setprecision(16);

  for (size_t iRow = 0; iRow < (size_t)matrix.rows(); iRow++) {
    for (size_t iCol = 0; iCol < (size_t)matrix.cols(); iCol++) {
      T val = matrix(iRow, iCol);
      outFile << val;
      if (iCol + 1 != (size_t)matrix.cols()) {
        outFile << " ";
      }
    }
    outFile << std::endl;
  }

  outFile.close();
}

void myCallback() {

  ImGui::PushItemWidth(100);

  ImGui::TextUnformatted("Intrinsic triangulation:");
  //ImGui::Text("  nVertices = %lu  nFaces = %lu", signpostTri->mesh.nVertices(), signpostTri->mesh.nFaces());

  if (ImGui::TreeNode("Output")) {

    ImGui::TreePop();
  }

  ImGui::PopItemWidth();
}

int main(int argc, char** argv) {

  // Configure the argument parser
  // clang-format off
  args::ArgumentParser parser("Demo for ");
  args::HelpFlag help(parser, "help", "Display this help message", {'h', "help"});
  args::Positional<std::string> inputFilename(parser, "mesh", "A surface mesh file (see geometry-central for valid formats)");

  args::Group algorithmOptions(parser, "algorithm options");
  args::ValueFlag<double> mollifyFactorArg(algorithmOptions, "mollifyFactor", "Amount of intrinsic mollification to perform, which gives robustness to degenerate triangles. Defined relative to the mean edge length. Default: 1e-6", {"mollifyFactor"}, 0);

  args::Group output(parser, "ouput");
  args::Flag noGUI(output, "noGUI", "exit after processing and do not open the GUI", {"noGUI"});
  args::ValueFlag<std::string> outputPrefixArg(output, "outputPrefix", "Prefix to prepend to output file paths. Default: tufted_", {"outputPrefix"}, "tufted_");
  args::Flag writeLaplace(output, "writeLaplace", "Write out the resulting (weak) Laplacian as a sparse matrix. name: 'laplacian.dmat'", {"writeLaplace"});
  args::Flag writeMass(output, "writeMass", "Write out the resulting diagonal lumped mass matrix sparse matrix. name: 'lumped_mass.dmat'", {"writeMass"});
  // clang-format on

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cout << parser;
    return EXIT_FAILURE;
  }

  // Set options
  withGUI = !noGUI;
  mollifyFactor = args::get(mollifyFactorArg);

  // Load mesh
  SimplePolygonMesh inputMesh(args::get(inputFilename));

  inputMesh.stripFacesWithDuplicateVertices();
  inputMesh.stripUnusedVertices();
  inputMesh.triangulate();
    

  std::tie(mesh, geometry) = makeGeneralHalfedgeAndGeometry(inputMesh.polygons, inputMesh.vertexCoordinates);

  if (withGUI) {

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;
 
    // The mesh
    polyscope::registerSurfaceMesh("input mesh", inputMesh.vertexCoordinates, inputMesh.polygons);

    /*
    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(args::get(inputFilename)),
                                            geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // Set vertex tangent spaces
    geometry->requireVertexTangentBasis();
    VertexData<Vector3> vBasisX(*mesh);
    for (Vertex v : mesh->vertices()) {
      vBasisX[v] = geometry->vertexTangentBasis[v][0];
    }
    polyscope::getSurfaceMesh()->setVertexTangentBasisX(vBasisX);

    // Set face tangent spaces
    geometry->requireFaceTangentBasis();
    FaceData<Vector3> fBasisX(*mesh);
    for (Face f : mesh->faces()) {
      fBasisX[f] = geometry->faceTangentBasis[f][0];
    }
    polyscope::getSurfaceMesh()->setFaceTangentBasisX(fBasisX);
    */
  }


  // Perform any operations requested

  // Give control to the polyscope gui
  if (withGUI) {
    polyscope::show();
  }

  return EXIT_SUCCESS;
}
