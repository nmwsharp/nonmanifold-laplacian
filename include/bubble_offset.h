#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/simple_polygon_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class BubbleOffset {

public:
  BubbleOffset(EmbeddedGeometryInterface& geom);

  // Parameters
  double relativeScale = 0.01;
  double dialate = 0.0;
  double normalOffset = 0.0;
  bool useEdgeScaling = true;

  // Methods
  Vector3 queryPoint(const SurfacePoint& p);


  // Members
  EmbeddedGeometryInterface& geom;
  EdgeData<Vector3> edgeNormals;
};

std::unique_ptr<SimplePolygonMesh> subdivideRounded(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                                                    int subdivLevel, double scale, double dialate, double normalOffset);
