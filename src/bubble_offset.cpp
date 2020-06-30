#include "bubble_offset.h"


BubbleOffset::BubbleOffset(EmbeddedGeometryInterface& geom_) : geom(geom_) {

  SurfaceMesh& mesh = geom.mesh;
  geom.requireVertexPositions();
  geom.requireFaceNormals();
  geom.requireEdgeLengths();
  geom.requireMeshLengthScale();

  // Compute edge normals
  edgeNormals = EdgeData<Vector3>(mesh, Vector3::zero());
  for (Edge e : mesh.edges()) {
    Vector3 n = Vector3::zero();
    n += geom.faceNormals[e.halfedge().face()];
    if (e.halfedge().twin().isInterior()) n += geom.faceNormals[e.halfedge().twin().face()];

    if (norm(n) < 1e-5) {
      Vector3 edgeDir =
          geom.vertexPositions[e.halfedge().twin().vertex()] - geom.vertexPositions[e.halfedge().vertex()];
      Vector3 faceN = geom.faceNormals[e.halfedge().face()];
      n = cross(edgeDir, faceN);
      n = unit(n);

      // n = Vector3::zero();
    } else {
      n = unit(n);
    }

    edgeNormals[e] = n;
  }
}


Vector3 BubbleOffset::queryPoint(const SurfacePoint& p) {

  // Position on the flat surface
  Vector3 pOrig = p.interpolate(geom.vertexPositions);

  // A point in some face
  SurfacePoint faceP = p.inSomeFace();
  Face origF = faceP.face;
  Vector3 bary = faceP.faceCoords;

  double scale = relativeScale * geom.meshLengthScale;

  // double minC = std::fmin(std::fmin(bary.x, bary.y), bary.z);
  // Vector3 offset = geom.faceNormals[origF] * scale * minC;

  // basis function coefs
  double u = bary.y;
  double v = bary.z;
  double w = bary.x;
  double phiI = (1 - u - v) * (1 - 2 * u - 2 * v);
  double phiJ = u * (2 * u - 1);
  double phiK = v * (2 * v - 1);
  double phiIJ = 4 * u * (1 - u - v);
  double phiJK = 4 * u * v;
  double phiKI = 4 * v * (1 - u - v);

  double lIJ = useEdgeScaling ? geom.edgeLengths[origF.halfedge().edge()] : 1.0;
  Vector3 offsetIJ = edgeNormals[origF.halfedge().edge()] * lIJ;

  double lJK = useEdgeScaling ? geom.edgeLengths[origF.halfedge().next().edge()] : 1.0;
  Vector3 offsetJK = edgeNormals[origF.halfedge().next().edge()] * lJK;

  double lKI = useEdgeScaling ? geom.edgeLengths[origF.halfedge().next().next().edge()] : 1.0;
  Vector3 offsetKI = edgeNormals[origF.halfedge().next().next().edge()] * lKI;

  Vector3 offset = offsetIJ * phiIJ + offsetJK * phiJK + offsetKI * phiKI;
  offset *= scale;
  offset += normalOffset * geom.faceNormals[origF];

  if (dialate > 0.0) {
    // pull towards center
    Vector3 initP = pOrig + offset;
    Vector3 faceCenter = Vector3::zero();
    for (Vertex v : origF.adjacentVertices()) {
      faceCenter += geom.vertexPositions[v] / 3;
    }
    return dialate * faceCenter + (1.0 - dialate) * initP;
  } else {
    return pOrig + offset;
  }
}


std::unique_ptr<SimplePolygonMesh> subdivideRounded(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                                                    int subdivLevel, double scale, double dialate,
                                                    double normalOffset) {

  // Make a copy
  geom.requireVertexPositions();
  geom.requireFaceNormals();
  geom.requireFaceAreas();
  geom.requireEdgeLengths();


  // == Good old-fashioned subdivision, preserving barycentric coords on to original triangle
  std::vector<std::array<size_t, 3>> faces;       // per-face
  std::vector<Vector3> vertexCoords;              // per-vert
  std::vector<std::array<Vector3, 3>> baryCoords; // per-face
  std::vector<Face> origFaceTri;                  // per-face

  { // Initialize
    for (Face f : mesh.faces()) {
      baryCoords.push_back({Vector3{1., 0., 0.}, Vector3{0., 1., 0.}, Vector3{0., 0., 1.}});
      origFaceTri.push_back(f);

      faces.push_back({vertexCoords.size(), vertexCoords.size() + 1, vertexCoords.size() + 2});

      // create distinct copies of the vertices
      vertexCoords.push_back(geom.vertexPositions[f.halfedge().vertex()]);
      vertexCoords.push_back(geom.vertexPositions[f.halfedge().next().vertex()]);
      vertexCoords.push_back(geom.vertexPositions[f.halfedge().next().next().vertex()]);
    }
  }

  // Iteratively subdivide
  for (int iSub = 0; iSub < subdivLevel; iSub++) {

    // New arrays
    std::vector<Vector3> vertexCoordsNew = vertexCoords; // per-vertex, copy existing
    // std::vector<Vector3> vertexOrigBaryCoordsNew = vertexOrigBaryCoords; // per-vertex, copy existing
    std::vector<std::array<size_t, 3>> facesNew;       // per-face, all new
    std::vector<Face> origFaceTriNew;                  // per-face, all new
    std::vector<std::array<Vector3, 3>> baryCoordsNew; // per-face, all new

    // Create new vertices along edges
    std::map<std::tuple<size_t, size_t>, size_t> edgeVert;
    std::vector<std::array<size_t, 3>> faceVert; // used if iSub == 0
    if (iSub == 0) {
      faceVert.resize(faces.size());
      // if (false) {
      // On the first iteration, separate the edges
      for (size_t iF = 0; iF < faces.size(); iF++) {
        for (size_t j = 0; j < 3; j++) {
          size_t vTail = faces[iF][j];
          size_t vTip = faces[iF][(j + 1) % 3];

          size_t newInd = vertexCoordsNew.size();
          Vector3 newPos = (vertexCoords[vTail] + vertexCoords[vTip]) / 2.;
          vertexCoordsNew.push_back(newPos);

          faceVert[iF][j] = newInd;
        }
      }
    } else {

      for (size_t iF = 0; iF < faces.size(); iF++) {
        for (size_t j = 0; j < 3; j++) {
          size_t vTail = faces[iF][j];
          size_t vTip = faces[iF][(j + 1) % 3];

          size_t newInd = vertexCoordsNew.size();
          Vector3 newPos = (vertexCoords[vTail] + vertexCoords[vTip]) / 2.;
          vertexCoordsNew.push_back(newPos);

          std::tuple<size_t, size_t> key{vTail, vTip};
          std::tuple<size_t, size_t> keyTwin{vTip, vTail};
          if (edgeVert.find(key) == edgeVert.end()) {
            edgeVert[key] = newInd;
            edgeVert[keyTwin] = newInd;
          }
        }
      }
    }


    // Create new faces
    for (size_t iF = 0; iF < faces.size(); iF++) {

      // Gather vertices
      size_t vA = faces[iF][0];
      size_t vB = faces[iF][1];
      size_t vC = faces[iF][2];

      size_t vAB, vBC, vCA;
      if (iSub == 0) {
        vAB = faceVert[iF][0];
        vBC = faceVert[iF][1];
        vCA = faceVert[iF][2];
      } else {

        if (edgeVert.find(std::tuple<size_t, size_t>{vA, vB}) == edgeVert.end())
          throw std::runtime_error("edge key " + std::to_string(vA) + " --- " + std::to_string(vB));
        if (edgeVert.find(std::tuple<size_t, size_t>{vB, vC}) == edgeVert.end())
          throw std::runtime_error("edge key " + std::to_string(vB) + " --- " + std::to_string(vC));
        if (edgeVert.find(std::tuple<size_t, size_t>{vC, vA}) == edgeVert.end())
          throw std::runtime_error("edge key " + std::to_string(vC) + " --- " + std::to_string(vA));


        vAB = edgeVert[std::tuple<size_t, size_t>{vA, vB}];
        vBC = edgeVert[std::tuple<size_t, size_t>{vB, vC}];
        vCA = edgeVert[std::tuple<size_t, size_t>{vC, vA}];
      }

      if (vAB == 0 || vBC == 0 || vCA == 0) {
        std::cout << "  face " << iF << " has verts " << vAB << " " << vBC << " " << vCA << std::endl;
      }

      Vector3 bA = baryCoords[iF][0];
      Vector3 bB = baryCoords[iF][1];
      Vector3 bC = baryCoords[iF][2];
      Vector3 bAB = 0.5 * (bA + bB);
      Vector3 bBC = 0.5 * (bB + bC);
      Vector3 bCA = 0.5 * (bC + bA);

      // Create new faces

      facesNew.push_back({vA, vAB, vCA});
      baryCoordsNew.push_back({bA, bAB, bCA});
      origFaceTriNew.push_back(origFaceTri[iF]);

      facesNew.push_back({vAB, vB, vBC});
      baryCoordsNew.push_back({bAB, bB, bBC});
      origFaceTriNew.push_back(origFaceTri[iF]);

      facesNew.push_back({vCA, vBC, vC});
      baryCoordsNew.push_back({bCA, bBC, bC});
      origFaceTriNew.push_back(origFaceTri[iF]);

      facesNew.push_back({vAB, vBC, vCA});
      baryCoordsNew.push_back({bAB, bBC, bCA});
      origFaceTriNew.push_back(origFaceTri[iF]);
    }

    // Swap in new arrays
    vertexCoords = vertexCoordsNew;
    faces = facesNew;
    origFaceTri = origFaceTriNew;
    baryCoords = baryCoordsNew;
  }

  BubbleOffset bubbleOffset(geom);
  bubbleOffset.relativeScale = scale;
  bubbleOffset.dialate = dialate;
  bubbleOffset.normalOffset = normalOffset;

  // Apply offsets
  // (this will process each vertex many times; that's fine)
  std::vector<Vector3> vertexCoordsOffset(vertexCoords.size());
  for (size_t iF = 0; iF < faces.size(); iF++) {
    for (size_t j = 0; j < 3; j++) {
      size_t vInd = faces[iF][j];
      Vector3 pOrig = vertexCoords[vInd];
      Face origF = origFaceTri[iF];
      Vector3 origBary = baryCoords[iF][j];

      SurfacePoint origP(origF, origBary);
      vertexCoordsOffset[vInd] = bubbleOffset.queryPoint(origP);
    }
  }


  // Store the output here
  std::unique_ptr<SimplePolygonMesh> outSoup(new SimplePolygonMesh({}, vertexCoordsOffset));
  for (size_t iF = 0; iF < faces.size(); iF++) {
    size_t vA = faces[iF][0];
    size_t vB = faces[iF][1];
    size_t vC = faces[iF][2];
    outSoup->polygons.push_back({vA, vB, vC});
  }

  return outSoup;
}
