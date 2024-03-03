#pragma once

#include "mesh.h"
#include "qslim_aligned.h"
#include "igl/qslim.h"

struct Window {
public:
  Mesh mesh;
  igl::opengl::glfw::Viewer viewer;

  Window(Mesh mesh) : mesh(mesh) {
    viewer.data().set_mesh(mesh.V, mesh.F);
    initialize_callbacks();
  };

  void initialize_callbacks();
  void visualizeFrameFields(RowVector3d maxCurvatureColor,
                            RowVector3d minCurvatureColor);
  void visualizeMetricOnEdges();
};
