#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/decimate.h>
#include <igl/principal_curvature.h>

#include "utils.h"

struct Mesh
{
public:
    MatrixXd V;
    MatrixXi F;
    MatrixXd PD1, PD2; // V * 3      Where V is the number of the Vertices
    VectorXd PV1, PV2; // V          Where V is the number of the Vertices
    std::vector<int> V2Fe; // adjacent faces to vertices
    VectorXd A; // frame field alignment cost vector

    igl::decimate_pre_collapse_callback custom_pre_collapse_callback;
    igl::decimate_post_collapse_callback custom_post_collapse_callback;
    igl::decimate_stopping_condition_callback stopping_condition_callback;

    Mesh(std::string filename)
    {
        igl::read_triangle_mesh(filename, V, F);
        get_frame_field();
        initialize_decimation_callbacks();
    };

    void initialize_decimation_callbacks();
    void decimate();
    void get_frame_field();
    void frame_field_alignment_data();
};