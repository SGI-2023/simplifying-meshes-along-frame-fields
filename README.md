# Simplifyin Meshes along Frame Fields

### Description
Mesh simplification is useful for adapting high-detail art assets for real-time rendering. However, if we use the common decimation techniques based on iterated local edge collapse, the result may be hard to edit further (as artists tend to work directly with quad meshes that have good edge loop structure), or have irregular shading / deformation artifacts. Local edge collapse / flip techniques tend not to preserve most edge loops in the original mesh.
On the other hand, simplification techniques designed for quad meshes, such as polychord collapse, can preserve more of the edge loop structure of the original mesh, but are more complex (non-local operations on sequences of quads (polychords), or even global optimization/reconstruction as in quad remeshing).
Is there a simpler way that allows us to preserve these edge loops throughout simplification but keep ourselves to local operations (edge collapses, flips, regularization)? The original quad mesh might provide a starting point—it seems to define a frame field at most points on the mesh that we can use to perhaps guide our local operations.
The potential benefit of this is we can obtain nicely tessellated quad meshes from either quad meshes or even a not-nice triangle mesh, without having to resort to more complex techniques such as quad polychord collapse or quad remeshing).

### Prerequisites
[Concept of a frame field](https://cims.nyu.edu/gcl/papers/frame-fields-2014.pdf)

Students can use the geometry library of their choice that their are most comfortable with, or (TBD) be supplied with an initial one based on geometry-central that can also do quad simplification (as a comparison point).

### Questions & Steps 
1. Can we derive a frame field from the original quad mesh? (Start w/ Spot)
2. Can we come up with a reliable measure of how well a mesh’s edges align to frame field directions?
3. With a quad mesh, and a frame field with 2 directions, this seems easy. What about for a triangulated quad mesh where the original quad edges aren’t super obvious?
4. Try decimating a triangulated quad mesh, and measuring how well the result aligns with frame field directions.
5. Can we use all local collapse / flip / tangent smoothing to obtain results similar to quad simplification, and show that it gives a better result for the measure than triangle decimation?
