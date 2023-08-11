categorized all 10K meshes in Thingi into the following:

-Too Nice Mesh Paths (728 files): These are triangle meshes that have No skinny faces. Also all of them are manifold, oriented, connected and < 20 MBs
-valid Mesh Paths: These are triangle meshes that have some skinny faces (some faces: C + B < A + (1e-3 * AvgEdgeLength)). Also all of them are manifold, oriented, connected and < 20 MBs
-Invalid Mesh Paths (5053 files): These are meshes that aren't easy to work with. i.e. Non-manifold, disconnected, non-oriented, not triangle or couldn't be loaded at all.
-too Large Mesh Paths (457 files): These files are > 20 MBs so I didn't parse them intentionally to save time. (We can do that later if needed)

These files ids are stored in a .txt file for each category, So one can easily take the id from there and get it from https://ten-thousand-models.appspot.com/
