# Robust Structure-based Shape Correspondence
Yanir Kleiman and Maks Ovsjanikov, CGF 2018

This code finds sturcture aware correspondences between regions on 3D shapes.
The input shapes can be provided in .off format.
The code allows matching each shape as triangular meshes or as point clouds.
Alternatively, the shapes can be provided as a matlab struct which contains 
a triangular meshes or a point cloud - for more details, see `ShapePairMapper.m`.

To run an example, run `run_example.m`. This file should generate images similar 
to the ones given in the `/results` folder. It is easy to switch the input format
of each shape from triangular mesh to point cloud (for more details see the code).

The main functions are detailed below (code taken from `run_example.m`):
```matlab
% Load shapes and generate shape graphs:
[M1, M2] = ShapePairMapper(filename1, filename2, ints_range, [], [], S1_opts, S2_opts);

% Compute the segment cor   respondences of the two shape graphs:
% use_val = 0 so the interval band id is not used for the matching of segments.
R = MatchShapes(M1, M2, 0); 

% Save results with visualization:
VisualizeMatching(R, save_name);
```

ShapePairMapper generates a graph structure for each input shape. The graphs are generated jointly.

MatchShapes matches the generated shape graphs and outputs the matching of each vertex.

VisualizeMatching displays the results and save both the image and the underlying data as a `.mat` file.

## Symmetry Breaking
The structure aware region correspondence outputs a symmetric matching between the regions;
a group of regions on the first shape (e.g. both hands) are mapped to a group of regions on the second shape.
We provide a simple heuristic to break these symmetries and output a one-to-one matching 
between the two shapes. To run this process run `run_breaksymmetry.m`.

## Prerequisites and credits
**This code requires MatlabBGL package which can be found here: http://dgleich.github.io/matlab-bgl/**

ANN MATLAB Wrapper package by Dahua Lin is included in a subfolder (`ann_mwrapper`) and provided as-is.

A geodesic computation package by Danil Kirsanov (https://code.google.com/archive/p/geodesic/) is included in a subfolder (`geodesics_matlab/geodesics`) and provided as-is.

Some of the code was written by Maks Ovsjanikov.

The repository includes four shapes out of the MPI FAUST shapes collection (http://faust.is.tue.mpg.de/).

--Yanir, 2018
