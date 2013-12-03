load vertex.mat;

find(isnan(vertex))
find(isinf(vertex))

build_params.algorithm='kdtree_single';
build_params.trees=1;
[index, parameters, speedup] = flann_build_index(vertex, build_params);
