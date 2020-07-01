# nonmanifold-laplacians

Robust Laplace operators for general (possibly nonmanifold) meshes. For details, see [A Laplacian for Nonmanifold Triangle Meshes](http://www.cs.cmu.edu/~kmcrane/Projects/NonmanifoldLaplace/NonmanifoldLaplace.pdf) by [Nicholas Sharp](http://nmwsharp.com) and [Keenan Crane](http://keenan.is/here) at SGP 2020.


This method efficiently generates a high-quality `V x V` Laplace matrix for any (possibly nonmanifold, with or without boundary) triangular 3D surface mesh. In particular, the resulting Laplacian will always satisfy the _maximum principle_, with all-positive edge weights between nodes. The method can also produce a similar Laplacian for a point cloud.

The core algorithm is implemented in [geometry-central](http://geometry-central.net). This is a simple application which loads a mesh or point cloud, builds our intrinsic tufted cover, and generates the resulting Laplace (and mass) matrix. 

As a command-line tool, the program will output matrices in a simple format which can be read in by other programs (including MATLAB, etc). It can also be run with the `--gui` flag, to create a window and visualize the relevant data---though be wary that this visualization is significantly more expensive and less robust than the main algorithm.

### Building and running

On unix-like machines, use:
```
git clone --recurse-submodules https://github.com/nmwsharp/nonmanifold-laplacian
cd nonmanifold-laplacian
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
./bin/tufted-idt /path/to/your/mesh.obj
```

The codebase also builds on Visual Studio 2017 & 2019, by using CMake to generate a Visual Studio solution file.

The input should be a mesh or point cloud; any inputs with no faces will be processed as point clouds. Use the `--gui` flag to load a 3D gui to inspect the results.

### Options

Use `--help` for defaults, etc.

| flag | purpose | 
| :------------- |:------------- |
| `--gui ` | Show the GUI| 
| `--mollifyFactor` | Amount of intrinsic mollification to apply, relative to the mesh length scale. Larger values will lend robustness to floating-point degeneracy, though very large values will distort geometry. Reasonable range is roughly 0 to 1e-3. Default: 1e-6 |
| `--nNeigh` | Number of nearest-neighbors to be used for point cloud Laplacian. The construction is not very sensitive to this parameter, it usually does not need to be tweaked. Default: 30 |
| `--outputPrefix` |  Prefix to prepend to all output file paths. Default: `tufted_` |
| `--writeLaplacian` | Write the resulting Laplace matrix. A sparse `VxV` matrix, holding the _weak_ Laplace matrix (that is, does not include mass matrix). Name: `laplacian.spmat` | |
| `--writeMass` | Write the resulting mass matrix. A sparse diagonal `VxV` matrix, holding lumped vertex areas. Name: `lumped_mass.spmat` | |


### Output formats

Sparse matrices are output as an ASCII file where each line one entry in the matrix, giving the row, column, and value. The row and column indices are **1-indexed** to make matlab happy. These files can be automatically loaded in matlab ([see here](https://www.mathworks.com/help/matlab/ref/spconvert.html)). Parsers in other environments should be straightforward.

### Direct Dependencies

- [geometry-central](http://geometry-central.net) for mesh data structures and 3D geometry
- [Polyscope](http://polyscope.run/) for 3D visualizations and rendering
- [jc_voronoi](https://github.com/JCash/voronoi) for generating 2D Delaunay triangulations for point clouds


(all are permissively licensed, and packaged with the repo)


### Citation

```bib
@article{Sharp:2020:LNT,
  author={Nicholas Sharp and Keenan Crane},
  title={{A Laplacian for Nonmanifold Triangle Meshes}},
  journal={Computer Graphics Forum (SGP)},
  volume={39},
  number={5},
  year={2020}
}
```
