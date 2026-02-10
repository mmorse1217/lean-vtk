# Lean-VTK ![build and test](https://github.com/mmorse1217/lean-vtk/workflows/build%20and%20test/badge.svg)

A barebones VTK writer for C++. It is implemented as a single header and source file for easy integration with other projects and only requires the C++ standard library as a dependency.

### Example usage
Writing a single triangle to a vtk file:
```cpp
    vector<double> points = {
             1.,  1., -1.,
             1., -1., 1.,
            -1., -1., 0.
        };
    vector<size_t> elements = { 0, 1, 2 };
    vector<double> scalar_field = { 0., 1., 2.  };
    vector<double> vector_field = points;;

    const size_t dim = 3;
    const size_t cell_size = 3;
    std::string filename = "single_tri.vtu";
    VTUWriter writer;
    
    writer.add_scalar_field("scalar_field", scalar_field);
    writer.add_vector_field("vector_field", vector_field, dim);

    writer.write_surface_mesh(filename, dim, cell_size, points, elements);
```
Other small examples for triangles, quads, hex and tet elements exist in `tests/test_lean_vtk.cpp`.

### Prerequisites
* CMake >=3.10
* A C/C++ compiler with C++11 (tested with gcc 7.5.0)
### Installation

#### Copy and paste
The easiest way to use Lean-VTK is to simply copy and paste `include/lean-vtk.hpp` and `src/lean-vtk.cpp` into your project source code explicitly. If you don't mind checking for updates periodically, this is the best way to go.

#### Subrepository
 1. Add Lean-VTK as a submodule of your project and update the submodules:
 ```bash
 git submodule add https://github.com/mmorse1217/lean-vtk
 git submodule init
 ```
 2a. If you are using CMake in your project already, add the following to your root-level `CMakeLists.txt`:
 ```cmake
 add_subdirectory(lean-vtk)
 ```
 and appropriate link Lean-VTK to the project's targets:
 ```cmake
 target_link_libraries(target LeanVTK)
 ```
 2b. If your project uses make, see the compilation steps below to compile the project, then appropriately link the library.

### Compiling, installing and running
To compile the project, run the following in the project root:
```
    mkdir build
    cd build
    cmake ..
    make
```
To install the project in `/usr/local/`, run the following in the `build/` directory created above:
```
    make install
```
To run unit tests via CTest, again run the following in the `build/` directory:
```
    make test
```
or 
```
    ctest
```

### Contributing
Please fork the repo and make a pull request for any future changes.
