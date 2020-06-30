#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <lean_vtk.hpp>
#include <string>
#include <vector>
#include <complex>
using std::cout;
using std::endl;
using std::vector;

TEST_CASE("Test single flat quad", "[single-quad]"){
    vector<double> points = {
         1.,  1., 0.,
        -1.,  1., 0.,
        -1., -1., 0.,
         1., -1., 0.
    };
    vector<int> quads = {
        0, 1, 2, 3
    };
    polyfem::VTUWriter writer;

    const int dim = 3;
    const int cell_size = 4;
    writer.write_tet_mesh("test.vtu", dim , cell_size, points, quads);

}
