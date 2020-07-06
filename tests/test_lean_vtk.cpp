#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <lean_vtk.hpp>
#include <string>
#include <vector>
#include <complex>
using std::cout;
using std::endl;
using std::vector;

TEST_CASE("Test single flat tri", "[single][tri]"){
    vector<double> points = {
         1.,  1., -1.,
         1., -1., 1.,
        -1., -1., 0.
    };
    vector<int> tris = {
        0,
        1,
        2
    };
    polyfem::VTUWriter writer;

    const int dim = 3;
    const int cell_size = 3;
    writer.write_surface_mesh("test_tri.vtu", dim, cell_size, points, tris);
}
TEST_CASE("Test single flat quad", "[single][quad]"){
    vector<double> points = {
         1.,  1., 1.,
         1., -1., 0.,
        -1., -1., 0.,
        -1.,  1., -1.
    };
    vector<int> quads = {
        0,  1,2,  3
    };
    polyfem::VTUWriter writer;

    const int dim = 3;
    const int cell_size = 4;
    writer.write_surface_mesh("test_quad.vtu", dim , cell_size, points, quads);

}


TEST_CASE("Test single flat hex", "[single][hex]"){
    vector<double> points = {
         1.,  1., 1.,
         1., -1., 1.,
        -1., -1., 1.,
        -1.,  1., 1.,
         1.,  1.,-1.,
         1., -1.,-1.,
        -1., -1.,-1.,
        -1.,  1.,-1.
    };
    vector<int> hexes = {0, 1, 2, 3, 4, 5, 6, 7};
    polyfem::VTUWriter writer;

    const int dim = 3;
    const int cell_size = 8;
    writer.write_volume_mesh("test_hex.vtu", dim, cell_size, points, hexes);
}

TEST_CASE("Test single flat tetrahedron", "[single][tet]"){
    vector<double> points = {
         0.,  1., 0.,
         1.,  0., 0.,
         0.,  0., 0.,
         0.,  0., 1.
    };
    vector<int> tets = {
        0,
        1,
        2,
        3,
    };
    polyfem::VTUWriter writer;

    const int dim = 3;
    const int cell_size = 4;
    writer.write_volume_mesh("test_tet.vtu", dim, cell_size, points, tets);
}
