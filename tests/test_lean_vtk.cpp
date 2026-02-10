#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <leanvtk.hpp>
#include <string>
#include <vector>
#include <complex>
using std::cout;
using std::endl;
using std::vector;
using  leanvtk::VTUWriter;
typedef double real_t;
TEST_CASE("Single elements", "[single]"){
    vector<real_t> points;
    vector<size_t> elements;
    vector<real_t> scalar_field;
    vector<real_t> vector_field;
    vector<real_t> cell_scalar_field;
    vector<real_t> cell_vector_field;
    size_t dim;
    size_t cell_size;
    std::string filename;
    bool is_volume_mesh = false;
    SECTION("Triangle 2D"){
        points = {
             1.,  1.,
             1., -1.,
            -1., -1.
        };
        elements = {
            0, 1, 2
        };
        VTUWriter writer;

        dim = 2;
        cell_size = 3;
        scalar_field = {
            0., 1., 2.
        };
        vector_field = points;
        cell_scalar_field = {1.};
        cell_vector_field = {1., 1.};

        filename = "single_tri_2D.vtu";
    }
    SECTION("Quad 2D"){
        points = {
             1.,  1., 
             1., -1., 
            -1., -1., 
            -1.,  1. 
        };
        elements = {
            0, 1, 2, 3
        };

        scalar_field = {
            0., 1., 2., 3.
        };
        vector_field = points;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1.};

        dim = 2;
        cell_size = 4;
        filename = "single_quad_2D.vtu";
    }
    SECTION("Triangle 3D"){
        points = {
             1.,  1., -1.,
             1., -1., 1.,
            -1., -1., 0.
        };
        elements = {
            0, 1, 2
        };
        VTUWriter writer;

        dim = 3;
        cell_size = 3;
        scalar_field = {
            0., 1., 2.
        };
        vector_field = points;;
        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        filename = "single_tri.vtu";
    }
    SECTION("Quad 3D"){
        points = {
             1.,  1., 1.,
             1., -1., 0.,
            -1., -1., 0.,
            -1.,  1., -1.
        };
        elements = {
            0,  1,2,  3
        };

        scalar_field = {
            0., 1., 2., 3.
        };
        vector_field = points;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        dim = 3;
        cell_size = 4;
        filename = "single_quad.vtu";
    }
    SECTION("Hex element"){
        points = {
             1.,  1., 1.,
             1., -1., 1.,
            -1., -1., 1.,
            -1.,  1., 1.,
             1.,  1.,-1.,
             1., -1.,-1.,
            -1., -1.,-1.,
            -1.,  1.,-1.
        };
        elements= {0, 1, 2, 3, 4, 5, 6, 7};
        
        scalar_field = { 0., 1., 2., 3., 4., 5., 6., 7.  };
        vector_field = points;;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        dim = 3;
        cell_size = 8;
        filename = "single_hex.vtu";
        is_volume_mesh = true;
    }

    SECTION("Tet element"){
        points = {
             0.,  1., 0.,
             1.,  0., 0.,
             0.,  0., 0.,
             0.,  0., 1.
        };
        elements = { 0, 1, 2, 3, };
        
        scalar_field = { 0., 1., 2., 3.  };
        vector_field = points;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        dim = 3;
        cell_size = 4;
        filename = "single_tet.vtu";
        is_volume_mesh = true;
    }
    VTUWriter writer;

    writer.add_scalar_field("scalar_field", scalar_field);
    writer.add_vector_field("vector_field", vector_field, dim);

    writer.add_cell_scalar_field("cell_scalar_field", cell_scalar_field);
    writer.add_cell_vector_field("cell_vector_field", cell_vector_field, dim);

    if (is_volume_mesh) {
        writer.write_volume_mesh(filename, dim, cell_size, points, elements);
    } else {
        writer.write_surface_mesh(filename, dim, cell_size, points, elements);
    }
    
}

TEST_CASE("Single binary elements", "[single_bin]"){
    vector<real_t> points;
    vector<size_t> elements;
    vector<real_t> scalar_field;
    vector<real_t> vector_field;
    vector<real_t> cell_scalar_field;
    vector<real_t> cell_vector_field;
    size_t dim;
    size_t cell_size;
    std::string filename;
    bool is_volume_mesh = false;
    SECTION("Triangle 2D"){
        points = {
             1.,  1.,
             1., -1.,
            -1., -1.
        };
        elements = {
            0, 1, 2
        };
        VTUWriter writer;

        dim = 2;
        cell_size = 3;
        scalar_field = {
            0., 1., 2.
        };
        vector_field = points;
        cell_scalar_field = {1.};
        cell_vector_field = {1., 1.};

        filename = "single_bin_tri_2D.vtu";
    }
    SECTION("Quad 2D"){
        points = {
             1.,  1., 
             1., -1., 
            -1., -1., 
            -1.,  1. 
        };
        elements = {
            0, 1, 2, 3
        };

        scalar_field = {
            0., 1., 2., 3.
        };
        vector_field = points;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1.};

        dim = 2;
        cell_size = 4;
        filename = "single_bin_quad_2D.vtu";
    }
    SECTION("Triangle 3D"){
        points = {
             1.,  1., -1.,
             1., -1., 1.,
            -1., -1., 0.
        };
        elements = {
            0, 1, 2
        };
        VTUWriter writer;

        dim = 3;
        cell_size = 3;
        scalar_field = {
            0., 1., 2.
        };
        vector_field = points;;
        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        filename = "single_bin_tri.vtu";
    }
    SECTION("Quad 3D"){
        points = {
             1.,  1., 1.,
             1., -1., 0.,
            -1., -1., 0.,
            -1.,  1., -1.
        };
        elements = {
            0,  1,2,  3
        };

        scalar_field = {
            0., 1., 2., 3.
        };
        vector_field = points;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        dim = 3;
        cell_size = 4;
        filename = "single_bin_quad.vtu";
    }
    SECTION("Hex element"){
        points = {
             1.,  1., 1.,
             1., -1., 1.,
            -1., -1., 1.,
            -1.,  1., 1.,
             1.,  1.,-1.,
             1., -1.,-1.,
            -1., -1.,-1.,
            -1.,  1.,-1.
        };
        elements= {0, 1, 2, 3, 4, 5, 6, 7};
        
        scalar_field = { 0., 1., 2., 3., 4., 5., 6., 7.  };
        vector_field = points;;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        dim = 3;
        cell_size = 8;
        filename = "single_bin_hex.vtu";
        is_volume_mesh = true;
    }

    SECTION("Tet element"){
        points = {
             0.,  1., 0.,
             1.,  0., 0.,
             0.,  0., 0.,
             0.,  0., 1.
        };
        elements = { 0, 1, 2, 3, };
        
        scalar_field = { 0., 1., 2., 3.  };
        vector_field = points;

        cell_scalar_field = {1.};
        cell_vector_field = {1., 1., 1.};

        dim = 3;
        cell_size = 4;
        filename = "single_bin_tet.vtu";
        is_volume_mesh = true;
    }
    VTUWriter writer;

    writer.add_scalar_field("scalar_field", scalar_field);
    writer.add_vector_field("vector_field", vector_field, dim);

    writer.add_cell_scalar_field("cell_scalar_field", cell_scalar_field);
    writer.add_cell_vector_field("cell_vector_field", cell_vector_field, dim);

    writer.set_binary();
    if (is_volume_mesh) {
        writer.write_volume_mesh(filename, dim, cell_size, points, elements);
    } else {
        writer.write_surface_mesh(filename, dim, cell_size, points, elements);
    }
    
}

TEST_CASE("Multiple elements", "[mesh]"){
    vector<real_t> points;
    vector<size_t> elements;
    vector<real_t> scalar_field;
    vector<real_t> vector_field;
    vector<real_t> cell_scalar_field;
    vector<real_t> cell_vector_field;
    size_t dim;
    size_t cell_size;
    std::string filename;
    bool is_volume_mesh = false;
    SECTION("Two triangle 2D"){
        points = {
            0., 0.,
            0., 1., 
            1., 0.,
            1., 1.
        };
        elements = {
            0, 1, 2,
            1, 3, 2 
        };
        VTUWriter writer;

        dim = 2;
        cell_size = 3;
        scalar_field = {
            0., 1., 2., 3.
        };
        vector_field = points;
        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0.,
                             0., 1.};
        filename = "mesh_tri_2D.vtu";
    }
    SECTION("Two quad 2D"){
        points = {
             1.,  1., 
             1., -1., 
            -1., -1., 
            -1.,  1.,
            3., 1.,
            3., -1.
        };
        elements = {
            0, 1, 2, 3,
            0, 4, 5, 1
        };

        scalar_field = {
            0., 1., 2., 3., 4., 5.
        };
        vector_field = points;

        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0.,
                             0., 1.};

        dim = 2;
        cell_size = 4;
        filename = "mesh_quad_2D.vtu";
    }
    SECTION("Two triangles 3D"){
        points = {
             1.,  1., -1.,
             1., -1., 1.,
            -1., -1., 0.,
            2., 1., 1.
        };
        elements = {
            0, 1, 2,
            1, 3, 2 
        };
        VTUWriter writer;

        dim = 3;
        cell_size = 3;
        scalar_field = {
            0., 1., 2., 3. 
        };
        vector_field = points;
        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};
        filename = "mesh_tri.vtu";
    }
    SECTION("Two quads 3D"){
        points = {
             1.,  1., 1.,
             1., -1., 0.,
            -1., -1., 0.,
            -1.,  1., -1.,
             3., 1., 2.,
             3., -1., 2.
        };
        elements = {
            0, 1, 2, 3,
            0, 4, 5, 1
        };

        scalar_field = {
            0., 1., 2., 3., 4., 5.
        };
        vector_field = points;

        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};

        dim = 3;
        cell_size = 4;
        filename = "mesh_quad.vtu";
    }
    SECTION("Two hex elements"){
        points = {
             1.,  1., 1.,
             1., -1., 1.,
            -1., -1., 1.,
            -1.,  1., 1.,
             1.,  1.,-1.,
             1., -1.,-1.,
            -1., -1.,-1.,
            -1.,  1.,-1.,
             1.,  1., 3.,
             1., -1., 3.,
            -1., -1., 3.,
            -1.,  1., 3.,

        };
        elements= {
            0, 1, 2, 3, 4, 5, 6, 7,
            8, 9, 10, 11,0, 1, 2, 3,
        };
        
        scalar_field = { 0., 1., 2., 3., 4., 5., 6., 7., 8.,9., 10., 11. };
        vector_field = points;

        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};

        dim = 3;
        cell_size = 8;
        filename = "mesh_hex.vtu";
        is_volume_mesh = true;
    }

    SECTION(" Two tet element"){
        points = {
             0.,  1., 0.,
             1.,  0., 0.,
             0.,  0., 0.,
             0.,  0., 1.,
             1., 1., 1.

        };
        elements = { 
            0, 1, 2, 3, 
            0, 1, 3, 4
        };
        
        scalar_field = { 0., 1., 2., 3.,4  };
        vector_field = points;
        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};

        dim = 3;
        cell_size = 4;
        filename = "mesh_tet.vtu";
        is_volume_mesh = true;
    }
    VTUWriter writer;

    writer.add_scalar_field("scalar_field", scalar_field);
    writer.add_vector_field("vector_field", vector_field, dim);
    writer.add_cell_scalar_field("cell_scalar_field", cell_scalar_field);
    writer.add_cell_vector_field("cell_vector_field", cell_vector_field, dim);
    if (is_volume_mesh) {
        writer.write_volume_mesh(filename, dim, cell_size, points, elements);
    } else {
        writer.write_surface_mesh(filename, dim, cell_size, points, elements);
    }
    
}

TEST_CASE("Multiple binary elements", "[mesh_bin]"){
    vector<real_t> points;
    vector<size_t> elements;
    vector<real_t> scalar_field;
    vector<real_t> vector_field;
    vector<real_t> cell_scalar_field;
    vector<real_t> cell_vector_field;
    size_t dim;
    size_t cell_size;
    std::string filename;
    bool is_volume_mesh = false;
    SECTION("Two triangle 2D"){
        points = {
            0., 0.,
            0., 1., 
            1., 0.,
            1., 1.
        };
        elements = {
            0, 1, 2,
            1, 3, 2 
        };
        VTUWriter writer;

        dim = 2;
        cell_size = 3;
        scalar_field = {
            0., 1., 2., 3.
        };
        vector_field = points;
        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0.,
                             0., 1.};
        filename = "mesh_bin_tri_2D.vtu";
    }
    SECTION("Two quad 2D"){
        points = {
             1.,  1., 
             1., -1., 
            -1., -1., 
            -1.,  1.,
            3., 1.,
            3., -1.
        };
        elements = {
            0, 1, 2, 3,
            0, 4, 5, 1
        };

        scalar_field = {
            0., 1., 2., 3., 4., 5.
        };
        vector_field = points;

        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0.,
                             0., 1.};

        dim = 2;
        cell_size = 4;
        filename = "mesh_bin_quad_2D.vtu";
    }
    SECTION("Two triangles 3D"){
        points = {
             1.,  1., -1.,
             1., -1., 1.,
            -1., -1., 0.,
            2., 1., 1.
        };
        elements = {
            0, 1, 2,
            1, 3, 2 
        };
        VTUWriter writer;

        dim = 3;
        cell_size = 3;
        scalar_field = {
            0., 1., 2., 3. 
        };
        vector_field = points;
        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};
        filename = "mesh_bin_tri.vtu";
    }
    SECTION("Two quads 3D"){
        points = {
             1.,  1., 1.,
             1., -1., 0.,
            -1., -1., 0.,
            -1.,  1., -1.,
             3., 1., 2.,
             3., -1., 2.
        };
        elements = {
            0, 1, 2, 3,
            0, 4, 5, 1
        };

        scalar_field = {
            0., 1., 2., 3., 4., 5.
        };
        vector_field = points;

        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};

        dim = 3;
        cell_size = 4;
        filename = "mesh_bin_quad.vtu";
    }
    SECTION("Two hex elements"){
        points = {
             1.,  1., 1.,
             1., -1., 1.,
            -1., -1., 1.,
            -1.,  1., 1.,
             1.,  1.,-1.,
             1., -1.,-1.,
            -1., -1.,-1.,
            -1.,  1.,-1.,
             1.,  1., 3.,
             1., -1., 3.,
            -1., -1., 3.,
            -1.,  1., 3.,

        };
        elements= {
            0, 1, 2, 3, 4, 5, 6, 7,
            8, 9, 10, 11,0, 1, 2, 3,
        };
        
        scalar_field = { 0., 1., 2., 3., 4., 5., 6., 7., 8.,9., 10., 11. };
        vector_field = points;

        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};

        dim = 3;
        cell_size = 8;
        filename = "mesh_bin_hex.vtu";
        is_volume_mesh = true;
    }

    SECTION(" Two tet element"){
        points = {
             0.,  1., 0.,
             1.,  0., 0.,
             0.,  0., 0.,
             0.,  0., 1.,
             1., 1., 1.

        };
        elements = { 
            0, 1, 2, 3, 
            0, 1, 3, 4
        };
        
        scalar_field = { 0., 1., 2., 3.,4  };
        vector_field = points;
        cell_scalar_field = {1., 2.};
        cell_vector_field = {1., 0., 0.,
                             0., 1., 0.};

        dim = 3;
        cell_size = 4;
        filename = "mesh_bin_tet.vtu";
        is_volume_mesh = true;
    }
    VTUWriter writer;

    writer.add_scalar_field("scalar_field", scalar_field);
    writer.add_vector_field("vector_field", vector_field, dim);
    writer.add_cell_scalar_field("cell_scalar_field", cell_scalar_field);
    writer.add_cell_vector_field("cell_vector_field", cell_vector_field, dim);
    writer.set_binary();
    if (is_volume_mesh) {
        writer.write_volume_mesh(filename, dim, cell_size, points, elements);
    } else {
        writer.write_surface_mesh(filename, dim, cell_size, points, elements);
    }
    
}

TEST_CASE("Test point clouds", "[point]"){
    vector<real_t> points;
    vector<real_t> scalar_field;
    vector<real_t> vector_field;
    size_t dim;
    size_t cell_size;
    std::string filename;
    SECTION("2D single point"){
        dim = 2;
        points = { 0., 1.};
        scalar_field = { 0.};
        vector_field = points;
        filename = "single_point_2d.vtu";

    }
    SECTION("2D multiple points"){
        points = {
            0., 0.,
            0., 1., 
            1., 0.,
            1., 1.
        };
        dim = 2;
        scalar_field = { 0., 1.,2.,3.};
        vector_field = points;
        filename = "multi_point_2d.vtu";

    }
   SECTION("3D single point"){
        points = {
            0., 0., 0.
        };
        dim = 3;
        scalar_field = { 0.};
        vector_field = points;
        filename = "single_point_3d.vtu";


    }
    SECTION("3D multiple points"){
        points = {
            0., 0., 1.,
            0., 1., 0.,
            1., 0., 1.,
            1., 1., 1.
        };
        dim = 3;
        scalar_field = { 0., 1.,2.,3.};
        vector_field = points;
        filename = "multi_point_3d.vtu";


    }
    VTUWriter writer;

    writer.add_scalar_field("scalar_field", scalar_field);
    writer.add_vector_field("vector_field", vector_field, dim);
    writer.write_point_cloud(filename, dim, points);
}
