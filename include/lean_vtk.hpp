
#ifndef VTU_WRITER_HPP
#define VTU_WRITER_HPP

#include <string>
#include <cassert>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace leanvtk {
inline int index(int N, int i, int j) {
  assert(N > 0);
  return i * N + j;
}
template <typename T> class VTKDataNode {

public:
  VTKDataNode() {}

  VTKDataNode(const std::string &name, const std::string &numeric_type,
              const std::vector<double> &data = std::vector<double>(),
              const int n_components = 1)
      : name_(name), numeric_type_(numeric_type), data_(data),
        n_components_(n_components) {}

  inline std::vector<double> &data() { return data_; }

  void initialize(const std::string &name, const std::string &numeric_type,
                  const std::vector<double> &data, const int n_components = 1) {
    name_ = name;
    numeric_type_ = numeric_type;
    data_ = data;
    n_components_ = n_components;
  }

  void write(std::ostream &os) const {
    // NOTE this writer implicitly assumes that 2D vectors will live in 2D
    // space. To decouple this, vtk_num_components must match n_components_
    // and the conditional if(n_components_==2){...} must be corrected for the
    // context.
    int vtk_num_components = n_components_ == 2 ? n_components_ + 1 : n_components_;
    os << "<DataArray type=\"" << numeric_type_ << "\" Name=\"" << name_
       << "\" NumberOfComponents=\"" << vtk_num_components
       << "\" format=\"ascii\">\n";
    if (n_components_ != 1) {
      std::cerr << "writing matrix in vtu file (check ordering of values)"
                << std::endl;
    }

    const int num_points = data_.size() / n_components_;

    for (int d = 0; d < num_points; ++d) {
      for (int i = 0; i < n_components_; ++i) {
        int idx = index(n_components_, d, i); 
        os << data_.at(idx);
        if (i < n_components_ - 1) {
          os << " ";
        }
      }

      if(n_components_==2)
         os << " 0";

      os << "\n";
    }

    os << "</DataArray>\n";
  }

  inline bool empty() const { return data_.size() <= 0; }

private:
  std::string name_;
  
  std::string numeric_type_;
  std::vector<double> data_;
  int n_components_;
};

class VTUWriter {
public:
  /**
   * Write surface mesh to a file
   * const string& path             filename to store vtk mesh (ending with .vtu)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell 
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. Format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and 
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size}, 
   *                                  ...  
   *                                  c_{cell_size,1}, c_{cell_size,2},..., c_{cell_size, cell_size}]
   *                                  (i.e. index c*i corresponds to the ith
   *                                  vertex in the cth cell in the mesh
   */
  bool write_surface_mesh(const std::string &path,
                          const int dim,
                          const int cell_size,
                          const std::vector<double> &points,
                          const std::vector<int> &elements);

  /**
   * Write surface mesh to an output stream
   * ostream &os                    output stream where to write vtk mesh (ending with .vtu)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. Format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size},
   *                                  ...
   *                                  c_{cell_size,1}, c_{cell_size,2},..., c_{cell_size, cell_size}]
   *                                  (i.e. index c*i corresponds to the ith
   *                                  vertex in the cth cell in the mesh
   */
  bool write_surface_mesh(std::ostream &os,
                          const int dim,
                          const int cell_size,
                          const std::vector<double> &points,
                          const std::vector<int> &elements);
  /**
   * Write volume mesh to a file
   *
   * const string& path             filename to store vtk mesh (ending with .vtu)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell 
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. If there are 
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and 
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size}, 
   *                                  ...  
   *                                  c_{m,1}, c_{m,2},..., c_{m, cell_size}]
   *                                if there are m cells
   *                                (i.e. index c*i corresponds to the ith
   *                                vertex in the cth cell in the mesh
   */
  bool write_volume_mesh(const std::string &path, 
                         const int dim,
                         const int cell_size, 
                         const std::vector<double> &points,
                         const std::vector<int> &elements);

  /**
   * Write volume mesh to an output stream
   *
   * ostream &os                    output stream where to write vtk mesh (ending with .vtu)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. If there are
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size},
   *                                  ...
   *                                  c_{m,1}, c_{m,2},..., c_{m, cell_size}]
   *                                if there are m cells
   *                                (i.e. index c*i corresponds to the ith
   *                                vertex in the cth cell in the mesh
   */
  bool write_volume_mesh(std::ostream &os,
                         const int dim,
                         const int cell_size,
                         const std::vector<double> &points,
                         const std::vector<int> &elements);

  /**
   * Write point cloud to a file
   *
   * const string& path             filename to store vtk mesh (ending with .vtu)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell 
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. If there are 
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and 
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size}, 
   *                                  ...  
   *                                  c_{m,1}, c_{m,2},..., c_{m, cell_size}]
   *                                if there are m cells
   *                                (i.e. index c*i corresponds to the ith
   *                                vertex in the cth cell in the mesh
   */
  bool write_point_cloud(const std::string &path, const int dim,
                                    const std::vector<double> &points); 

  /**
   * Write point cloud to a file
   *
   * ostream &os                    output stream where to write vtk mesh (ending with .vtp)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell 
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. If there are 
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and 
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size}, 
   *                                  ...  
   *                                  c_{m,1}, c_{m,2},..., c_{m, cell_size}]
   *                                if there are m cells
   *                                (i.e. index c*i corresponds to the ith
   *                                vertex in the cth cell in the mesh
   */
  bool write_point_cloud(std::ostream &os, const int dim,
                                    const std::vector<double> &points); 

  /**
   * Add a general field to the mesh
   * const string& name             name of the field to store vtk mesh 
   * const vector<double>& data     list of field values. There must be dimension 
   *                                values for each point in the mesh to be written.
   *                                Format of the vector is 
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension}, 
   *                                  ...  
   *                                  f_{n,1}, f_{n,2},..., f_{n, dimension}]
   *                                if there are n points in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  void add_field(const std::string &name, 
                 const std::vector<double> &data,
                 const int &dimension);

  /**
   * Add a general cell/element field to the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be dimension
   *                                values for each cell in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension},
   *                                  ...
   *                                  f_{m,1}, f_{m,2},..., f_{m, dimension}]
   *                                if there are m cells in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  void add_cell_field(const std::string &name,
                 const std::vector<double> &data,
                 const int &dimension);

  /**
   * Add a scalar field to the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be one
   *                                value for each point in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_1, f_2,..., f_n]
   *                                if there are n points in the mesh
   */
  void add_scalar_field(const std::string &name,
                        const std::vector<double> &data);

  /**
   * Add a scalar field to cells/elements of the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be one
   *                                value for each cell in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_1, f_2,..., f_m]
   *                                if there are m cells in the mesh
   */
  void add_cell_scalar_field(const std::string &name,
                        const std::vector<double> &data);

  /**
   * Add a vector field to the mesh
   * const string& name             name of the field to store vtk mesh 
   * const vector<double>& data     list of field values. There must be dimension 
   *                                values for each point in the mesh to be written.
   *                                Format of the vector is 
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension}, 
   *                                  ...  
   *                                  f_{n,1}, f_{n,2},..., f_{n, dimension}]
   *                                if there are n points in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  void add_vector_field(const std::string &name,
                        const std::vector<double> &data, 
                        const int &dimension);

  /**
   * Add a vector field to cells/elements of the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be dimension
   *                                values for each cell in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension},
   *                                  ...
   *                                  f_{m,1}, f_{m,2},..., f_{m, dimension}]
   *                                if there are m cells in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  void add_cell_vector_field(const std::string &name,
                        const std::vector<double> &data,
                        const int &dimension);

  // Remove all fields and initialized data from the writer.
  void clear();

private:
  std::vector<VTKDataNode<double>> point_data_;
  std::vector<VTKDataNode<double>> cell_data_;
  std::string current_scalar_point_data_;
  std::string current_vector_point_data_;
  std::string current_scalar_cell_data_;
  std::string current_vector_cell_data_;

  
  void write_point_data(std::ostream &os);

  void write_cell_data(std::ostream &os);

  void write_header(const int n_vertices, const int n_elements,
                    std::ostream &os);

  void write_footer(std::ostream &os);
  
  bool write_mesh(std::ostream &os, const int dim, const int cell_size,
                  const std::vector<double> &points, const std::vector<int> &tets, 
                  bool is_volume_mesh=true);

  bool write_mesh(const std::string &path, const int dim, const int cell_size,
                  const std::vector<double> &points, const std::vector<int> &tets, 
                  bool is_volume_mesh=true);

  void write_points(const int num_points, const std::vector<double> &points,
                    std::ostream &os, bool is_volume_mesh = true);

  void write_cells(const int n_vertices, const std::vector<int> &tets,
                   std::ostream &os, bool is_volume_mesh = true);
  };
} // namespace leanvtk

#endif // VTU_WRITER_HPP

