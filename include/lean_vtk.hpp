
#ifndef VTU_WRITER_HPP
#define VTU_WRITER_HPP

#include <string>
#include <cassert>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>

namespace leanvtk {
inline int index(int N, int i, int j) {
  assert(N > 0);
  return i * N + j;
}

static const char base64_encoding_table[] = {
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
  'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
  'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
  'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
  'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
  'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
  'w', 'x', 'y', 'z', '0', '1', '2', '3',
  '4', '5', '6', '7', '8', '9', '+', '/'};

static const int base64_mod_table[] = {0, 2, 1};

inline std::string base64_encode(const unsigned char *data,
                                 size_t input_length) {
    size_t output_length = 4 * ((input_length + 2) / 3);
    std::string encoded_data;
    encoded_data.resize(output_length);

    for (int i = 0, j = 0; i < input_length;) {
        uint32_t octet_a = i < input_length ? (unsigned char)data[i++] : 0;
        uint32_t octet_b = i < input_length ? (unsigned char)data[i++] : 0;
        uint32_t octet_c = i < input_length ? (unsigned char)data[i++] : 0;

        uint32_t triple = (octet_a << 0x10) + (octet_b << 0x08) + octet_c;

        encoded_data[j++] = base64_encoding_table[(triple >> 3 * 6) & 0x3F];
        encoded_data[j++] = base64_encoding_table[(triple >> 2 * 6) & 0x3F];
        encoded_data[j++] = base64_encoding_table[(triple >> 1 * 6) & 0x3F];
        encoded_data[j++] = base64_encoding_table[(triple >> 0 * 6) & 0x3F];
    }

    for (int i = 0; i < base64_mod_table[input_length % 3]; i++)
        encoded_data[output_length - 1 - i] = '=';

    return encoded_data;
}


class VTKDataNodeBase {
public:
  VTKDataNodeBase(const std::string &name="") 
      : name_(name)
      , binary_(false)
  {}

  virtual void write(std::ostream &os) const{};

  /// Set the format to binary
  inline void set_binary() { binary_ = true; }

  /// Set the format to ASCII
  inline void set_ascii() { binary_ = false; }

  /// Set whether binary format is considered or not
  inline void set_binary(bool enable) { binary_ = enable; }

  /// Get if binary format is enabled
  inline bool is_binary() { return binary_; }
protected:
  std::string name_;
  bool binary_;
};

template <typename T>
class VTKDataNode : public VTKDataNodeBase {

public:
  VTKDataNode()
      : VTKDataNodeBase()
  {}

  VTKDataNode(const std::string &name, const std::string &numeric_type,
              const std::vector<T> &data = std::vector<double>(),
              const int n_components = 1)
      : VTKDataNodeBase(name)
      , numeric_type_(numeric_type)
      , data_(data)
      , n_components_(n_components)
  {}

  inline std::vector<T> &data() { return data_; }

  void initialize(const std::string &name, const std::string &numeric_type,
                  const std::vector<T> &data, const int n_components = 1) {
    name_ = name;
    numeric_type_ = numeric_type;
    data_ = data;
    n_components_ = n_components;
  }

  void write(std::ostream &os) const {
    os << "<DataArray type=\"" << numeric_type_ << "\" Name=\"" << name_
       << "\" NumberOfComponents=\"" << n_components_
       << "\" format=\"" << (binary_ ? "binary" : "ascii") << "\">\n";
    if (n_components_ != 1) {
      std::cerr << "writing matrix in vtu file (check ordering of values)"
                << std::endl;
    }


    if (binary_) {
      os << base64_encode((unsigned char*)data_.data(),
                          sizeof(T) * data_.size());
    } else {
      const int num_points = data_.size() / n_components_;
      for (int d = 0; d < num_points; ++d) {
        for (int i = 0; i < n_components_; ++i) {
          int idx = index(n_components_, d, i); 
          os << data_.at(idx);
          if (i < n_components_ - 1) {
            os << " ";
          }
        }
        os << "\n";
      }
    }

    os << "</DataArray>\n";
  }

  inline bool empty() const { return data_.size() <= 0; }

private:
  std::string numeric_type_;
  std::vector<T> data_;
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
   * const vector<double>& points   list of point locations. If there are 
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and 
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
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
  template <typename T>
  inline void add_field(const std::string &name, 
                        const std::vector<T> &data,
                        const int &dimension)
  {
    if (dimension == 1)
      add_scalar_field<T>(name, data);
    else
      add_vector_field<T>(name, data, dimension);
  }

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
  template <typename T>
  void add_cell_field(const std::string &name,
                 const std::vector<T> &data,
                 const int &dimension)
  {
    if (dimension == 1)
      add_cell_scalar_field<T>(name, data);
    else
      add_cell_vector_field<T>(name, data, dimension);
  }

  /**
   * Add a scalar field to the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be one
   *                                value for each point in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_1, f_2,..., f_n]
   *                                if there are n points in the mesh
   */
  template <typename T>
  void add_scalar_field(const std::string &name,
                        const std::vector<T> &data);

  /**
   * Add a scalar field to cells/elements of the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be one
   *                                value for each cell in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_1, f_2,..., f_m]
   *                                if there are m cells in the mesh
   */
  template <typename T>
  void add_cell_scalar_field(const std::string &name,
                        const std::vector<T> &data);

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
  template <typename T>
  void add_vector_field(const std::string &name,
                        const std::vector<T> &data, 
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
   *                                if there are m bool binary = falsecells in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  template <typename T>
  void add_cell_vector_field(const std::string &name,
                        const std::vector<T> &data,
                        const int &dimension);

  // Remove all fields and initialized data from the writer.
  void clear();

  /// Set the format to binary
  inline void set_binary() { binary_ = true; }

  /// Set the format to ASCII
  inline void set_ascii() { binary_ = false; }

  /// Set whether binary format is considered or not
  inline void set_binary(bool enable) { binary_ = enable; }

  /// Get if binary format is enabled
  inline bool is_binary() { return binary_; }
private:
  std::vector<VTKDataNodeBase> point_data_;
  std::vector<VTKDataNodeBase> cell_data_;
  std::string current_scalar_point_data_;
  std::string current_vector_point_data_;
  std::string current_scalar_cell_data_;
  std::string current_vector_cell_data_;
  bool binary_;

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

  template <typename T>
  inline VTKDataNode<T> make_data_node(const std::string &name,
                                       const std::vector<T> &data,
                                       std::string num_type="Float",
                                       const int dimension=1)
  {
    VTKDataNode<T> node;
    node.initialize(name,
                    num_type + std::to_string(8 * sizeof(T)),
                    data,
                    dimension);
    return node;
  }
};

} // namespace leanvtk

#endif // VTU_WRITER_HPP

