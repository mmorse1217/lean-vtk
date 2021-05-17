#include <lean_vtk.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
std::string a_library_function(){
    return std::string("a function specified in the source code");
}

namespace leanvtk {

static const int VTK_TETRA = 10;
static const int VTK_TRIANGLE = 5;
static const int VTK_QUAD = 9;
static const int VTK_HEXAHEDRON = 12;
static const int VTK_POLYGON = 7;

inline static int VTKTagVolume(const int n_vertices) {
  switch (n_vertices) {
  case 4:
    return VTK_TETRA;
  case 8:
    return VTK_HEXAHEDRON;
  default:
    // element type not supported. To add it
    // (http://www.vtk.org/VTK/img/file-formats.pdf)
    cerr << n_vertices << " not supported, " << endl;
    assert(false);
    return -1;
  }
}

inline static int VTKTagPlanar(const int n_vertices) {
  switch (n_vertices) {
  case 3:
    return VTK_TRIANGLE;
  case 4:
    return VTK_QUAD;
  default:
    // element type not supported. To add it
    // (http://www.vtk.org/VTK/img/file-formats.pdf)
    cerr << "{} not supported, " << n_vertices << endl;
    assert(false);
    return -1;
  }
}

void VTUWriter::write_point_data(std::ostream &os) {
  if (current_scalar_point_data_.empty() && current_vector_point_data_.empty())
    return;

  os << "<PointData ";
  if (!current_scalar_point_data_.empty())
    os << "Scalars=\"" << current_scalar_point_data_ << "\" ";
  if (!current_vector_point_data_.empty())
    os << "Vectors=\"" << current_vector_point_data_ << "\" ";
  os << ">\n";

  for (auto it = point_data_.begin(); it != point_data_.end(); ++it) {
    it->write(os);
  }

  os << "</PointData>\n";
}

void VTUWriter::write_header(const int n_vertices, const int n_elements,
                             std::ostream &os) {
  os << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
  os << "<UnstructuredGrid>\n";
  os << "<Piece NumberOfPoints=\"" << n_vertices << "\" NumberOfCells=\""
     << n_elements << "\">\n";
}

void VTUWriter::write_footer(std::ostream &os) {
  os << "</Piece>\n";
  os << "</UnstructuredGrid>\n";
  os << "</VTKFile>\n";
}

void VTUWriter::write_points(const int num_points, const vector<double> &points,
                             std::ostream &os, bool is_volume_mesh) {
  os << "<Points>\n";
  os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        "format=\"ascii\">\n";
  const int dim = points.size() / num_points;
  assert(double(dim) == double(points.size()) / double(num_points));

  for (int d = 0; d < num_points; ++d) {
    for (int i = 0; i < dim; ++i) {
      int idx = index(dim, d, i); 
      os << points.at(idx);
      if (i < dim - 1) {
        os << " ";
      }
    }

    if (!is_volume_mesh && dim == 2)
      os << " 0";

    os << "\n";
  }

  os << "</DataArray>\n";
  os << "</Points>\n";
}

void VTUWriter::write_cells(const int n_vertices, const vector<int> &tets,
                            std::ostream &os, bool is_volume_mesh) {
  const int n_cells = tets.size() / n_vertices;
  os << "<Cells>\n";
  /////////////////////////////////////////////////////////////////////////////
  // List vertex id's i=0, ..., n_vertices associated with each cell c
  os << "<DataArray type=\"Int64\" Name=\"connectivity\" "
        "format=\"ascii\">\n";
  for (int c = 0; c < n_cells; ++c) {
    for (int i = 0; i < n_vertices; ++i) {
      int idx = index(n_vertices, c, i);
      const int v_index = tets.at(idx);
      os << v_index;
      if (i < n_vertices - 1) {
        os << " ";
      }
    }
    os << "\n";
  }

  os << "</DataArray>\n";
  /////////////////////////////////////////////////////////////////////////////
  // List the VTK cell type for each mesh element.
  // This assumes a uniform cell type the entire mesh; to generalize, pass
  // or compute the number of vertices per cell and recompute the cell type
  int cell_type =
      is_volume_mesh ? VTKTagVolume(n_vertices) : VTKTagPlanar(n_vertices);

  os << "<DataArray type=\"Int8\" Name=\"types\" format=\"ascii\" "
        "RangeMin=\""
     << cell_type << "\" RangeMax=\"" << cell_type << "\">\n";
  for (int i = 0; i < n_cells; ++i) {
    os << cell_type << "\n";
  }
  os << "</DataArray>\n";

  /////////////////////////////////////////////////////////////////////////////
  // List offsets to access the vertex indices of the ith cell. Non-trivial
  // if the mesh is a general polyognal mesh.
  os << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" "
        "RangeMin=\""
     << n_vertices << "\" RangeMax=\"" << n_cells * n_vertices << "\">\n";

  int acc = n_vertices;
  for (int i = 0; i < n_cells; ++i) {
    os << acc << "\n";
    acc += n_vertices;
  }

  os << "</DataArray>\n";
  /////////////////////////////////////////////////////////////////////////////
  os << "</Cells>\n";
}

void VTUWriter::write_cell_data(std::ostream &os) {
  if (current_scalar_cell_data_.empty() && current_vector_cell_data_.empty())
    return;

  os << "<CellData ";
  if (!current_scalar_cell_data_.empty())
    os << "Scalars=\"" << current_scalar_cell_data_ << "\" ";
  if (!current_vector_cell_data_.empty())
    os << "Vectors=\"" << current_vector_cell_data_ << "\" ";
  os << ">\n";

  for (auto it = cell_data_.begin(); it != cell_data_.end(); ++it) {
    it->write(os);
  }

  os << "</CellData>\n";
}
void VTUWriter::clear() {
  point_data_.clear();
  cell_data_.clear();
}

void VTUWriter::add_field(const std::string &name, const vector<double> &data,
                          const int &dimension) {
  using std::abs;

  vector<double> tmp(data.size(), 0.);

  for (long i = 0; i < data.size(); ++i)
    tmp[i] = std::abs(data[i]) < 1e-16 ? 0 : data[i];

  if (dimension == 1)
    add_scalar_field(name, tmp);
  else
    add_vector_field(name, tmp, dimension);
}

void VTUWriter::add_scalar_field(const std::string &name,
                                 const vector<double> &data) {
  point_data_.push_back(VTKDataNode<double>());
  point_data_.back().initialize(name, "Float64", data);
  current_scalar_point_data_ = name;
}

void VTUWriter::add_vector_field(const std::string &name,
                                 const vector<double> &data,
                                 const int &dimension) {
  point_data_.push_back(VTKDataNode<double>());

  point_data_.back().initialize(name, "Float64", data, dimension);
  current_vector_point_data_ = name;
}

void VTUWriter::add_cell_scalar_field(const std::string &name,
                                 const vector<double> &data) {
  cell_data_.push_back(VTKDataNode<double>());
  cell_data_.back().initialize(name, "Float64", data);
  current_scalar_cell_data_ = name;
}

void VTUWriter::add_cell_vector_field(const std::string &name,
                                 const vector<double> &data,
                                 const int &dimension) {
  cell_data_.push_back(VTKDataNode<double>());

  cell_data_.back().initialize(name, "Float64", data, dimension);
  current_vector_cell_data_ = name;
}

void VTUWriter::add_cell_field(const std::string &name, const vector<double> &data,
                          const int &dimension) {
  using std::abs;

  vector<double> tmp(data.size(), 0.);

  for (long i = 0; i < data.size(); ++i)
    tmp[i] = std::abs(data[i]) < 1e-16 ? 0 : data[i];

  if (dimension == 1)
    add_cell_scalar_field(name, tmp);
  else
    add_cell_vector_field(name, tmp, dimension);
}

bool VTUWriter::write_mesh(std::ostream &os, const int dim, const int cell_size,
              const vector<double> &points, const vector<int> &tets, bool is_volume_mesh){
  assert(dim > 1);
  assert(cell_size > 2);

  int num_points = points.size() / dim;
  int num_cells = tets.size() / cell_size;

  write_header(num_points, num_cells, os);
  write_points(num_points, points, os, is_volume_mesh);
  write_point_data(os);
  write_cells(cell_size, tets, os, is_volume_mesh);
  write_cell_data(os);
  write_footer(os);
  clear();
  return true;
}

bool VTUWriter::write_mesh(const std::string &path, const int dim, const int cell_size,
              const vector<double> &points, const vector<int> &tets, bool is_volume_mesh){


  std::ofstream os;
  os.open(path.c_str());
  if (!os.good()) {
    os.close();
    return false;
  }

  write_mesh(os, dim, cell_size, points, tets, is_volume_mesh);

  os.close();
  return true;
}

bool VTUWriter::write_surface_mesh(const std::string &path, const int dim,
                                   const int cell_size,
                                   const vector<double> &points,
                                   const vector<int> &tets) {
  
    return write_mesh(path, dim, cell_size, points, tets, false);
}

bool VTUWriter::write_volume_mesh(const std::string &path, const int dim,
                                  const int cell_size,
                                  const vector<double> &points,
                                  const vector<int> &tets) {
  
    return write_mesh(path, dim, cell_size, points, tets, true);
}

bool VTUWriter::write_surface_mesh(std::ostream &os, const int dim,
                                   const int cell_size,
                                   const vector<double> &points,
                                   const vector<int> &tets) {
  
    return write_mesh(os, dim, cell_size, points, tets, false);
}

bool VTUWriter::write_volume_mesh(std::ostream &os, const int dim,
                                  const int cell_size,
                                  const vector<double> &points,
                                  const vector<int> &tets) {
  
    return write_mesh(os, dim, cell_size, points, tets, true);
}

bool VTUWriter::write_point_cloud(std::ostream &os, const int dim,
                                  const vector<double> &points) {
  
  assert(dim > 1);

  int num_points = points.size() / dim;

  write_header(num_points, 0, os);
  write_points(num_points, points, os, false);
  write_point_data(os);
  os << "<Cells>\n";
  /////////////////////////////////////////////////////////////////////////////
  // List vertex id's i=0, ..., n_vertices associated with each cell c
  os << "<DataArray type=\"Int64\" Name=\"connectivity\" "
        "format=\"ascii\">\n";

  os << "</DataArray>\n";
  os << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" "
        "RangeMin=\""
     <<1e+299 << "\" RangeMax=\"" << -1e+299 << "\">\n";
  os << "</DataArray>\n";
  os << "</Cells>\n";
  write_footer(os);
  clear();
  return true;
}

bool VTUWriter::write_point_cloud(const std::string &path, const int dim,
                                  const vector<double> &points) {
  
  std::ofstream os;
  os.open(path.c_str());
  if (!os.good()) {
    os.close();
    return false;
  }
  write_point_cloud(os, dim, points);
  os.close();
  return true;
}

} // namespace leanvtk
