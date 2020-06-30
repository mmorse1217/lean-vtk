#include <cassert>
#include <iostream>
#include <lean_vtk.hpp>
#include <vector>
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
std::string a_library_function(){
    return std::string("a function specified in the source code");
}

namespace polyfem
{
    namespace
    {
        static const int VTK_TETRA = 10;
        static const int VTK_TRIANGLE = 5;
        static const int VTK_QUAD = 9;
        static const int VTK_HEXAHEDRON = 12;
        static const int VTK_POLYGON = 7;

        inline static int VTKTagVolume(const int n_vertices)
        {
            switch (n_vertices) {
                case 4:
                return VTK_TETRA;
                case 8:
                return VTK_HEXAHEDRON;
                default:
                //element type not supported. To add it (http://www.vtk.org/VTK/img/file-formats.pdf)
                cerr << n_vertices << " not supported, " << endl;
                assert(false);
                return -1;
            }
        }

        inline static int VTKTagPlanar(const int n_vertices)
        {
            switch (n_vertices) {
                case 3:
                return VTK_TRIANGLE;
                case 4:
                return VTK_QUAD;
                default:
                //element type not supported. To add it (http://www.vtk.org/VTK/img/file-formats.pdf)
                cerr <<"{} not supported, " << n_vertices << endl;
                assert(false);
                return -1;
            }
        }
    }

    void VTUWriter::write_point_data(std::ostream &os)
    {
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

    void VTUWriter::write_header(const int n_vertices, const int n_elements, std::ostream &os)
    {
        os << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        os << "<UnstructuredGrid>\n";
        os << "<Piece NumberOfPoints=\"" << n_vertices << "\" NumberOfCells=\"" << n_elements << "\">\n";
    }

    void VTUWriter::write_footer(std::ostream &os)
    {
        os << "</Piece>\n";
        os << "</UnstructuredGrid>\n";
        os << "</VTKFile>\n";
    }

    void VTUWriter::write_points(const int num_points, const vector<double>& points, std::ostream &os)
    {
        os << "<Points>\n";
        os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        const int dim = points.size()/num_points;
        assert(double(dim) == double(points.size())/double(num_points));

        //for (int d = 0; d < points.rows(); ++d)
        for (int d = 0; d < num_points; ++d)
        {
            //for (int i = 0; i < points.cols(); ++i)
            for (int i = 0; i < dim; ++i)
            {
                //os << points(d, i);
                int idx = index(dim, d, i); // TODO fix
                os << points.at(idx);
                if (i < dim - 1)
                {
                    os << " ";
                }
            }

            if(!is_volume_)
                os << " 0";

            os << "\n";
        }

        os << "</DataArray>\n";
        os << "</Points>\n";
    }

    void VTUWriter::write_cells(const int n_vertices, const vector<int>& tets, std::ostream &os)
    {
        //const int n_cells = tets.rows();
        const int n_cells = tets.size()/n_vertices;
        os << "<Cells>\n";
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        os << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";

        //const int n_vertices = tets.cols();

        for(int c = 0; c < n_cells; ++c)
        {
            for (int i = 0; i < n_vertices; ++i)
            {
                int idx = index(n_vertices, c, i);
                //const int v_index = tets(c,i);
                const int v_index = tets.at(idx);
                os << v_index;
                if (i < n_vertices - 1) {
                    os << " ";
                }
            }
            os << "\n";
        }

        os << "</DataArray>\n";
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        int min_tag, max_tag;
        if (!is_volume_) {
            min_tag = VTKTagPlanar(n_vertices);
            max_tag = VTKTagPlanar(n_vertices);
        } else
        {
            min_tag = VTKTagVolume(n_vertices);
            max_tag = VTKTagVolume(n_vertices);
        }

        os << "<DataArray type=\"Int8\" Name=\"types\" format=\"ascii\" RangeMin=\"" << min_tag << "\" RangeMax=\"" << max_tag << "\">\n";
        for (int i = 0; i < n_cells; ++i)
        {
            if (is_volume_)
                os << VTKTagVolume(n_vertices) << "\n";
            else
                os << VTKTagPlanar(n_vertices) << "\n";
        }
        os << "</DataArray>\n";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        os << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"" << n_vertices << "\" RangeMax=\"" << n_cells *n_vertices << "\">\n";

        int acc = n_vertices;
        for (int i = 0; i < n_cells; ++i) {
            os << acc << "\n";
            acc += n_vertices;
        }

        os << "</DataArray>\n";
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        os << "</Cells>\n";
    }

    void VTUWriter::clear()
    {
        point_data_.clear();
        cell_data_.clear();
    }

    void VTUWriter::add_field(const std::string &name, const vector<double>& data, const int& dimension)
    {
        using std::abs;

        //Eigen::MatrixXd tmp;
        vector<double> tmp(data.size(), 0.);
        //tmp.resizeLike(data);

        for(long i = 0; i < data.size(); ++i)
            tmp[i] = abs(data[i]) < 1e-16 ? 0 : data[i];

        if(dimension == 1)
            add_scalar_field(name, tmp);
        else
            add_vector_field(name, tmp, dimension);
    }

    void VTUWriter::add_scalar_field(const std::string &name, const vector<double>& data)
    {
        point_data_.push_back(VTKDataNode<double>());
        point_data_.back().initialize(name, "Float64", data);
        current_scalar_point_data_ = name;
    }

    void VTUWriter::add_vector_field(const std::string &name, const vector<double>& data, const int& dimension)
    {
        point_data_.push_back(VTKDataNode<double>());

        //FIXME?
        // if (data.cols() == 2)
        // {
        //  express::BlockEigen::MatrixXd<typename Eigen::MatrixXd::EntryType> data3((data.rows() * 3) / 2, data.columns(), 3, data.columns());
        //  data3.allSet(0);
        //  for (int i = 0; i < data3.nBlockRows(); ++i) {
        //      data3.setBlockAt(i, 0, data.rowRange(i * n_components, (i + 1) * n_components));
        //  }
        //  point_data_.back().initialize(name, "Float64", data3, 3);
        // } else

        point_data_.back().initialize(name, "Float64", data, dimension);
        current_vector_point_data_ = name;
    }

    bool VTUWriter::write_tet_mesh(const std::string &path, const int dim, const int cell_size,
            const vector<double>& points, const vector<int>& tets)
    {
      assert(dim > 1);
      assert(cell_size > 2);

      std::ofstream os;
      os.open(path.c_str());
      if (!os.good()) {
        os.close();
        return false;
      }
        //const static int DIM = 3;
        int num_points = points.size()/dim;
        int num_cells = tets.size()/cell_size;
        cout << "num_points: " << num_points << endl;
        cout << "num_cells: " << num_cells << endl;
        //is_volume_ = points.cols() == 3;
        is_volume_ = dim == 3;

        write_header(num_points, num_cells, os);
        write_points(num_points, points, os);
        write_point_data(os);
        write_cells(cell_size, tets, os);

        write_footer(os);
        os.close();
        clear();
        return true;
    }
}
