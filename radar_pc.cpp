#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

//types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Normal;
typedef CGAL::Point_set_3<Point, Normal> Point_set;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

class My_point_property_map
{
  const CGAL::Point_set_3<Point, Normal>& points;
public:
  typedef Point value_type;
  typedef const value_type& reference;
  typedef std::size_t key_type;
  typedef boost::lvalue_property_map_tag category;
  My_point_property_map(const CGAL::Point_set_3<Point, Normal>& pts):points(pts){}
  reference operator[](key_type k) const {return points.point(k);}
  friend reference get(const My_point_property_map& ppmap,key_type i)
  {return ppmap[i];}
};


int main(int argc, char*argv[])
{
    //reading a point set file in points.
    const std::string filename = argc > 1 ? argv[1] : CGAL::data_file_path("points_3/sphere_1k.xyz");

    Point_set point_set;
    if(!CGAL::IO::read_point_set(filename, point_set))
    {
        std::cerr << "Can't read input file " << filename << std::endl;
        return EXIT_FAILURE;
    }

    //estimating normal directions
    const int k = 20; 
    CGAL::pca_estimate_normals<Concurrency_tag>(point_set, k);

    // Orients normals.
    Point_set::iterator unoriented_points_begin = CGAL::mst_orient_normals(point_set, k);

    //delete points with an unoriented normal
    point_set.remove(unoriented_points_begin, point_set.end());


    //performing a k nearest neighbor search for a chosen point
    typedef CGAL::Search_traits_3<Kernel>                                        Traits_base;
    typedef CGAL::Search_traits_adapter<std::size_t,My_point_property_map,Traits_base> Traits;
    typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_neighbor_search;
    typedef K_neighbor_search::Tree                                         Tree;
    typedef Tree::Splitter                                                  Splitter;
    typedef K_neighbor_search::Distance                                     Distance;

    My_point_property_map ppmap(point_set);

    Point query_point = point_set.point(2);

    Tree tree(boost::counting_iterator<std::size_t>(0),
              boost::counting_iterator<std::size_t>(point_set.size()),
              Splitter(),
              Traits(ppmap));

    Distance tr_dist(ppmap);
    K_neighbor_search search(tree, query_point, k,0,true,tr_dist);

    //adding color to neighbors 
    typedef std::array<unsigned char, 3> Color;
    typedef Point_set::Property_map<Color> Color_map; 

    Color black = {{ 0, 0, 0 }};
    bool success = false;
    Color_map color;

    boost::tie (color, success) = point_set.add_property_map<Color> ("color", black);
    assert (success);

    for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++){
    std::cout << " d(q, nearest neighbor)=  "<< tr_dist.inverse_of_transformed_distance(it->second) << " index:" << it->first << std::endl;

    auto ps_it = point_set.begin();
    std::advance(ps_it, it->first);
    Color red = {{255,0,0}};
    color[*ps_it] = red;
  }

    // Write to ply file
    const std::string output_filename = "output.ply";
    if (!CGAL::IO::write_point_set(output_filename, point_set))
    {
        std::cerr << "Can't create output file " << output_filename << std::endl;
        return EXIT_FAILURE;
    }

  return 0;
}