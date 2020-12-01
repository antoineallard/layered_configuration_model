/*
 *
 * This class provides the functions to shuffle an edgelist according to the LCCM.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++11 my_code.cpp
 *
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    Dec. 2017
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */

// Standard template library
#include <algorithm> // std::max_element
#include <cmath>     // std::floor
#include <cstdlib>   // std::terminate
#include <ctime>     // std::time(NULL)
#include <fstream>   // std::fstream, std::fstream::in, std::fstream::out
#include <iomanip>   // std::setw
#include <iostream>  // std::cerr, std::endl, std::ws
#include <map>       // std::map
#include <set>       // std::set
#include <sstream>   // std::stringstream
#include <string>    // std::string, std::getline
#include <random>    // std::mt19937, std::uniform_real_distribution
#include <utility>   // std::swap, std::pair, std::make_pair
#include <vector>    // std::vector


class lccm_rewiring_t
{
  public:
    // Custom type to indicate the type of correlation to enforce.
    enum ec_t {ec_none, ec_KxK, ec_LKxLK};
  private:
    // Name and ID conversion.
    std::map<std::string, int> Name2ID;
    std::vector<std::string> ID2Name;
    // Edge correlation to enforce.
    ec_t edge_correlation;
    // Number of vertices.
    int nb_vertices;
    // Number of edges.
    int nb_edges;
    // Core of layers.
    std::vector<int> c;
    // Core of vertices.
    std::vector<int> coreness;
    // Layer of vertices.
    std::vector<int> layer;
    // Degree of vertices.
    std::vector<int> degree;
    // Edgelist.
    std::set< std::pair<int, int> > list_of_unique_edges;
    std::vector< std::vector<int> > edgelist;
    // List of edge id based on the degree of its vertices.
    std::vector<int> nb_edges_per_vertex_class_KxK;
    std::vector<int> vertex_classes_KxK;
    std::vector< std::vector<int> > edges_per_vertex_class_KxK;
    // List of edge id based on the degree/layer of its vertices.
    std::vector<int> nb_edges_per_vertex_class_LKxLK;
    std::vector< std::pair<int, int> > vertex_classes_LKxLK;
    std::vector< std::vector<int> > edges_per_vertex_class_LKxLK;
    // List of edges to account for multiedges.
    std::vector< std::set<int> > adjacency_list;
    // Conditions ensuring the position in the OD.
    std::vector< std::vector<int> > stubs;
    // Random number generator
    std::mt19937 engine;
    std::uniform_real_distribution<double> uniform_01;
    std::discrete_distribution<int> get_vertex_class;
    // Loads the edgelist.
    void load_edgelist(std::string edgelist_filename);
    // Utilities.
    void build_auxiliary_objects();
    bool check_local_rules(int v1, int l1, int nb_stubs_0, int nb_stubs_1);
    bool is_layer_preserved(int v1, int l1, int l2, int l3);
    bool existing_edge(int v1, int v2);
    void onion_decomposition();
    void update_stubs(int v1, int l1, int old_neigh_layer, int new_neigh_layer);
    // Shuffles without enforcing correlations.
    void shuffle_edges_ec_none(double nb_times_nb_edges);
    // Shuffles while enforcing degree-degree correlations.
    void shuffle_edges_ec_KxK(double nb_times_nb_edges);
    // Shuffles while enforcing degree-layer-degree-layer correlations.
    void shuffle_edges_ec_LKxLK(double nb_times_nb_edges);
  public:
    // Constructor.
    lccm_rewiring_t(std::string edgelist_filename, ec_t = ec_LKxLK);
    // Rewires the edges acording to the LCCM.
    void shuffle_edges(double nb_times_nb_edges = 10);
    // Writes the shuffled edgelist into a file.
    void write_edgelist(std::string filename, int width = 15);
};


// ================================================================================================
// ================================================================================================
lccm_rewiring_t::lccm_rewiring_t(std::string edgelist_filename, ec_t ec)
{
  // Initializes the random number generator.
  engine.seed(std::time(NULL));
  // Edge correlations to enforce.
  edge_correlation = ec;
  // Loads the data.
  load_edgelist(edgelist_filename);
  // Builds auxiliary objects.
  build_auxiliary_objects();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void lccm_rewiring_t::build_auxiliary_objects()
{
  // Variables.
  int v1, v2, k1, k2, l1, l2, n, m;
  int e = 0;
  // Number of edges.
  nb_edges = list_of_unique_edges.size();
  // Resizes the containers.
  stubs.resize(nb_vertices, std::vector<int>(2, 0));
  degree.resize(nb_vertices, 0);
  adjacency_list.resize(nb_vertices);
  edgelist.resize(nb_edges, std::vector<int>(2, 0));
  // Fills the containers.
  for(auto el : list_of_unique_edges)
  {
    // Identifies the vertices.
    v1 = el.first;
    v2 = el.second;
    // Fills the list of edges.
    edgelist[e][0] = v1;
    edgelist[e][1] = v2;
    ++e;
    // Fills the adjacency list.
    adjacency_list[v1].insert(v2);
    adjacency_list[v2].insert(v1);
    // Updates the degree.
    degree[v1] += 1;
    degree[v2] += 1;
  }
  // Clears unecessary objects.
  list_of_unique_edges.clear();
  // Extract the onion decomposition.
  onion_decomposition();
  // Counts the number of stubs to ensure that the OD is preserved.
  for(e = 0; e<nb_edges; ++e)
  {
    // Identifies the vertices and their layers.
    v1 = edgelist[e][0];
    v2 = edgelist[e][1];
    l1 = layer[v1];
    l2 = layer[v2];
    // Registers the conditions.
    if(l2 >= (l1 - 1))
      stubs[v1][0] += 1;
    if(l2 >= l1)
      stubs[v1][1] += 1;
    if(l1 >= (l2 - 1))
      stubs[v2][0] += 1;
    if(l1 >= l2)
      stubs[v2][1] += 1;
  }
  // Fills the lists of edges based on the degree of vertices.
  if(edge_correlation == ec_KxK)
  {
    // Map.
    std::map< int, std::multiset<int> > tmp_edges_per_vertex_class;
    for(e = 0; e<nb_edges; ++e)
    {
      // Gets the degree of the vertices.
      k1 = degree[edgelist[e][0]];
      k2 = degree[edgelist[e][1]];
      // Inserts the id of the edge.
      tmp_edges_per_vertex_class[k1].insert(e);
      tmp_edges_per_vertex_class[k2].insert(e);
    }
    // Builds the objects to select edges during rewiring.
    n = 0;
    vertex_classes_KxK.resize(tmp_edges_per_vertex_class.size());
    nb_edges_per_vertex_class_KxK.resize(tmp_edges_per_vertex_class.size());
    edges_per_vertex_class_KxK.resize(tmp_edges_per_vertex_class.size());
    for(auto el1 : tmp_edges_per_vertex_class)
    {
      vertex_classes_KxK[n] = el1.first;
      nb_edges_per_vertex_class_KxK[n] = el1.second.size();
      edges_per_vertex_class_KxK[n].resize(el1.second.size());
      m = 0;
      for(auto el2 : el1.second)
      {
        edges_per_vertex_class_KxK[n][m] = el2;
        ++m;
      }
      ++n;
    }
    // Sets the discrete distribution.
    std::discrete_distribution<int> tmp_get_vertex_class(nb_edges_per_vertex_class_KxK.begin(), nb_edges_per_vertex_class_KxK.end());
    get_vertex_class.param(tmp_get_vertex_class.param());
  }
  // Fills the lists of edges based on the degree/layer of vertices.
  if(edge_correlation == ec_LKxLK)
  {
    // Pairs.
    std::pair<int, int> p1, p2;
    // Map.
    std::map< std::pair<int, int>, std::multiset<int> > tmp_edges_per_vertex_class;
    // Surveys edges in terms of the class of their vertices.
    for(e = 0; e<nb_edges; ++e)
    {
      // Gets the degree of the vertices.
      k1 = degree[edgelist[e][0]];
      k2 = degree[edgelist[e][1]];
      // Gets the layer of the vertices.
      l1 = layer[edgelist[e][0]];
      l2 = layer[edgelist[e][1]];
      // Builds the pair objects.
      p1 = std::make_pair(l1, k1);
      p2 = std::make_pair(l2, k2);
      // Inserts the id of the edge.
      tmp_edges_per_vertex_class[p1].insert(e);
      tmp_edges_per_vertex_class[p2].insert(e);
    }
    // Builds the objects to select edges during rewiring.
    n = 0;
    vertex_classes_LKxLK.resize(tmp_edges_per_vertex_class.size());
    nb_edges_per_vertex_class_LKxLK.resize(tmp_edges_per_vertex_class.size());
    edges_per_vertex_class_LKxLK.resize(tmp_edges_per_vertex_class.size());
    for(auto el1 : tmp_edges_per_vertex_class)
    {
      vertex_classes_LKxLK[n] = el1.first;
      nb_edges_per_vertex_class_LKxLK[n] = el1.second.size();
      edges_per_vertex_class_LKxLK[n].resize(el1.second.size());
      m = 0;
      for(auto el2 : el1.second)
      {
        edges_per_vertex_class_LKxLK[n][m] = el2;
        ++m;
      }
      ++n;
    }
    // Sets the discrete distribution.
    std::discrete_distribution<int> tmp_get_vertex_class(nb_edges_per_vertex_class_LKxLK.begin(), nb_edges_per_vertex_class_LKxLK.end());
    get_vertex_class.param(tmp_get_vertex_class.param());
  }
}


// ================================================================================================
// ================================================================================================
bool lccm_rewiring_t::check_local_rules(int v1, int l1, int nb_stubs_0, int nb_stubs_1)
{
  // Checks if in the first layer of the core.
  if(c[l1 - 1] != c[l1])
  {
    return (nb_stubs_1 == coreness[v1]);
  }
  else
  {
    int k = coreness[v1];
    return ((nb_stubs_0 >= (k + 1)) && (nb_stubs_1 <= k));
  }
}


// ================================================================================================
// ================================================================================================
bool lccm_rewiring_t::existing_edge(int v1, int v2)
{
  // Assumes that the adjacency list is symmetric.
  if(adjacency_list[v1].find(v2) != adjacency_list[v1].end())
    return true;
  else
    return false;
}


// ================================================================================================
// ================================================================================================
bool lccm_rewiring_t::is_layer_preserved(int v1, int l1, int old_neigh_layer, int new_neigh_layer)
{
  if(old_neigh_layer >= l1)
  {
    if(new_neigh_layer >= l1)
    {
      return true;
    }
    else if(new_neigh_layer < (l1 - 1))
    {
      return check_local_rules(v1, l1, stubs[v1][0] - 1, stubs[v1][1] - 1);
    }
    else if(new_neigh_layer == (l1 - 1))
    {
      return check_local_rules(v1, l1, stubs[v1][0], stubs[v1][1] - 1);
    }
    else
    {
      std::cerr << "ERROR: This should never happen (#1)." << std::endl;
      std::terminate();
    }
  }
  else if(old_neigh_layer < (l1 - 1))
  {
    if(new_neigh_layer >= l1)
    {
      return check_local_rules(v1, l1, stubs[v1][0] + 1, stubs[v1][1] + 1);
    }
    else if(new_neigh_layer < (l1 - 1))
    {
      return true;
    }
    else if(new_neigh_layer == (l1 - 1))
    {
      return check_local_rules(v1, l1, stubs[v1][0] + 1, stubs[v1][1]);
    }
    else
    {
      std::cerr << "ERROR: This should never happen (#2)." << std::endl;
      std::terminate();
    }
  }
  else if(old_neigh_layer == (l1 - 1))
  {
    if(new_neigh_layer >= l1)
    {
      return check_local_rules(v1, l1, stubs[v1][0], stubs[v1][1] + 1);
    }
    else if(new_neigh_layer < (l1 - 1))
    {
      return check_local_rules(v1, l1, stubs[v1][0] - 1, stubs[v1][1]);
    }
    else if(new_neigh_layer == (l1 - 1))
    {
      return true;
    }
    else
    {
      std::cerr << "ERROR: This should never happen (#3)." << std::endl;
      std::terminate();
    }
  }
  else
  {
    std::cerr << "ERROR: This should never happen (#4)." << std::endl;
    std::terminate();
  }
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void lccm_rewiring_t::load_edgelist(std::string edgelist_filename)
{
  // Stream objects.
  std::fstream edgelist_file;
  std::stringstream one_line;
  // Variables.
  int v1, v2;
  // Initializes the number of vertices.
  nb_vertices = 0;
  // String objects.
  std::string full_line, name1_str, name2_str;
  // Iterators.
  std::map< std::string, int >::iterator name_it;
  // Opens the stream and terminates if the operation did not succeed.
  edgelist_file.open(edgelist_filename.c_str(), std::fstream::in);
  if( !edgelist_file.is_open() )
  {
    std::cerr << "ERROR: Could not open file: " << edgelist_filename << "." << std::endl;
    std::terminate();
  }
  // Reads the edgelist file line by line.
  while( !edgelist_file.eof() )
  {
    // Reads a line of the file.
    std::getline(edgelist_file, full_line);
    edgelist_file >> std::ws;
    one_line.str(full_line);
    one_line >> std::ws;
    one_line >> name1_str >> std::ws;
    // Skips a line of comment.
    if(name1_str == "#")
    {
      one_line.clear();
      continue;
    }
    one_line >> name2_str >> std::ws;
    one_line.clear();
    // Ignores self-loops.
    if(name1_str != name2_str)
    {
      // Is name1 new?
      name_it = Name2ID.find(name1_str);
      if(name_it == Name2ID.end())
      {
        Name2ID[name1_str] = nb_vertices;
        v1 = nb_vertices;
        ++nb_vertices;
      }
      else
      {
        v1 = name_it->second;
      }
      // Is name2 new?
      name_it = Name2ID.find(name2_str);
      if(name_it == Name2ID.end())
      {
        Name2ID[name2_str] = nb_vertices;
        v2 = nb_vertices;
        ++nb_vertices;
      }
      else
      {
        v2 = name_it->second;
      }
      // Adds the edge (multiedges are ignored).
      if(v1 > v2)
      {
        std::swap(v1, v2);
      }
      list_of_unique_edges.insert(std::make_pair(v1, v2));
    }
  }
  // Closes the stream.
  edgelist_file.close();
}


// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void lccm_rewiring_t::onion_decomposition()
{
  // Resizes the containers.
  coreness.resize(nb_vertices);
  layer.resize(nb_vertices);
  // Builds two lists (std::vector, std::set) of the degree of the vertices.
  std::vector<int> DegreeVec(nb_vertices);
  std::set<std::pair<int, int> > DegreeSet;
  for(int v(0); v<nb_vertices; ++v)
  {
    DegreeSet.insert(std::make_pair(degree[v], v));
    DegreeVec[v] = degree[v];
  }
  // Determines the coreness and the layer based on the modified algorithm of Batagelj and
  //   Zaversnik by Hébert-Dufresne, Grochow and Allard.
  int v1, v2, d1, d2, c_max;
  int id_layer = 0;
  std::set< std::pair<int, int> > LayerSet;
  std::set< std::pair<int, int> >::iterator m_it;
  while(DegreeSet.size() > 0)
  {
    // Increases the layer id.
    id_layer += 1;
    // Populates the set containing the vertices belonging to the same layer.
    m_it = DegreeSet.begin();
    d1 = m_it->first;
    // Sets the coreness and the layer the vertices with the same degree.
    while(m_it->first == d1 && m_it != DegreeSet.end())
    {
      // Sets the coreness and the layer.
      v1 = m_it->second;
      coreness[v1] = d1;
      layer[v1] = id_layer;
      // Looks at the next vertex.
      ++m_it;
    }
    // Adds the vertices of the layer to the set.
    LayerSet.insert(DegreeSet.begin(), m_it);
    // Removes the vertices of the current layer.
    DegreeSet.erase(DegreeSet.begin(), m_it);
    // Modifies the "effective" degree of the neighbors of the vertices in the layer.
    while(LayerSet.size() > 0)
    {
      // Gets information about the next vertex of the layer.
      v1 = LayerSet.begin()->second;
      // Reduces the "effective" degree of its neighbours.
      // for(int n(0), nn(adjacency_list[v1].size()); n<nn; ++n)
      for(auto el : adjacency_list[v1])
      {
        // Identifies the neighbor.
        v2 = el;
        // v2 = adjacency_list[v1][n];
        d2 = DegreeVec[v2];
        // Finds the neighbor in the list "effective" degrees.
        m_it = DegreeSet.find(std::make_pair(d2, v2));
        if(m_it != DegreeSet.end())
        {
          if(d2 > d1)
          {
            DegreeVec[v2] = d2 - 1;
            DegreeSet.erase(m_it);
            DegreeSet.insert(std::make_pair(d2 - 1, v2));
          }
        }
      }
      // Removes the vertices from the LayerSet.
      LayerSet.erase(LayerSet.begin());
    }
  }
  // Builds the layer-to-coreness vector.
  c_max = *std::max_element(layer.begin(), layer.end());
  c.resize(c_max + 1);
  for(int v(0); v<nb_vertices; ++v)
  {
    c[layer[v]] = coreness[v];
  }
}


// ================================================================================================
// ================================================================================================
void lccm_rewiring_t::shuffle_edges(double nb_times_nb_edges)
{
  if(edge_correlation == ec_LKxLK)
    shuffle_edges_ec_LKxLK(nb_times_nb_edges);
  else if(edge_correlation == ec_KxK)
    shuffle_edges_ec_KxK(nb_times_nb_edges);
  else if(edge_correlation == ec_none)
    shuffle_edges_ec_none(nb_times_nb_edges);
  else
    std::cerr << "ERROR: Invalid option for edge correlation." << std::endl;
}

// ================================================================================================
// ================================================================================================
void lccm_rewiring_t::shuffle_edges_ec_none(double nb_times_nb_edges)
{
  // Counting variables.
  int nb_swaps = 0;
  int nb_attempts = 0;
  int nb_remaining_attempts = nb_edges * (1 + uniform_01(engine));
  int nb_swaps_to_achieve = nb_times_nb_edges * nb_edges;
  // Variables.
  int e12, e34, v1, v2, v3, v4, l1, l2, l3, l4;
  // Boolean indicators.
  std::vector<bool> visited(nb_edges, false);
  bool keep_going = true;
  while(keep_going)
  {
    ++nb_attempts;
    if(nb_attempts%1000==0)
    {
      std::cout << "\r" << std::setw(6) << nb_swaps << " swaps out of " << std::setw(6) << nb_swaps_to_achieve << " after " << std::setw(9) << nb_attempts << " attempts."/* << std::endl*/;
      std::cout.flush();
    }
    // Checks if enough swaps/attempts have been made.
    if(nb_swaps > nb_swaps_to_achieve)
    {
      --nb_remaining_attempts;
      if(nb_remaining_attempts == 0)
      {
        std::cout << "\r" << std::setw(6) << nb_swaps << " swaps out of " << std::setw(6) << nb_swaps_to_achieve << " after " << std::setw(9) << nb_attempts << " attempts." << std::endl;
        std::cout.flush();
        keep_going = false;
      }
    }
    // Chooses the first edge.
    e12 = static_cast<int>(std::floor(nb_edges*uniform_01(engine)));
    // Gets information about the nodes involved in the first edge.
    v1 = edgelist[e12][0];
    v2 = edgelist[e12][1];
    l1 = layer[v1];
    l2 = layer[v2];
    // Chooses the second edge without any discrimination.
    e34 = std::floor(nb_edges*uniform_01(engine));
    while(e12 == e34)
      e34 = std::floor(nb_edges*uniform_01(engine));
    // Gets information about the nodes involved in the second edge.
    v3 = edgelist[e34][0];
    v4 = edgelist[e34][1];
    l3 = layer[v3];
    l4 = layer[v4];
    // Determines which new edge will be atempted.
    if(uniform_01(engine) < 0.5)
    {
      std::swap(v3, v4);
      std::swap(l3, l4);
    }
    // Checks wether the swap preserves the layer of nodes.
    if(!is_layer_preserved(v1, l1, l2, l3))
      continue;
    if(!is_layer_preserved(v2, l2, l1, l4))
      continue;
    if(!is_layer_preserved(v3, l3, l4, l1))
      continue;
    if(!is_layer_preserved(v4, l4, l3, l2))
      continue;
    // Checks for self-links.
    if(v1 == v3)
      continue;
    if(v2 == v4)
      continue;
    // Checks for multi-links.
    if(existing_edge(v1, v3))
      continue;
    if(existing_edge(v2, v4))
      continue;
    // Effectuates the swap since it is valid.
    edgelist[e12][0] = v1;
    edgelist[e12][1] = v3;
    edgelist[e34][0] = v2;
    edgelist[e34][1] = v4;
    adjacency_list[v1].erase(v2);
    adjacency_list[v2].erase(v1);
    adjacency_list[v3].erase(v4);
    adjacency_list[v4].erase(v3);
    adjacency_list[v1].insert(v3);
    adjacency_list[v3].insert(v1);
    adjacency_list[v2].insert(v4);
    adjacency_list[v4].insert(v2);
    update_stubs(v1, l1, l2, l3);
    update_stubs(v2, l2, l1, l4);
    update_stubs(v3, l3, l4, l1);
    update_stubs(v4, l4, l3, l2);
    visited[e12] = true;
    visited[e34] = true;
    ++nb_swaps;
  }
  // Checks whether all edges have been swapped at least once.
  bool all_visited = false;
  for(int e(0); e<nb_edges; ++e)
  {
    if(!visited[e])
    {
      std::cerr << "WARNING: Not all vertices have been swapped at least once." << std::endl;
      break;
    }
  }
}


// ================================================================================================
// ================================================================================================
void lccm_rewiring_t::shuffle_edges_ec_LKxLK(double nb_times_nb_edges)
{
  // Counting variables.
  int nb_attempts = 0;
  int nb_swaps = 0;
  int nb_remaining_attempts = nb_edges * (1 + uniform_01(engine));
  int nb_swaps_to_achieve = nb_times_nb_edges * nb_edges;
  // Variables.
  int e12, n_e12, e34, n_e34, v1, v2, v3, v4, k1, k2, k3, k4, l1, l2, l3, l4, vc, l_vc, k_vc;
  // Boolean indicators.
  std::vector<bool> visited(nb_edges, false);
  bool keep_going = true;
  while(keep_going)
  {
    ++nb_attempts;
    if(nb_attempts%1000==0)
    {
      std::cout << "\r" << std::setw(6) << nb_swaps << " swaps out of " << std::setw(6) << nb_swaps_to_achieve << " after " << std::setw(9) << nb_attempts << " attempts."/* << std::endl*/;
      std::cout.flush();
    }
    // Checks if enough swaps/attempts have been made.
    if(nb_swaps > nb_swaps_to_achieve)
    {
      --nb_remaining_attempts;
      if(nb_remaining_attempts == 0)
      {
        std::cout << "\r" << std::setw(6) << nb_swaps << " swaps out of " << std::setw(6) << nb_swaps_to_achieve << " after " << std::setw(9) << nb_attempts << " attempts." << std::endl;
        std::cout.flush();
        keep_going = false;
      }
    }
    // Selects the edge class.
    vc = get_vertex_class(engine);
    if(nb_edges_per_vertex_class_LKxLK[vc] == 1)
      continue;
    // Extracts the info about the vertex class.
    l_vc = vertex_classes_LKxLK[vc].first;
    k_vc = vertex_classes_LKxLK[vc].second;
    // Chooses the edges.
    n_e12 = std::floor(nb_edges_per_vertex_class_LKxLK[vc]*uniform_01(engine));
    n_e34 = std::floor(nb_edges_per_vertex_class_LKxLK[vc]*uniform_01(engine));
    while(n_e12 == n_e34)
    {
      n_e34 = std::floor(nb_edges_per_vertex_class_LKxLK[vc]*uniform_01(engine));
    }
    // Gets information about the nodes involved in the first edge.
    e12 = edges_per_vertex_class_LKxLK[vc][n_e12];
    v1 = edgelist[e12][0];
    v2 = edgelist[e12][1];
    k1 = degree[v1];
    k2 = degree[v2];
    l1 = layer[v1];
    l2 = layer[v2];
    // Gets information about the nodes involved in the second edge.
    e34 = edges_per_vertex_class_LKxLK[vc][n_e34];
    v3 = edgelist[e34][0];
    v4 = edgelist[e34][1];
    k3 = degree[v3];
    k4 = degree[v4];
    l3 = layer[v3];
    l4 = layer[v4];
    // If necessary, swaps the labels to make sure that the right vertices are swapped.
    if(!(k1 == k_vc && l1 == l_vc))
    {
      std::swap(v1, v2);
      // std::swap(k1, k2);
      std::swap(l1, l2);
    }
    if(!(k3 == k_vc && l3 == l_vc))
    {
      std::swap(v3, v4);
      // std::swap(k3, k4);
      std::swap(l3, l4);
    }
    // Checks wether the swap preserves the layer of nodes.
    if(!is_layer_preserved(v1, l1, l2, l4))
      continue;
    if(!is_layer_preserved(v2, l2, l1, l3))
      continue;
    if(!is_layer_preserved(v3, l3, l4, l2))
      continue;
    if(!is_layer_preserved(v4, l4, l3, l1))
      continue;
    // Checks for self-links.
    if(v1 == v4)
      continue;
    if(v2 == v3)
      continue;
    // Checks for multi-links.
    if(existing_edge(v1, v4))
      continue;
    if(existing_edge(v2, v3))
      continue;
    // Effectuates the swap since it is valid.
    //   (swaps the id of the edges to avoid modifying the edges_per_vertex_class object)
    edgelist[e34][0] = v1;
    edgelist[e34][1] = v4;
    edgelist[e12][0] = v2;
    edgelist[e12][1] = v3;
    adjacency_list[v1].erase(v2);
    adjacency_list[v2].erase(v1);
    adjacency_list[v3].erase(v4);
    adjacency_list[v4].erase(v3);
    adjacency_list[v1].insert(v4);
    adjacency_list[v2].insert(v3);
    adjacency_list[v3].insert(v2);
    adjacency_list[v4].insert(v1);
    update_stubs(v1, l1, l2, l4);
    update_stubs(v2, l2, l1, l3);
    update_stubs(v3, l3, l4, l2);
    update_stubs(v4, l4, l3, l1);
    visited[e12] = true;
    visited[e34] = true;
    ++nb_swaps;
  }
  // Checks whether all edges have been swapped at least once.
  bool all_visited = false;
  for(int e(0); e<nb_edges; ++e)
  {
    if(!visited[e])
    {
      std::cerr << "WARNING: Not all vertices have been swapped at least once." << std::endl;
      break;
    }
  }
}


// ================================================================================================
// ================================================================================================
void lccm_rewiring_t::shuffle_edges_ec_KxK(double nb_times_nb_edges)
{
  // Counting variables.
  int nb_attempts = 0;
  int nb_swaps = 0;
  int nb_remaining_attempts = nb_edges * (1 + uniform_01(engine));
  int nb_swaps_to_achieve = nb_times_nb_edges * nb_edges;
  // Variables.
  int e12, n_e12, e34, n_e34, v1, v2, v3, v4, k1, k2, k3, k4, l1, l2, l3, l4, vc, l_vc, k_vc;
  // Boolean indicators.
  std::vector<bool> visited(nb_edges, false);
  bool keep_going = true;
  while(keep_going)
  {
    ++nb_attempts;
    if(nb_attempts%1000==0)
    {
      std::cout << "\r" << std::setw(6) << nb_swaps << " swaps out of " << std::setw(6) << nb_swaps_to_achieve << " after " << std::setw(9) << nb_attempts << " attempts."/* << std::endl*/;
      std::cout.flush();
    }
    // Checks if enough swaps/attempts have been made.
    if(nb_swaps > nb_swaps_to_achieve)
    {
      --nb_remaining_attempts;
      if(nb_remaining_attempts == 0)
      {
        std::cout << "\r" << std::setw(6) << nb_swaps << " swaps out of " << std::setw(6) << nb_swaps_to_achieve << " after " << std::setw(9) << nb_attempts << " attempts." << std::endl;
        std::cout.flush();
        keep_going = false;
      }
    }
    // Selects the edge class.
    vc = get_vertex_class(engine);
    if(nb_edges_per_vertex_class_KxK[vc] == 1)
      continue;
    // Extracts the info about the vertex class.
    k_vc = vertex_classes_KxK[vc];
    // Chooses the edges.
    n_e12 = std::floor(nb_edges_per_vertex_class_KxK[vc]*uniform_01(engine));
    n_e34 = std::floor(nb_edges_per_vertex_class_KxK[vc]*uniform_01(engine));
    while(n_e12 == n_e34)
    {
      n_e34 = std::floor(nb_edges_per_vertex_class_KxK[vc]*uniform_01(engine));
    }
    // Gets information about the nodes involved in the first edge.
    e12 = edges_per_vertex_class_KxK[vc][n_e12];
    v1 = edgelist[e12][0];
    v2 = edgelist[e12][1];
    k1 = degree[v1];
    k2 = degree[v2];
    l1 = layer[v1];
    l2 = layer[v2];
    // Gets information about the nodes involved in the second edge.
    e34 = edges_per_vertex_class_KxK[vc][n_e34];
    v3 = edgelist[e34][0];
    v4 = edgelist[e34][1];
    k3 = degree[v3];
    k4 = degree[v4];
    l3 = layer[v3];
    l4 = layer[v4];
    // If necessary, swaps the labels to make sure that the right vertices are swapped.
    if(k1 != k_vc)
    {
      std::swap(v1, v2);
      // std::swap(k1, k2);
      std::swap(l1, l2);
    }
    if(k3 != k_vc)
    {
      std::swap(v3, v4);
      // std::swap(k3, k4);
      std::swap(l3, l4);
    }
    // Checks wether the swap preserves the layer of nodes.
    if(!is_layer_preserved(v1, l1, l2, l4))
      continue;
    if(!is_layer_preserved(v2, l2, l1, l3))
      continue;
    if(!is_layer_preserved(v3, l3, l4, l2))
      continue;
    if(!is_layer_preserved(v4, l4, l3, l1))
      continue;
    // Checks for self-links.
    if(v1 == v4)
      continue;
    if(v2 == v3)
      continue;
    // Checks for multi-links.
    if(existing_edge(v1, v4))
      continue;
    if(existing_edge(v2, v3))
      continue;
    // Effectuates the swap since it is valid.
    //   (swaps the id of the edges to avoid modifying the edges_per_vertex_class object)
    edgelist[e34][0] = v1;
    edgelist[e34][1] = v4;
    edgelist[e12][0] = v2;
    edgelist[e12][1] = v3;
    adjacency_list[v1].erase(v2);
    adjacency_list[v2].erase(v1);
    adjacency_list[v3].erase(v4);
    adjacency_list[v4].erase(v3);
    adjacency_list[v1].insert(v4);
    adjacency_list[v2].insert(v3);
    adjacency_list[v3].insert(v2);
    adjacency_list[v4].insert(v1);
    update_stubs(v1, l1, l2, l4);
    update_stubs(v2, l2, l1, l3);
    update_stubs(v3, l3, l4, l2);
    update_stubs(v4, l4, l3, l1);
    visited[e12] = true;
    visited[e34] = true;
    ++nb_swaps;
  }
  // Checks whether all edges have been swapped at least once.
  bool all_visited = false;
  for(int e(0); e<nb_edges; ++e)
  {
    if(!visited[e])
    {
      std::cerr << "WARNING: Not all vertices have been swapped at least once." << std::endl;
      break;
    }
  }
}


// ================================================================================================
// ================================================================================================
void lccm_rewiring_t::update_stubs(int v1, int l1, int old_neigh_layer, int new_neigh_layer)
{
  if(old_neigh_layer >= l1)
  {
    if(new_neigh_layer >= l1)
    {
      // Nothing to do.
    }
    else if(new_neigh_layer < (l1 - 1))
    {
      stubs[v1][0] -= 1;
      stubs[v1][1] -= 1;
    }
    else if(new_neigh_layer == (l1 - 1))
    {
      stubs[v1][1] -= 1;
    }
    else
    {
      std::cerr << "ERROR: This should never happen (#5)." << std::endl;
      std::terminate();
    }
  }
  else if(old_neigh_layer < (l1 - 1))
  {
    if(new_neigh_layer >= l1)
    {
      stubs[v1][0] += 1;
      stubs[v1][1] += 1;
    }
    else if(new_neigh_layer < (l1 - 1))
    {
      // Nothing to do.
    }
    else if(new_neigh_layer == (l1 - 1))
    {
      stubs[v1][0] += 1;
    }
    else
    {
      std::cerr << "ERROR: This should never happen (#6)." << std::endl;
      std::terminate();
    }
  }
  else if(old_neigh_layer == (l1 - 1))
  {
    if(new_neigh_layer >= l1)
    {
      stubs[v1][1] += 1;
    }
    else if(new_neigh_layer < (l1 - 1))
    {
      stubs[v1][0] -= 1;
    }
    else if(new_neigh_layer == (l1 - 1))
    {
      // Nothing to do.
    }
    else
    {
      std::cerr << "ERROR: This should never happen (#7)." << std::endl;
      std::terminate();
    }
  }
  else
  {
    std::cerr << "ERROR: This should never happen (#8)." << std::endl;
    std::terminate();
  }
}


// =================================================================================================
// =================================================================================================
void lccm_rewiring_t::write_edgelist(std::string filename, int width)
{
  // Builds the ID2Name vector.
  ID2Name.resize(nb_vertices);
  for(auto el : Name2ID)
  {
    ID2Name[el.second] = el.first;
  }
  // Variables.
  int v1, v2;
  // Stream accessing the file.
  std::fstream output_file(filename.c_str(), std::fstream::out);
  if(!output_file.is_open())
  {
    std::cerr << "ERROR: Could not open file: " << filename << "." << std::endl;
    std::terminate();
  }
  // Writes the shuffled edgelist into the file.
  output_file << "#" << std::setw(width - 1) << "Vertex1" << " ";
  output_file        << std::setw(width)     << "Vertex2" << " ";
  output_file << std::endl;
  for(int e(0); e<nb_edges; ++e)
  {
    v1 = edgelist[e][0];
    v2 = edgelist[e][1];
    if(v1 > v2)
    {
      std::swap(v1, v2);
    }
    output_file << std::setw(width) << ID2Name[v1] << " ";
    output_file << std::setw(width) << ID2Name[v2] << " ";
    output_file << std::endl;
  }
  // Closes the stream.
  output_file.close();
}
