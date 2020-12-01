/*
 *
 * This code shuffles an edgelist according to various null models.
 *
 * Compilation requires the c++11 standard.
 *   Example: g++ -O3 -std=c++11 shuffle_edges.cpp -o shuffle_edges
 *
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    Nov. 2020
 *
 */

// Standard template library
#include <string>  // std::string
// PGL
#include "../src/lccm_rewiring_t.hpp"

int main(int argc, char** argv)
{
  // Filename of the original edgelist.
  std::string edgelist_filename = "example00_edgelist.dat";
  std::string rootname = edgelist_filename.substr(0, edgelist_filename.find("_edgelist.dat"));

  // Number of times each edge should be rewired (on average).
  int average_number_of_rewirings_per_edge = 1000;

  // Shuffles the edges according to the LCCM (no correlation).
  lccm_rewiring_t lccm_none_network(edgelist_filename, lccm_rewiring_t::ec_none);
  lccm_none_network.shuffle_edges(average_number_of_rewirings_per_edge);
  lccm_none_network.write_edgelist(rootname + "_lccm_none_edgelist.dat");

  // Shuffles the edges according to the LCCM and enforcing the degree-degree  correlations.
  lccm_rewiring_t lccm_KxK_network(edgelist_filename, lccm_rewiring_t::ec_KxK);
  lccm_KxK_network.shuffle_edges(average_number_of_rewirings_per_edge);
  lccm_KxK_network.write_edgelist(rootname + "_lccm_KxK_edgelist.dat");

  // Shuffles the edges according to the LCCM and enforcing the degree-layer-degree -layer correlations.
  lccm_rewiring_t lccm_LKxLK_network(edgelist_filename, lccm_rewiring_t::ec_LKxLK);
  lccm_LKxLK_network.shuffle_edges(average_number_of_rewirings_per_edge);
  lccm_LKxLK_network.write_edgelist(rootname + "_lccm_LKxLK_edgelist.dat");

  // Exits the program successfully.
  return 0;
}
