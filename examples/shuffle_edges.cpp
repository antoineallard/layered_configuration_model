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
  std::string edgelist_filename = "example04_edgelist.dat";

  // Shuffles the edges according to the LCCM (no correlation).
  lccm_rewiring_t lccm_none_network(edgelist_filename, lccm_rewiring_t::ec_none);
  lccm_none_network.shuffle_edges(1000);
  lccm_none_network.write_edgelist("example04_lccm_none_edgelist.dat");

  // Shuffles the edges according to the LCCM and enforcing the degree-degree  correlations.
  lccm_rewiring_t lccm_KxK_network(edgelist_filename, lccm_rewiring_t::ec_KxK);
  lccm_KxK_network.shuffle_edges(1000);
  lccm_KxK_network.write_edgelist("example04_lccm_KxK_edgelist.dat");

  // Shuffles the edges according to the LCCM and enforcing the degree-layer-degree -layer correlations.
  lccm_rewiring_t lccm_LKxLK_network(edgelist_filename, lccm_rewiring_t::ec_LKxLK);
  lccm_LKxLK_network.shuffle_edges(1000);
  lccm_LKxLK_network.write_edgelist("example04_lccm_LKxLK_edgelist.dat");

  // Exits the program successfully.
  return 0;
}
