/*
//@HEADER
// *****************************************************************************
//
// PuLP: Multi-Objective Multi-Constraint Partitioning Using Label Propagation
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  George M. Slota (gmslota@sandia.gov)
//                      Siva Rajamanickam (srajama@sandia.gov)
//
// *****************************************************************************
//@HEADER
*/

#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <stdint.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>

#include "pulp.h"

#include "rand.cpp"
#include "init_nonrandom.cpp"
#include "label_prop.cpp"
#include "label_balance_verts.cpp"
#include "label_balance_edges.cpp"
#include "label_balance_edges_maxcut.cpp"

int seed;

int* connectivity_bfs(pulp_graph_t& g, int num_parts, int* parts)
{
  int* conn = new int[g.n]; // connectivity assignments  
  for (int v = 0; v < g.n; ++v)
    conn[v] = -1;

  // number of connected components per part
  int* part_conns = (int*)malloc(num_parts*sizeof(int));
  for (int p = 0; p < num_parts; ++p)
    part_conns[p] = 0;

  // Store the size of each component at the root vertex
  int* conn_sizes = (int*) malloc(g.n*sizeof(int));
  for (int v = 0; v < g.n; ++v) {
    conn_sizes[v] = -1;
  }

  int* queue = (int*)malloc(g.n*sizeof(int));
  int* next_queue = (int*)malloc(g.n*sizeof(int));
  int queue_size = 0;
  int next_size = 0;
  int conn_size = 0;

  bool* visited = (bool*)malloc(g.n*sizeof(int));
  for (int v = 0; v < g.n; ++v)
    visited[v] = false;

  int num_comps = 0;

  for (int v = 0; v < g.n; ++v) {
    if (!visited[v]) {
      visited[v] = true;
      conn[v] = v;
      queue[0] = v;
      queue_size = 1;
      next_size = 0;
      conn_size = 1;

      while (queue_size) {
        for (int i = 0; i < queue_size; ++i) {
          int vert = queue[i];

          for (int j = 0; j < out_degree(g, vert); ++j) {
	    int adj = (out_vertices(g, vert))[j];
            if (!visited[adj] && (parts[adj]==parts[v])) {
              visited[adj] = true;
              next_queue[next_size++] = adj;
              conn[adj] = v;
	      ++conn_size;
            }
          }
        }

        int* temp = queue;
        queue = next_queue;
        next_queue = temp;
        queue_size = next_size;
        next_size = 0;
      }

      conn_sizes[v] = conn_size;
      ++part_conns[parts[v]];
      ++num_comps;
    }
  }

  int** part_conn_sizes = (int**) malloc(num_parts*sizeof(int*));
  for (int p = 0; p < num_parts; ++p) {
    part_conn_sizes[p] = (int*) malloc(part_conns[p]*sizeof(int));
  }

  int* part_conn_iter = (int*) malloc(num_parts*sizeof(int));
  for (int p = 0; p < num_parts; ++p) {
    part_conn_iter[p] = 0;
  }

  for (int v = 0; v < g.n; ++v) {
    if (conn[v] == v) {
      int p = parts[v];
      part_conn_sizes[p][part_conn_iter[p]] = conn_sizes[v];
      ++part_conn_iter[p];
    } 
  }

  printf("Total number of components: %d\n", num_comps);
  for (int p = 0; p < num_parts; ++p) {
    printf("  Partition %2d has %4d component(s)\n", p, part_conns[p]);
    for (int i = 0; i < part_conns[p]; ++i) {
      printf("%4d  ", part_conn_sizes[p][i]);
    }
    printf("\n");
  }

  free(visited);
  free(queue);
  free(next_queue);
  free(part_conns);
  free(conn_sizes);

  return conn;
}


extern "C" int pulp_run(pulp_graph_t* g, pulp_part_control_t* ppc, 
          int* parts, int num_parts)
{
  srand(time(0));
  double vert_balance = ppc->vert_balance;
  double vert_balance_lower = 0.25;
  double edge_balance = ppc->edge_balance;
  double do_label_prop = ppc->do_lp_init;
  double do_nonrandom_init = ppc->do_bfs_init;
  double verbose = ppc->verbose_output;
  bool do_vert_balance = true;
  bool do_edge_balance = ppc->do_edge_balance;
  bool do_maxcut_balance = ppc->do_maxcut_balance;
  int balance_outer_iter = 1;
  int label_prop_iter = 3;
  int vert_outer_iter = 3;
  int vert_balance_iter = 5;
  int vert_refine_iter = 10;
  int edge_outer_iter = 3;
  int edge_balance_iter = 5;
  int edge_refine_iter = 10;
  seed = ppc->pulp_seed;

  double elt, elt2, elt3;
  elt = timer();

  if (do_label_prop && 
        g->vertex_weights == NULL && g->edge_weights == NULL)
  {
    if (verbose) printf("\tDoing label prop stage with %d parts\n", num_parts);
    elt2 = timer();
    label_prop(*g, num_parts, parts, label_prop_iter, vert_balance_lower);
    elt2 = timer() - elt2;
    if (verbose) printf("done: %9.6lf(s)\n", elt2);
  }
  else if (do_label_prop && 
           (g->vertex_weights != NULL || g->edge_weights != NULL))
  {
    if (verbose) printf("\tDoing (weighted) label prop stage with %d parts\n", num_parts);
    elt2 = timer();
    label_prop_weighted(*g, num_parts, parts, 
      label_prop_iter, vert_balance_lower);
    elt2 = timer() - elt2;
    if (verbose) printf("done: %9.6lf(s)\n", elt2);
  }
  else if (do_nonrandom_init)
  {
    if (verbose) printf("\tDoing bfs init stage with %d parts\n", num_parts);
    elt2 = timer();
    init_nonrandom(*g, num_parts, parts);
    elt2 = timer() - elt2;
    if (verbose) printf("done: %9.6lf(s)\n", elt2);
  }


  if (verbose) printf("\tBeginning vertex (and edge) refinement\n");
  for (int boi = 0; boi < balance_outer_iter; ++boi)
  {
    elt2 = timer();
    if (do_vert_balance && 
        g->vertex_weights == NULL && g->edge_weights == NULL)
    {
      if (verbose) printf("\t\tDoing vert balance and refinement stage\n");
      elt3 = timer(); 
      label_balance_verts(*g, num_parts, parts,
        vert_outer_iter, vert_balance_iter, vert_refine_iter,
        vert_balance);
      elt3 = timer() - elt3;
      if (verbose) printf("done: %9.6lf(s)\n", elt3);
    }
    else if (do_vert_balance)
    {
      if (verbose) printf("\t\tDoing (weighted) vert balance and refinement stage\n");
      elt3 = timer(); 
      label_balance_verts_weighted(*g, num_parts, parts,
        vert_outer_iter, vert_balance_iter, vert_refine_iter,
        vert_balance);
      elt3 = timer() - elt3;
      if (verbose) printf("done: %9.6lf(s)\n", elt3);
    }

    if (do_edge_balance && !do_maxcut_balance && 
        g->vertex_weights == NULL && g->edge_weights == NULL)
    {
      if (verbose) printf("\t\tDoing edge balance and refinement stage\n");
      elt3 = timer();
      label_balance_edges(*g, num_parts, parts,
        edge_outer_iter, edge_balance_iter, edge_refine_iter,
        vert_balance, edge_balance);
      elt3 = timer() - elt3;
      if (verbose) printf("done: %9.6lf(s)\n", elt3);
    }
    else if (do_edge_balance && do_maxcut_balance && 
        g->vertex_weights == NULL && g->edge_weights == NULL)
    {
      if (verbose) printf("\t\tDoing maxcut balance and refinement stage\n");
      elt3 = timer();
      label_balance_edges_maxcut(*g, num_parts, parts,
        edge_outer_iter, edge_balance_iter, edge_refine_iter,
        vert_balance, edge_balance);
      elt3 = timer() - elt3;
      if (verbose) printf("done: %9.6lfs\n", elt3);
    }
    else if (do_edge_balance && !do_maxcut_balance && 
             (g->vertex_weights != NULL || g->edge_weights != NULL))
    {
      if (verbose) printf("\t\tDoing (weighted) edge balance and refinement stage\n");
      elt3 = timer();
      label_balance_edges_weighted(*g, num_parts, parts,
        edge_outer_iter, edge_balance_iter, edge_refine_iter,
        vert_balance, edge_balance);
      elt3 = timer() - elt3;
      if (verbose) printf("done: %9.6lf(s)\n", elt3);
    }    
    else if (do_edge_balance && do_maxcut_balance && 
             (g->vertex_weights != NULL || g->edge_weights != NULL))
    {
      if (verbose) printf("\t\tDoing (weighted) maxcut balance and refinement stage\n");
      elt3 = timer();
      label_balance_edges_maxcut_weighted(*g, num_parts, parts,
        edge_outer_iter, edge_balance_iter, edge_refine_iter,
        vert_balance, edge_balance);
      elt3 = timer() - elt3;
      if (verbose) printf("done: %9.6lfs\n", elt3);
    }

    elt2 = timer() - elt2;
    if (verbose) printf("\tFinished outer loop iter %d: %9.6lf(s)\n", (boi+1), elt2);
  
    // Show connectivity details after each iteration
    connectivity_bfs(*g, num_parts, parts);
  }

  elt = timer() - elt;
  if (verbose) printf("Partitioning finished: %9.6lf(s)\n", elt);

  return 0;
}


double timer() 
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}

void evaluate_quality(pulp_graph_t& g, int num_parts, int* parts)
{
  for (int i = 0; i < g.n; ++i)
    if (parts[i] < 0)
    {
      printf("invalid part: %d %d\n", i, parts[i]);
      exit(0);
    }

  double comms_frac = 0.0;

  int num_verts = g.n;
  unsigned num_edges = g.m;
  unsigned num_comms = 0;
  bool** neighborhoods = new bool*[num_parts];
  bool** comms = new bool*[num_parts];
  long* part_sizes = new long[num_parts];
  int* num_comms_out = new int[num_parts];
  long* edge_cuts = new long[num_parts];
  int* boundary_verts = new int[num_parts];
  bool** part_to_part = new bool*[num_parts];
  unsigned* edges_per_part = new unsigned[num_parts];
  bool weighted = (g.vertex_weights_sum > 0);

  for (int i = 0; i < num_parts; ++i)
  {
    part_sizes[i] = 0;
    num_comms_out[i] = 0;
    edge_cuts[i] = 0;
    edges_per_part[i] = 0;
    boundary_verts[i] = 0;

    neighborhoods[i] = new bool[num_verts];
    comms[i] = new bool[num_verts];
    for (int j = 0; j < num_verts; ++j)
    {
      neighborhoods[i][j] = false;
      comms[i][j] = false;
    }

    part_to_part[i] = new bool[num_parts];
    for (int j = 0; j < num_parts; ++j)
      part_to_part[i][j] = false;
  }

  for (int v = 0; v < num_verts; ++v)
  { 
    if (weighted)
      part_sizes[parts[v]] += g.vertex_weights[v];
    else
      ++part_sizes[parts[v]];

    int part = parts[v];
    neighborhoods[part][v] = true;
    bool boundary = false;

    int out_degree = out_degree(g, v);
    int* outs = out_vertices(g, v);
    int* weights;
    if (weighted)
      weights = out_weights(g, v);
    for (int j = 0; j < out_degree; ++j)
    {
      int out = outs[j];
      neighborhoods[part][out] = true;

      int out_part = parts[out];
      if (out_part != part)
      {
        comms[part][out] = true;
        part_to_part[part][out_part] = true;

        if (weighted)
          edge_cuts[part] += weights[j];
        else
          ++edge_cuts[part];

        boundary = true;
      }
      ++edges_per_part[part];
    }
    if (boundary)
      ++boundary_verts[part];
  }

  for (int i = 0; i < num_parts; ++i)
  {
    for (int j = 0; j < num_verts; ++j)
    {
      if (comms[i][j])
      {
        ++num_comms_out[i];
        ++num_comms;
      }
    }
    for (int j = 0; j < num_parts; ++j)
      if (part_to_part[i][j])
        ++comms_frac;
  }

  long quality = 0;
  long edge_cut = 0;
  long max_vert_size = 0;
  unsigned max_edge_size = 0.0;
  int max_comm_size = 0;
  unsigned max_edge_cut = 0;
  int max_bound = 0;
  unsigned num_bound = 0;
  for (int i = 0; i < num_parts; ++i)
  {
#if VERBOSE
    printf("p: %d, v: %li, e: %u, com: %d, cut: %li, bound: %d\n", 
      i, part_sizes[i], 
      edges_per_part[i], num_comms_out[i], edge_cuts[i], boundary_verts[i]);
#endif
    quality += num_comms_out[i];
    edge_cut += edge_cuts[i];
    num_bound += boundary_verts[i];

    if (edge_cuts[i] > max_edge_cut)
      max_edge_cut = edge_cuts[i];
    if (part_sizes[i] > max_vert_size)
      max_vert_size = part_sizes[i];
    if (edges_per_part[i] > max_edge_size)
      max_edge_size = edges_per_part[i];
    if (num_comms_out[i] > max_comm_size)
      max_comm_size = num_comms_out[i];
    if (boundary_verts[i] > max_bound)
      max_bound = boundary_verts[i];
  }

  comms_frac = comms_frac / (double)(num_parts*(num_parts-1));

  long avg_size_vert;
  if (weighted) 
    avg_size_vert = g.vertex_weights_sum / (long)num_parts;
  else 
    avg_size_vert = num_verts / (unsigned)num_parts;
  unsigned avg_size_edge = num_edges / (unsigned)num_parts;
  unsigned avg_comm_size = num_comms / (unsigned)num_parts;
  unsigned avg_edge_cut = edge_cut / (unsigned)num_parts;
  unsigned avg_bound = num_bound / (unsigned)num_parts;
  double max_overweight_v = (double)max_vert_size/(double)avg_size_vert;
  double max_overweight_e = (double)max_edge_size/(double)avg_size_edge;
  double max_overweight_cv = (double)max_comm_size/(double)avg_comm_size;
  double max_overweight_ec = (double)max_edge_cut/(double)avg_edge_cut;
  double max_overweight_b = (double)max_bound/(double)avg_bound;
  edge_cut /= 2;
  long unsigned comVol = (long unsigned)quality;
  double comVolRatio = (double)comVol / (double)(num_edges/2);
  long unsigned edgeCut = (long unsigned)edge_cut;
  double edgeCutRatio = (double)edgeCut / (double)(num_edges/2);
  double boundVertRatio = (double)num_bound / (double)(num_verts);

  printf("Edge Cut: %lu\n", edgeCut);
  printf("Max Cut: %u\n", max_edge_cut);
  printf("Comm Vol: %lu\n", comVol);
  printf("Num boundary verts: %u\n", num_bound);
  printf("Comm ratio: %9.3lf\n", comVolRatio);
  printf("Edge ratio: %9.3lf\n", edgeCutRatio);
  printf("Boundary ratio: %9.3lf\n", boundVertRatio);
  printf("Vert overweight: %9.3lf\n", max_overweight_v);
  printf("Edge overweight: %9.3lf\n", max_overweight_e);
  printf("Boundary overweight: %9.3lf\n", max_overweight_b);
  printf("CommVol overweight: %9.3lf, max: %u\n", max_overweight_cv, max_comm_size);
  printf("EdgeCut overweight: %9.3lf, max: %u\n", max_overweight_ec, max_edge_cut);

  for (int i = 0; i < num_parts; ++i)
  {
    delete [] neighborhoods[i];
    delete [] comms[i];
    delete [] part_to_part[i];
  }
  delete [] neighborhoods;
  delete [] comms;
  delete [] part_to_part;
  delete [] part_sizes;
  delete [] num_comms_out;
  delete [] edge_cuts;
  delete [] boundary_verts;
  delete [] edges_per_part;
}


