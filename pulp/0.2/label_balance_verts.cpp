/*
//@HEADER
// *****************************************************************************
//
// PULP: Multi-Objective Multi-Constraint Partitioning Using Label Propagation
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
#include<fstream> // for component merging stats output to file
using namespace std;

void get_cut_overweight(pulp_graph_t& g, int num_parts, int* parts,
                        long* cut, double* overweight)
{
  int num_verts = g.n;
  long* part_sizes = new long[num_parts];
  long* edge_cuts = new long[num_parts];
  bool weighted = (g.vertex_weights_sum > 0);

  for (int i = 0; i < num_parts; ++i)
  {
    part_sizes[i] = 0;
    edge_cuts[i] = 0;
  }

  for (int v = 0; v < num_verts; ++v)
  { 
    if (weighted)
      part_sizes[parts[v]] += g.vertex_weights[v];
    else
      ++part_sizes[parts[v]];

    int part = parts[v];

    int out_degree = out_degree(g, v);
    int* outs = out_vertices(g, v);
    int* weights;
    if (weighted)
      weights = out_weights(g, v);
    for (int j = 0; j < out_degree; ++j)
    {
      int out = outs[j];
      int out_part = parts[out];
      if (out_part != part)
      {
        if (weighted)
          edge_cuts[part] += weights[j];
        else
          ++edge_cuts[part];
      }
    }
  }

  long edge_cut = 0;
  long max_vert_size = 0;
  for (int i = 0; i < num_parts; ++i)
  {
    edge_cut += edge_cuts[i];

    if (part_sizes[i] > max_vert_size)
      max_vert_size = part_sizes[i];
  }

  long avg_size_vert;
  if (weighted) 
    avg_size_vert = g.vertex_weights_sum / (long)num_parts;
  else 
    avg_size_vert = num_verts / (unsigned)num_parts;
  double max_overweight_v = (double)max_vert_size/(double)avg_size_vert;
  edge_cut /= 2;
  long unsigned edgeCut = (long unsigned)edge_cut;

  *cut = edgeCut;
  *overweight = max_overweight_v;

  delete [] part_sizes;
  delete [] edge_cuts;
}

/*
 * Find the true connectivity of vertex v
 * Which may not be equal to just conn[v] due to
 *   being in the middle of merging.
 */
int true_conn(int v, int* conn)
{
  while (v != conn[v]) {
    v = conn[v];
  }
  return v;
}

/*
 *  Find a neighboring component to component c_id,
 *      which will be merged together.
 *
 *  c_id is also the root vertex of this component
 *      (will be used as the root for BFS)
 *
 *  Return the number of neighboring components that will be merged.
 *  The IDs of these comps will be returns via int* merging_comp_ids
 */
int find_neighbor_part(pulp_graph_t& g, int num_parts, int* parts, int* part_sizes,
    int* conn_sizes, int c_id, int part, int* conn, int** merging_comp_ids)
{
  /*  We want a more intelligent way of choosing the neighboring component
   *    to merge, as opposed to just finding the closest one via BFS
   *  
   *  We'll do BFS on the entire component, keeping track of incident
   *    edges/vertices that do not belong to this component/part.
   *
   *  If we see a lot of neighboring components that are different, but
   *    belong to the same part, it can be beneficial to turn this
   *    component to belong to that part so that it connects all of them.
   *
   *  We should also consider with part within the entire graph is of the
   *    smallest current size, as merging with that part will best keep
   *    the vertex-overweightness low.
   *
   *  If we see a lot of incident edges for a certain part or component,
   *    it makes sense to merge with that part/component.
   *
   *  So, stats we want to collect before choosing who to merge with:
   *    Each neighboring component:
   *        size of that component
   *        how many times we are incident to that component
   *    Each part other than this part:
   *        how many vertices in the graph are in this part
   *
   */
  
  // neighboring component ids, per part
  int** neighbor_ids_per_part = (int**) malloc(num_parts * sizeof(int*));
  // Number of neighboring components per part
  int* num_neighbor_comps = (int*) malloc(num_parts * sizeof(int));
  int* max_num_neighbor_comps = (int*) malloc(num_parts * sizeof(int));

  // total sizes of neighboring components per part
  int* total_neighbor_part_sizes = (int*) malloc(num_parts * sizeof(int));
  
  // total number of neighboring incident edges (edge cut) per part
  int* edge_cut_per_part = (int*) malloc(num_parts * sizeof(int));

  // init
  for (int p = 0; p < num_parts; ++p) {
    max_num_neighbor_comps[p] = 2;
    num_neighbor_comps[p] = 0;
    neighbor_ids_per_part[p] = (int*) malloc(max_num_neighbor_comps[p] * sizeof(int));
    total_neighbor_part_sizes[p] = 0;
    edge_cut_per_part[p] = 0;
  }
  
  // Run BFS
  bool visited[g.n];
  for (int v = 0; v < g.n; ++v) {
    visited[v] = false;
  }
  visited[c_id] = true;
  
  int* queue = (int*) malloc(g.n * sizeof(int));
  queue[0] = c_id;
  int queue_size = 1;
  
  int* next_queue = (int*) malloc(g.n * sizeof(int));
  int next_size = 0;
  
  while (queue_size) {
    for (int i = 0; i < queue_size; ++i) {
      int vert = queue[i];

      for (int j = 0; j < out_degree(g, vert); ++j) {
        int adj = (out_vertices(g, vert))[j];
        if (!visited[adj]) {
          if (part != parts[adj]) {
            // adj belongs to a different component
            ++edge_cut_per_part[parts[adj]];
            
            int adj_c_id = true_conn(adj, conn);

            // See if we've already found this component
            bool found = false;
            for (int k = 0; k < num_neighbor_comps[parts[adj]]; ++k) {
              if (neighbor_ids_per_part[parts[adj]][k] == adj_c_id)
                found = true;
            }
            // We haven't found this component before, add stats about it
            if (!found) {
              // Resize neighboring ID list if needed
              if (max_num_neighbor_comps[parts[adj]] == num_neighbor_comps[parts[adj]]) {
                // Resize
                max_num_neighbor_comps[parts[adj]] *= 2;
                neighbor_ids_per_part[parts[adj]] = 
                  (int*) realloc(neighbor_ids_per_part[parts[adj]],
                      max_num_neighbor_comps[parts[adj]] * sizeof(int));
              }
              // Add this component to the ID list
              neighbor_ids_per_part[parts[adj]][num_neighbor_comps[parts[adj]]++] = adj_c_id;
              // Add this component's size to the total
              total_neighbor_part_sizes[parts[adj]] += conn_sizes[adj_c_id];
            } 
          } else {
            next_queue[next_size++] = adj;
          }
          visited[adj] = true;
        }
      }
    }

    int* temp = queue;
    queue = next_queue;
    next_queue = temp;
    queue_size = next_size;
    next_size = 0;
  }
 
  // Now we have a bunch of stats based on the parts neighboring this one
  // Decide which part to merge with and return the corresponding neighboring comp ids
  int chosen_part = -1;
  int num_merging_comps = -1;

  /*
   * One method:
   *    just pick part of the lowest size (still has to have a neighboring comp) 
   */  
  int min_size = g.n;
  for (int p = 0; p < num_parts; ++p) {
    if (min_size > part_sizes[p] && p != part) {
      min_size = part_sizes[p];
      chosen_part = p;
    }
  }

  /*
   * Second method:
   *    pick the part that has the highest number of neighboring components
   *
  int max_comps = 0;
  for (int p = 0; p < num_parts; ++p) {
    if (max_comps < num_neighbor_comps[p] && p != part) {
      max_comps = num_neighbor_comps[p];
      chosen_part = p;
    }
  }
  */

  /*
   * Third method:
   *    pick the part with the highest number of incident edges
   *
  int max_edges = 0;
  for (int p = 0; p < num_parts; ++p) {
    if (max_edges < edge_cut_per_part[p] && p != part) {
      max_edges = edge_cut_per_part[p];
      chosen_part = p;
    }
  }
  */

  // Return the component ids that will be merged with the chosen part
  if (chosen_part != -1) {
    *merging_comp_ids = neighbor_ids_per_part[chosen_part];
    num_merging_comps = num_neighbor_comps[chosen_part];
  }

  // Free memory
  for (int p = 0; p < num_parts; ++p) {
    if (p == chosen_part)
      continue;
    free(neighbor_ids_per_part[p]);
  }
  free(neighbor_ids_per_part);
  free(max_num_neighbor_comps);
  free(num_neighbor_comps);
  free(total_neighbor_part_sizes);
  free(edge_cut_per_part);
  free(queue);
  free(next_queue);
  return num_merging_comps;
}

void merge_small_components(pulp_graph_t& g, int num_parts, int* parts, int* part_sizes)
{
  bool debug = false;

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

  // Track size of largest component in each part
  int* max_conn_sizes = (int*) malloc(num_parts*sizeof(int));
  for (int p = 0; p < num_parts; ++p) {
    max_conn_sizes[p] = 0;
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

      // Update max component size for the part
      if (conn_size > max_conn_sizes[parts[v]])
        max_conn_sizes[parts[v]] = conn_size;
      // increment number of components for the part
      ++part_conns[parts[v]];

      conn_sizes[v] = conn_size;
      ++num_comps;
    }
  }

  // Per part, store component IDs and sizes
  int** part_conn_ids = (int**) malloc(num_parts*sizeof(int*));
  int** part_conn_sizes = (int**) malloc(num_parts*sizeof(int*));
  for (int p = 0; p < num_parts; ++p) {
    part_conn_ids[p] = (int*) malloc(part_conns[p]*sizeof(int));
    part_conn_sizes[p] = (int*) malloc(part_conns[p]*sizeof(int));
  }

  int* part_conn_iter = (int*) malloc(num_parts*sizeof(int));
  for (int p = 0; p < num_parts; ++p) {
    part_conn_iter[p] = 0;
  }

  for (int v = 0; v < g.n; ++v) {
    if (conn[v] == v) {
      int p = parts[v];
      part_conn_ids[p][part_conn_iter[p]] = v;
      part_conn_sizes[p][part_conn_iter[p]] = conn_sizes[v];
      ++part_conn_iter[p];
    } 
  }

  if (debug) {
    printf("BEFORE -- Total number of components: %d\n", num_comps);
    for (int p = 0; p < num_parts; ++p) {
      printf("  Partition %2d has %4d component(s) -- max size: %4d\n", 
             p, part_conns[p], max_conn_sizes[p]);
      for (int i = 0; i < part_conns[p]; ++i) {
        if (part_conn_sizes[p][i] > 0)
          printf("%4d  ", part_conn_sizes[p][i]);
      }
      printf("\n");
    }
  }

  long edge_cut_before;
  double vert_overweight_before;
  get_cut_overweight(g, num_parts, parts, &edge_cut_before, &vert_overweight_before);
  
  // Merge all small components until desired number of parts is reached
  int max_iters = 10;
  int iter = 0;
  while (num_comps > num_parts && iter < max_iters) {
    // Iterate over each component in each part
    for (int p = 0; p < num_parts; ++p) {
      for (int c = 0; c < part_conns[p]; ++c) {
        // Only merge components smaller than the largest one in their part
        //  and if they have not already been merged
        if (part_conn_sizes[p][c] < max_conn_sizes[p] &&
            part_conn_sizes[p][c] > 0) {
          int c_id = part_conn_ids[p][c];
          int part = parts[c_id];

          int* merging_comp_ids = NULL;
          int num_merging_comps = find_neighbor_part(g, num_parts, parts, part_sizes, conn_sizes, c_id, part, conn, &merging_comp_ids);

          // No neighbor for this component. Merging not possible
          if (num_merging_comps <= 0 || merging_comp_ids == NULL)
            continue;

          /* Merge c with these components, which all belong to the same part
           *    All these neighboring components will effectively merge with
           *    the main component. So, each of them will be relabeled to have
           *    connectivity c_id, and the size of c_id will increase by the
           *    sizes of all these neighbors' sizes.
           *    The number of components will decrease by num_merging_comps.
           */

          int new_part = parts[merging_comp_ids[0]];
          parts[c_id] = new_part;

          // BFS to change parts[]
          queue[0] = c_id;
          queue_size = 1;
          next_size = 0;

          while (queue_size) {
            for (int i = 0; i < queue_size; ++i) {
              int vert = queue[i];
              
              for (int j = 0; j < out_degree(g, vert); ++j) {
                int adj = (out_vertices(g, vert))[j];
                
                if (parts[adj] == part) {
                  parts[adj] = new_part;
                  next_queue[next_size++] = adj;
                }
              }
            }
            int* temp = queue;
            queue = next_queue;
            next_queue = temp;
            queue_size = next_size;
            next_size = 0;
          }

          // Total up size of the new merged component
          int num_orig_verts = part_conn_sizes[p][c];
          int total_merged_size = num_orig_verts;
          
          /* For each merging component,
           *    update part_conn_ids[] and part_conn_sizes[]
           *    and update true connectivity
           */
          for (int mc = 0; mc < num_merging_comps; ++mc) {
            // mc_id is the root of this component
            int mc_id = merging_comp_ids[mc];

            // Find this id in part_conn_ids[] to update part_conn_sizes[]
            for (int cc = 0; cc < part_conns[new_part]; ++cc) {
              if (part_conn_ids[new_part][cc] == mc_id) {
                total_merged_size += part_conn_sizes[new_part][cc];
                
                // If this is the last merging component, have its index
                //  be taken over by the overall merged component
                if (mc == num_merging_comps - 1) {
                  part_conn_sizes[new_part][cc] = total_merged_size;
                  part_conn_sizes[p][c] = 0;
                  part_conn_ids[new_part][cc] = c_id;
                  part_conn_ids[p][c] = mc_id;

                  // Update max comp size if necessary
                  if (total_merged_size > max_conn_sizes[new_part]) {
                    max_conn_sizes[new_part] = total_merged_size;
                  }
                
                } else {
                  part_conn_sizes[new_part][cc] = 0;
                }
                
                break;
              }
            }

            // Set connectivity id of mc_id to be c_id
            conn[true_conn(mc_id, conn)] = c_id;
          }

          /* 
           * Merging is done, update number of comps
           *  and part_sizes
           */
          num_comps -= num_merging_comps;
          part_sizes[part] -= num_orig_verts;
          part_sizes[new_part] += num_orig_verts;
        
          free(merging_comp_ids);
        } // end If 
      } // end c loop
    } // end p loop
    ++iter;
  } // end while

  if (debug) {
    printf("AFTER  -- Total number of components: %d\n", num_comps);
    for (int p = 0; p < num_parts; ++p) {
      printf("  Partition %2d has %4d component(s)\n", p, part_conns[p]);
      for (int i = 0; i < part_conns[p]; ++i) {
        if (part_conn_sizes[p][i] > 0)
          printf("%4d  ", part_conn_sizes[p][i]);
      }
      printf("\n");
    }
  }
  
  long edge_cut_after;
  double vert_overweight_after;
  get_cut_overweight(g, num_parts, parts, &edge_cut_after, &vert_overweight_after);
  printf("\nEDGE_CUT:          %8ld -> %8ld  | %12ld\n", 
         edge_cut_before, edge_cut_after, edge_cut_after-edge_cut_before);
  printf("VERT_OVERWEIGHT:   %8f -> %8f  | %12f\n", 
         vert_overweight_before, vert_overweight_after, 
         vert_overweight_after-vert_overweight_before);

  // Output edge cut and vert overweightness to file
  ofstream output;
  output.open("merge_stats.csv", ios::app);
  output << edge_cut_before << ", " << edge_cut_after << ", " <<
          edge_cut_after-edge_cut_before << ", " << vert_overweight_before <<
          ", " << vert_overweight_after << ", " << 
          vert_overweight_after-vert_overweight_before << "\n";
  output.close();

  free(visited);
  free(queue);
  free(next_queue);
  free(part_conns);
  free(conn_sizes);
}

/*
'########:::::'###::::'##::::::::::'##::::'##:'########:'########::'########:
 ##.... ##:::'## ##::: ##:::::::::: ##:::: ##: ##.....:: ##.... ##:... ##..::
 ##:::: ##::'##:. ##:: ##:::::::::: ##:::: ##: ##::::::: ##:::: ##:::: ##::::
 ########::'##:::. ##: ##:::::::::: ##:::: ##: ######::: ########::::: ##::::
 ##.... ##: #########: ##::::::::::. ##:: ##:: ##...:::: ##.. ##:::::: ##::::
 ##:::: ##: ##.... ##: ##:::::::::::. ## ##::: ##::::::: ##::. ##::::: ##::::
 ########:: ##:::: ##: ########::::::. ###:::: ########: ##:::. ##:::: ##::::
........:::..:::::..::........::::::::...:::::........::..:::::..:::::..:::::
*/
void label_balance_verts(pulp_graph_t& g, int num_parts, int* parts,
  int vert_outer_iter, int vert_balance_iter, int vert_refine_iter,
  double vert_balance)
{
  int num_verts = g.n;
  int* part_sizes = new int[num_parts];
 
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  double avg_size = num_verts / num_parts;
  int num_swapped_1 = 0;
  int num_swapped_2 = 0;
  double max_v;
  double running_max_v = (double)num_verts;

  int* queue = new int[num_verts*QUEUE_MULTIPLIER];
  int* queue_next = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int queue_size = num_verts;
  int next_size = 0;
  int t = 0;
  int num_tries = 0;

#pragma omp parallel
{
  int* part_sizes_thread = new int[num_parts];
  for (int i = 0; i < num_parts; ++i) 
    part_sizes_thread[i] = 0;

#pragma omp for schedule(static) nowait
  for (int i = 0; i < num_verts; ++i)
    ++part_sizes_thread[parts[i]];

  for (int i = 0; i < num_parts; ++i) 
#pragma omp atomic
    part_sizes[i] += part_sizes_thread[i];

  delete [] part_sizes_thread;


  double* part_counts = new double[num_parts];
  double* part_weights = new double[num_parts];

  int thread_queue[ THREAD_QUEUE_SIZE ];
  int thread_queue_size = 0;
  int thread_start;

  for (int p = 0; p < num_parts; ++p)
  {        
    part_weights[p] = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
    if (part_weights[p] < 0.0)
      part_weights[p] = 0.0;
  }

while(t < vert_outer_iter)
{

#pragma omp for schedule(static) nowait
  for (int i = 0; i < num_verts; ++i)
    queue[i] = i;
#pragma omp for schedule(static)
  for (int i = 0; i < num_verts; ++i)
    in_queue_next[i] = false;

#pragma omp single
{
  num_swapped_1 = 0;
  queue_size = num_verts;
  next_size = 0;
}

  int num_iter = 0;
  while (/*swapped &&*/ num_iter < vert_balance_iter)
  {
#pragma omp for schedule(guided) reduction(+:num_swapped_1) nowait
    for (int i = 0; i < queue_size; ++i)
    {
      int v = queue[i];
      in_queue[v] = false;
      int part = parts[v];
      for (int p = 0; p < num_parts; ++p)
        part_counts[p] = 0.0;
  
      unsigned out_degree = out_degree(g, v);
      int* outs = out_vertices(g, v);
      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out = outs[j];
        int part_out = parts[out];
        part_counts[part_out] += out_degree(g, out);
        //part_counts[part_out] += 1.0;//out_degree(g, out);
      }
      
      int max_part = part;
      double max_val = 0.0;
      for (int p = 0; p < num_parts; ++p)
      {
        part_counts[p] *= part_weights[p];
        
        if (part_counts[p] > max_val)
        {
          max_val = part_counts[p];
          max_part = p;
        }
      }

      if (max_part != part)
      {
        parts[v] = max_part;
        ++num_swapped_1;
    #pragma omp atomic
        --part_sizes[part];
    #pragma omp atomic
        ++part_sizes[max_part];
        
        part_weights[part] = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
        part_weights[max_part] = vert_balance * avg_size / (double)part_sizes[max_part]  - 1.0;   
        
        if (part_weights[part] < 0.0)
          part_weights[part] = 0.0;
        if (part_weights[max_part] < 0.0)
          part_weights[max_part] = 0.0;

        if (!in_queue_next[v])
        {
          in_queue_next[v] = true;
          thread_queue[thread_queue_size++] = v;

          if (thread_queue_size == THREAD_QUEUE_SIZE)
          {
#pragma omp atomic capture
            thread_start = next_size += thread_queue_size;
            
            thread_start -= thread_queue_size;
            for (int l = 0; l < thread_queue_size; ++l)
              queue_next[thread_start+l] = thread_queue[l];
            thread_queue_size = 0;
          }
        }
        for (unsigned j = 0; j < out_degree; ++j)
        {
          if (!in_queue_next[outs[j]])
          {
            in_queue_next[outs[j]] = true;
            thread_queue[thread_queue_size++] = outs[j];

            if (thread_queue_size == THREAD_QUEUE_SIZE)
            {
#pragma omp atomic capture
              thread_start = next_size += thread_queue_size;
              
              thread_start -= thread_queue_size;
              for (int l = 0; l < thread_queue_size; ++l)
                queue_next[thread_start+l] = thread_queue[l];
              thread_queue_size = 0;
            }
          }
        }
      }
    }

#pragma omp atomic capture
    thread_start = next_size += thread_queue_size;
    
    thread_start -= thread_queue_size;
    for (int l = 0; l < thread_queue_size; ++l)
      queue_next[thread_start+l] = thread_queue[l];
    thread_queue_size = 0;

#pragma omp barrier

    ++num_iter;
#pragma omp single
{
#if VERBOSE
    printf("%d\n", num_swapped_1);
#endif
    int* temp = queue;
    queue = queue_next;
    queue_next = temp;
    bool* temp_b = in_queue;
    in_queue = in_queue_next;
    in_queue_next = temp_b;
    queue_size = next_size;
    next_size = 0;

    num_swapped_1 = 0;

#if OUTPUT_STEP
  evaluate_quality_step(g, "VertBalance", parts, num_parts);
#endif
}
  } // end while

#pragma omp for schedule(static)
  for (int i = 0; i < num_verts; ++i)
    queue[i] = i;

#pragma omp single
{
  num_swapped_2 = 0;
  queue_size = num_verts;
  next_size = 0;
}

  num_iter = 0;
  while (/*swapped &&*/ num_iter < vert_refine_iter)
  {
#pragma omp for schedule(guided) reduction(+:num_swapped_2) nowait  
    for (int i = 0; i < queue_size; ++i)
    {
      int v = queue[i];
      in_queue[v] = false;
      for (int p = 0; p < num_parts; ++p)
        part_counts[p] = 0;

      int part = parts[v];
      unsigned out_degree = out_degree(g, v);
      int* outs = out_vertices(g, v);
      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out = outs[j];
        int part_out = parts[out];
        part_counts[part_out]++;
      }

      int max_part = -1;
      int max_count = -1;
      for (int p = 0; p < num_parts; ++p)
        if (part_counts[p] > max_count)
        {
          max_count = part_counts[p];
          max_part = p;
        }

      if (max_part != part)
      {
        double new_max_imb = (double)(part_sizes[max_part] + 1) / avg_size;
        if ( new_max_imb < vert_balance)
        {
          ++num_swapped_2;
          parts[v] = max_part;
      #pragma omp atomic
          ++part_sizes[max_part];
      #pragma omp atomic
          --part_sizes[part];

          if (!in_queue_next[v])
          {
            in_queue_next[v] = true;
            thread_queue[thread_queue_size++] = v;

            if (thread_queue_size == THREAD_QUEUE_SIZE)
            {
#pragma omp atomic capture
              thread_start = next_size += thread_queue_size;
              
              thread_start -= thread_queue_size;
              for (int l = 0; l < thread_queue_size; ++l)
                queue_next[thread_start+l] = thread_queue[l];
              thread_queue_size = 0;
            }
          }
          for (unsigned j = 0; j < out_degree; ++j) 
          {
            if (!in_queue_next[outs[j]])
            {
              in_queue_next[outs[j]] = true;
              thread_queue[thread_queue_size++] = outs[j];

              if (thread_queue_size == THREAD_QUEUE_SIZE)
              {
#pragma omp atomic capture
                thread_start = next_size += thread_queue_size;
                
                thread_start -= thread_queue_size;
                for (int l = 0; l < thread_queue_size; ++l)
                  queue_next[thread_start+l] = thread_queue[l];
                thread_queue_size = 0;
              }
            }
          }
        }
      }
    }   

#pragma omp atomic capture
    thread_start = next_size += thread_queue_size;
    
    thread_start -= thread_queue_size;
    for (int l = 0; l < thread_queue_size; ++l)
      queue_next[thread_start+l] = thread_queue[l];
    thread_queue_size = 0;

#pragma omp barrier

    ++num_iter;
#pragma omp single
{
#if VERBOSE
    printf("%d\n", num_swapped_2);
#endif
    int* temp = queue;
    queue = queue_next;
    queue_next = temp;
    bool* temp_b = in_queue;
    in_queue = in_queue_next;
    in_queue_next = temp_b;
    queue_size = next_size;
    next_size = 0;

    num_swapped_2 = 0;

    max_v = 0.0;
    for (int p = 0; p < num_parts; ++p)
    {
      if ((double)part_sizes[p] / avg_size > max_v)
        max_v = (double)part_sizes[p] / avg_size;
    }
#if OUTPUT_STEP
  evaluate_quality_step(g, "VertRefine", parts, num_parts);
#endif
}
  } // end while

#pragma omp single
{
  if (max_v > vert_balance*1.01 && t == vert_outer_iter-1 && num_tries < 3)
  {
    --t;
    if (max_v < running_max_v)
    {
      running_max_v = max_v;
      printf("Vertex balance missed, attempting further iterations: (%2.3lf)\n", max_v);
    }
    else
      ++num_tries;
  }
  else
    ++t;

  // Check connectivity of parts and merge small components
  merge_small_components(g, num_parts, parts, part_sizes);
}
} // end for

  delete [] part_counts;
  delete [] part_weights;

} // end par


  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}





/*
'########:::::'###::::'##::::::::::'##::::'##:'########:'########::'########:
 ##.... ##:::'## ##::: ##:::::::::: ##:::: ##: ##.....:: ##.... ##:... ##..::
 ##:::: ##::'##:. ##:: ##:::::::::: ##:::: ##: ##::::::: ##:::: ##:::: ##::::
 ########::'##:::. ##: ##:::::::::: ##:::: ##: ######::: ########::::: ##::::
 ##.... ##: #########: ##::::::::::. ##:: ##:: ##...:::: ##.. ##:::::: ##::::
 ##:::: ##: ##.... ##: ##:::::::::::. ## ##::: ##::::::: ##::. ##::::: ##::::
 ########:: ##:::: ##: ########::::::. ###:::: ########: ##:::. ##:::: ##::::
........:::..:::::..::........::::::::...:::::........::..:::::..:::::..:::::
*/
void label_balance_verts_weighted(
  pulp_graph_t& g, int num_parts, int* parts,
  int vert_outer_iter, int vert_balance_iter, int vert_refine_iter,
  double vert_balance)
{
  int num_verts = g.n;
  long* part_sizes = new long[num_parts];

  bool has_vwgts = (g.vertex_weights != NULL);
  bool has_ewgts = (g.edge_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;
 
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  double avg_size = (double)g.vertex_weights_sum / (double)num_parts;
  int num_swapped_1 = 0;
  int num_swapped_2 = 0;
  double max_v;
  double running_max_v = (double)num_verts;

  int* queue = new int[num_verts*QUEUE_MULTIPLIER];
  int* queue_next = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int queue_size = num_verts;
  int next_size = 0;
  int t = 0;
  int num_tries = 0;

#pragma omp parallel
{
  long* part_sizes_thread = new long[num_parts];
  for (int i = 0; i < num_parts; ++i) 
    part_sizes_thread[i] = 0;

#pragma omp for schedule(static) nowait
  for (int i = 0; i < num_verts; ++i)
    if (has_vwgts)
      part_sizes_thread[parts[i]] += g.vertex_weights[i];
    else
      ++part_sizes_thread[parts[i]];

  for (int i = 0; i < num_parts; ++i) 
#pragma omp atomic
    part_sizes[i] += part_sizes_thread[i];

  delete [] part_sizes_thread;
#pragma omp barrier

  double* part_counts = new double[num_parts];
  double* part_weights = new double[num_parts];

  int thread_queue[ THREAD_QUEUE_SIZE ];
  int thread_queue_size = 0;
  int thread_start;

  for (int p = 0; p < num_parts; ++p)
  {        
    part_weights[p] = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
    if (part_weights[p] < 0.0)
      part_weights[p] = 0.0;
  }

while(t < vert_outer_iter)
{

#pragma omp for schedule(static) nowait
  for (int i = 0; i < num_verts; ++i)
    queue[i] = i;
#pragma omp for schedule(static)
  for (int i = 0; i < num_verts; ++i)
    in_queue_next[i] = false;

#pragma omp single
{
  num_swapped_1 = 0;
  queue_size = num_verts;
  next_size = 0;
}  

  int num_iter = 0;
  while (/*swapped &&*/ num_iter < vert_balance_iter)
  {
#pragma omp for schedule(guided) reduction(+:num_swapped_1) nowait
    for (int i = 0; i < queue_size; ++i)
    {
      int v = queue[i];
      in_queue[v] = false;
      int part = parts[v];
      int v_weight = 1;
      if (has_vwgts) v_weight = g.vertex_weights[v];

      for (int p = 0; p < num_parts; ++p)
        part_counts[p] = 0.0;

      unsigned out_degree = out_degree(g, v);
      int* outs = out_vertices(g, v);
      int* weights = out_weights(g, v);
      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out = outs[j];
        int part_out = parts[out];
        double weight_out = 1.0;
        if (has_ewgts) weight_out = (double)weights[j];
        part_counts[part_out] += (double)out_degree(g, out)*weight_out;
      }
      
      int max_part = part;
      double max_val = 0.0;
      for (int p = 0; p < num_parts; ++p)
      {
        part_counts[p] *= part_weights[p];
        
        if (part_counts[p] > max_val)
        {
          max_val = part_counts[p];
          max_part = p;
        }
      }

      if (max_part != part)
      {
        parts[v] = max_part;
        ++num_swapped_1;
    #pragma omp atomic
        part_sizes[max_part] += v_weight;
    #pragma omp atomic
        part_sizes[part] -= v_weight;
        
        part_weights[part] = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
        part_weights[max_part] = vert_balance * avg_size / (double)part_sizes[max_part]  - 1.0;   
        
        if (part_weights[part] < 0.0)
          part_weights[part] = 0.0;
        if (part_weights[max_part] < 0.0)
          part_weights[max_part] = 0.0;

        if (!in_queue_next[v])
        {
          in_queue_next[v] = true;
          thread_queue[thread_queue_size++] = v;

          if (thread_queue_size == THREAD_QUEUE_SIZE)
          {
#pragma omp atomic capture
            thread_start = next_size += thread_queue_size;
            
            thread_start -= thread_queue_size;
            for (int l = 0; l < thread_queue_size; ++l)
              queue_next[thread_start+l] = thread_queue[l];
            thread_queue_size = 0;
          }
        }
        for (unsigned j = 0; j < out_degree; ++j)
        {
          if (!in_queue_next[outs[j]])
          {
            in_queue_next[outs[j]] = true;
            thread_queue[thread_queue_size++] = outs[j];

            if (thread_queue_size == THREAD_QUEUE_SIZE)
            {
#pragma omp atomic capture
              thread_start = next_size += thread_queue_size;
              
              thread_start -= thread_queue_size;
              for (int l = 0; l < thread_queue_size; ++l)
                queue_next[thread_start+l] = thread_queue[l];
              thread_queue_size = 0;
            }
          }
        }
      }
    }

#pragma omp atomic capture
    thread_start = next_size += thread_queue_size;
    
    thread_start -= thread_queue_size;
    for (int l = 0; l < thread_queue_size; ++l)
      queue_next[thread_start+l] = thread_queue[l];
    thread_queue_size = 0;

#pragma omp barrier

    ++num_iter;
#pragma omp single
{
#if VERBOSE
    printf("%d\n", num_swapped_1);
#endif
    int* temp = queue;
    queue = queue_next;
    queue_next = temp;
    bool* temp_b = in_queue;
    in_queue = in_queue_next;
    in_queue_next = temp_b;
    queue_size = next_size;
    next_size = 0;

    num_swapped_1 = 0;

#if OUTPUT_STEP
  evaluate_quality_step(g, "VertBalance", parts, num_parts);
#endif
}
  } // end while

#pragma omp for schedule(static)
  for (int i = 0; i < num_verts; ++i)
    queue[i] = i;

#pragma omp single
{
  num_swapped_2 = 0;
  queue_size = num_verts;
  next_size = 0;
}

  num_iter = 0;
  while (/*swapped &&*/ num_iter < vert_refine_iter)
  {
#pragma omp for schedule(guided) reduction(+:num_swapped_2) nowait  
    for (int i = 0; i < queue_size; ++i)
    {
      int v = queue[i];
      in_queue[v] = false;      
      int part = parts[v];
      int v_weight = 1;
      if (has_vwgts) v_weight = g.vertex_weights[v];

      for (int p = 0; p < num_parts; ++p)
        part_counts[p] = 0;

      unsigned out_degree = out_degree(g, v);
      int* outs = out_vertices(g, v);
      int* weights = out_weights(g, v);
      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out = outs[j];
        int part_out = parts[out];
        int out_weight = 1;
        if (has_ewgts) out_weight = weights[j];
        part_counts[part_out] += out_weight;
      }

      int max_part = -1;
      int max_count = -1;
      for (int p = 0; p < num_parts; ++p)
        if (part_counts[p] > max_count)
        {
          max_count = part_counts[p];
          max_part = p;
        }

      if (max_part != part)
      {
        double new_max_imb = (double)(part_sizes[max_part] + v_weight) / avg_size;
        if (new_max_imb < vert_balance)
        {
          ++num_swapped_2;
          parts[v] = max_part;
      #pragma omp atomic
          part_sizes[max_part] += v_weight;
      #pragma omp atomic
          part_sizes[part] -= v_weight;

          if (!in_queue_next[v])
          {
            in_queue_next[v] = true;
            thread_queue[thread_queue_size++] = v;

            if (thread_queue_size == THREAD_QUEUE_SIZE)
            {
#pragma omp atomic capture
              thread_start = next_size += thread_queue_size;
              
              thread_start -= thread_queue_size;
              for (int l = 0; l < thread_queue_size; ++l)
                queue_next[thread_start+l] = thread_queue[l];
              thread_queue_size = 0;
            }
          }
          for (unsigned j = 0; j < out_degree; ++j) 
          {
            if (!in_queue_next[outs[j]])
            {
              in_queue_next[outs[j]] = true;
              thread_queue[thread_queue_size++] = outs[j];

              if (thread_queue_size == THREAD_QUEUE_SIZE)
              {
#pragma omp atomic capture
                thread_start = next_size += thread_queue_size;
                
                thread_start -= thread_queue_size;
                for (int l = 0; l < thread_queue_size; ++l)
                  queue_next[thread_start+l] = thread_queue[l];
                thread_queue_size = 0;
              }
            }
          }
        }
      }
    }   

#pragma omp atomic capture
    thread_start = next_size += thread_queue_size;
    
    thread_start -= thread_queue_size;
    for (int l = 0; l < thread_queue_size; ++l)
      queue_next[thread_start+l] = thread_queue[l];
    thread_queue_size = 0;

#pragma omp barrier

    ++num_iter;
#pragma omp single
{
#if VERBOSE
    printf("%d\n", num_swapped_2);
#endif
    int* temp = queue;
    queue = queue_next;
    queue_next = temp;
    bool* temp_b = in_queue;
    in_queue = in_queue_next;
    in_queue_next = temp_b;
    queue_size = next_size;
    next_size = 0;

    num_swapped_2 = 0;

    max_v = 0.0;
    for (int p = 0; p < num_parts; ++p)
    {
      if ((double)part_sizes[p] / avg_size > max_v)
        max_v = (double)part_sizes[p] / avg_size;
    }
#if OUTPUT_STEP
  evaluate_quality_step(g, "VertRefine", parts, num_parts);
#endif
}
  } // end while

#pragma omp single
{
  if (max_v > vert_balance*1.01 && t == vert_outer_iter-1 && num_tries < 3)
  {
    --t;
    if (max_v < running_max_v)
    {
      running_max_v = max_v;
      printf("Vertex balance missed, attempting further iterations: (%2.3lf)\n", max_v);
    }
    else
      ++num_tries;
  }
  else
    ++t;
}
} // end for

  delete [] part_counts;
  delete [] part_weights;

} // end par


  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}

