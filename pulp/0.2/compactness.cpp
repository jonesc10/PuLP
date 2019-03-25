// Features for Area and Perimeter in QGIS:
//  Shape_Area
//  Shape_Leng

#include <fstream>
#include <string>
#include <sstream>

// File where each line is a perimeter and an area, separated by a space
#define perimeter_area_file "perimeter_area.tsv"

void calculate_outside_lengths(pulp_graph_t& g,
                               long* perimeter,
                               long* area,
                               long* outside_perimeter)
{
  int n = g.n;

  for (int v = 0; v < n; ++v)
  {
    long total_edge_weights = 0;
    int d = out_degree(g, v);
    int* weights = out_weights(g, v);
    
    for (int j = 0; j < d; ++j)
    {
      total_edge_weights += (long) weights[j];
    }

    outside_perimeter[v] = perimeter[v] - total_edge_weights;
  }
}

float compactness(pulp_graph_t& g, int num_parts, int* parts)
{
  float avg_compactness = 0.0;
  
  ifstream infile;
  string line;
  string val;

  int n = g.n;

  long perimeter[n] = {0};
  long area[n] = {0};

  infile.open(perimeter_area_file);
  
  getline(infile, line);

  int v = 0;
  while (getline(infile, line))
  {
    stringstream ss(line);
    getline(ss, val, ' ');
    perimeter[v] = atol(val.c_str());
    getline(ss, val, ' ');
    area[v] = atol(val.c_str());
    ++v;
  }
  infile.close();

  long outside_perimeter[n] = {0};

  calculate_outside_lengths(g, perimeter, area, outside_perimeter);

  // Run BFS to get total perimeter and area for each part
  long total_perimeter[num_parts] = {0};
  long total_area[num_parts] = {0};

  int* queue = (int*) malloc(n*sizeof(int));
  int* next_queue = (int*) malloc(n*sizeof(int));
  int queue_size = 0;
  int next_size = 0;
  bool* visited = (bool*) malloc(n*sizeof(int));
  for (int v = 0; v < n; ++v)
  {
    visited[v] = false;
  }
  
  for (int v = 0; v < n; ++v)
  {
    if (!visited[v])
    {
      visited[v] = true;
      
      queue[0] = v;
      queue_size = 1;
      next_size = 1;
      
      while (queue_size)
      {
        for (int i = 0; i < queue_size; ++i)
        {
          int vert = queue[i];
          total_area[parts[vert]] += area[vert];
          total_perimeter[parts[vert]] += outside_perimeter[vert];

          for (int j = 0; j < out_degree(g, vert); ++j)
          {
            int adj = (out_vertices(g, vert))[j];
            if (parts[adj] != parts[v])
            {
              total_perimeter[parts[v]] += (long) (out_weights(g, vert))[j];
            }
            if (!visited[adj] && (parts[adj] == parts[v]))
            {
              visited[adj] = true;
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
    }
  }
  free(queue);
  free(next_queue);
  free(visited);

  // Now we have area and perimeter for each part
  float compactness_value[num_parts] = {0};
  for (int p = 0; p < num_parts; ++p)
  {
    compactness_value[p] = ((float) total_area[p])/((float) (total_perimeter[p]*total_perimeter[p]));
    avg_compactness += compactness_value[p] / (float)num_parts;
    printf("Part %2d Compactness = %12.4f\n", p, compactness_value[p]);
  }

  return avg_compactness;
}




