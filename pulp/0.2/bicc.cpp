/*Adapted from java code:*/
/*************************************************************************
 *  Compilation:  javac Biconnected.java
 *  Dependencies: Graph.java 
 *
 *  Identify articulation points and print them out.
 *  This can be used to decompose a graph into biconnected components.
 *  Runs in O(E + V) time.
 *
 *  http://www.cs.brown.edu/courses/cs016/book/slides/Connectivity2x2.pdf
 *
 ****/

class biconnected
{
public:

  biconnected(pulp_graph_t& g, int* parts)
  {
    int num_verts = g.n;
    g_parts = parts;

    cnt = 0;
    articulation_count = 0;

    low = new int[num_verts];
    pre = new int[num_verts];
    articulation = new bool[num_verts];
    for (int v = 0; v < num_verts; ++v) low[v] = -1;
    for (int v = 0; v < num_verts; ++v) pre[v] = -1;
    for (int v = 0; v < num_verts; ++v) articulation[v] = false;

    for (int v = 0; v < num_verts; ++v)
      if (pre[v] == -1)
        dfs(g, v, v);
  }

  ~biconnected()
  {
    delete [] low;
    delete [] pre;
    delete [] articulation;
  }

  int get_articulation_count(pulp_graph_t& g)
  {
    articulation_count = 0;
    int num_verts = g.n;
    for (int i = 0; i < num_verts; ++i)
      if (articulation[i])
        {
          ++articulation_count;
        }

    return articulation_count;
  }

  bool is_articulation(int v)
  {
    return articulation[v];
  }

private:
  void dfs(pulp_graph_t& g, int u, int v) 
  {
    int children = 0;
    pre[v] = cnt++;
    low[v] = pre[v];

    int* outs = out_vertices(g, v);
    unsigned out_degree = out_degree(g, v);
    int out;    
    for (unsigned i = 0; i < out_degree; ++i)
    {
      out = outs[i];

      // Only consider neighbors in the same part
      if (g_parts[v] != g_parts[out])
        continue;

      if (pre[out] == -1) 
      {
        children++;
        dfs(g, v, out);

        // update low number
        low[v] = low[v] < low[out] ? low[v] : low[out];

        // non-root of DFS is an articulation point if low[out] >= pre[v]
        if (low[out] >= pre[v] && u != v) 
        {
          articulation[v] = true;
          articulation_count++;
        }
      }

      // update low number - ignore reverse of edge leading to v
      else if (out != u)
        low[v] = low[v] < pre[out] ? low[v] : pre[out];
    }

    // root of DFS is an articulation point if it has more than 1 child
    if (u == v && children > 1)
    {
      articulation[v] = true;
      articulation_count++;
    }
  }

  int* g_parts;
  int* low;
  int* pre;
  int cnt;
  int articulation_count;
  bool* articulation;
  bool* visited;
};

