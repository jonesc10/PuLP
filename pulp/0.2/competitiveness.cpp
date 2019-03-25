// Features for Numbers of Democrats and Republicans in QGIS:
//   REG10G_D
//   REG10G_R

#include <fstream>
#include <string>
#include <sstream>

// File where each line has number of dems registered and reps registered,
//   separated by a space.
#define filename "demrep.tsv"

// Note that a lower value of competitiveness is desired:
// A value of zero indicates dem. and rep. registration is the same in all districts
//
// The competitiveness value is normalized to [0,1], so you could just reverse
//   it by adding one to the negative value
float competitiveness(pulp_graph_t& g, int num_parts, int* parts)
{
  float competitiveness = -1.0;
 
  ifstream infile;
  string line;
  string val;

  int n = g.n;

  int dem[n] = {0};
  int rep[n] = {0};

  infile.open(filename);

  getline(infile, line);

  int v = 0;
  while (getline(infile, line))
  {
    stringstream ss(line);
    getline(ss, val, ' ');
    dem[v] = atoi(val.c_str());
    getline(ss, val, ' ');
    rep[v] = atoi(val.c_str());
    ++v;
  }

  int total_dem[num_parts] = {0};
  int total_rep[num_parts] = {0};

  for (int v = 0; v < n; ++v)
  {
    total_dem[parts[v]] += dem[v];
    total_rep[parts[v]] += rep[v];
  }

  // Formula for competitiveness taken from PEAR paper
  float T_e = 0.0;
  float T_p = 0.0;
  float B_R = 0.0;
  for (int p = 0; p < num_parts; ++p)
  {
    if (total_rep[p] > total_dem[p]) B_R += 1.0;
    float deviation = (float)(total_rep[p])/(float)(total_dem[p]+total_rep[p]) - 0.5;
    if (deviation < 0) deviation = -deviation;
    T_p += deviation;
  }
  T_p = T_p / num_parts;
  T_e = B_R / num_parts - 0.5;
  if (T_e < 0) T_e = -T_e;

  competitiveness = T_p * (1.0 + T_e) * 4.0 / 3.0;

  return competitiveness;
}
