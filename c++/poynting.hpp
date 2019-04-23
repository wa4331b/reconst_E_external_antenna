#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void poyntingMatrix(int N_RWG,
                    double RWGNumber_trianglesCoord_tmp[],
                    int T,
                    int triangle_NRWG[],
                    int triangle_RWGNumber_tmp[],
                    int triangle_signInRWG_tmp[],
                    int triangle_surfaces[],
                    int N_target_surface,
                    int target_surface_tmp[],
                    complex<double> poyntingMatrix_tmp[]);
