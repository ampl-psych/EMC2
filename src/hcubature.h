#include <iostream>
#include <iomanip>
#include <vector>

using std::setw;
using std::cout;
using std::vector;



static const double xd7[8] = { -9.9145537112081263920685469752598e-01,
                               -9.4910791234275852452618968404809e-01,
                               -8.6486442335976907278971278864098e-01,
                               -7.415311855993944398638647732811e-01,
                               -5.8608723546769113029414483825842e-01,
                               -4.0584515137739716690660641207707e-01,
                               -2.0778495500789846760068940377309e-01,
                               0.0 };

static const double wd7[8] = { 2.2935322010529224963732008059913e-02,
                               6.3092092629978553290700663189093e-02,
                               1.0479001032225018383987632254189e-01,
                               1.4065325971552591874518959051021e-01,
                               1.6900472663926790282658342659795e-01,
                               1.9035057806478540991325640242055e-01,
                               2.0443294007529889241416199923466e-01,
                               2.0948214108472782801299917489173e-01 };

static const double gwd7[4] = { 1.2948496616886969327061143267787e-01,
                                2.797053914892766679014677714229e-01,
                                3.8183005050511894495036977548818e-01,
                                4.1795918367346938775510204081658e-01 };

struct one_d {
  double	result;
  double	err;
  int kdivide = 0;
};

struct GenzMalik {
  vector<vector<double>> p[4];
  double w[5];
  double wd[4];
};

void combination(int* c, int n, int p, int x);

int choose(int n, int k);

void gauss_kronrod(double a, double b, one_d& out, void* pars, int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval));

void combos(int k, double lambda, int n, vector<vector<double>>& p);

void increment(vector<bool>& index, int k, double lambda, int n, int* c, vector<double>& temp);

void signcombos(int k, double lambda, int n, vector<vector<double>>& p);

void make_GenzMalik(int n, GenzMalik& g);

void integrate_GenzMalik(GenzMalik g, int n, const double* a, const double* b, one_d& out, void* pars, double integrand(double* x, void* pars));

int hcubature(int integrand(unsigned dim, const double* x, void* p, unsigned fdim, double* retval), void* pars, unsigned n, const double* a, const double* b,
              size_t maxEval, double reqAbsError, double reqRelError, double* val, double* err);
