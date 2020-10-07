#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <cmath>
#include <cfloat>

#if defined(_OPENMP)
#include <omp.h>
extern const bool parallelism_enabled = true;
#else
extern const bool parallelism_enabled = false;
#endif

#include "Point.h"

const double epsilon = 0.00000001;
const double pi = 3.141592654;
const std::complex<double> I (0,1.0);

struct cp {
    std::complex<double> p;
    bool active;
    
    cp(std::complex<double> P)
    {
        p = P;
        active = true;
    }
};

class simplex {
private:
    std::complex<double> pnts[2];
    bool active;
    size_t left;
    size_t right;
    
public:
    simplex(size_t left, size_t right);
    simplex(std::complex<double> p0, std::complex<double> p1);
    void points(std::vector<cp> &points);
    void subdivide(std::vector<simplex> &simplices, std::vector<cp> &points, const double delta);
    std::complex<double> integrate(const double mu, const double nu, const int N) const;
    
    bool isActive();
    std::complex<double> p0();
    std::complex<double> p1();
    size_t L();
    size_t R();
};

void externalVariables(std::vector<double> &mus, const double muMin, const double muMax, int NMu);
void initialize(std::vector<simplex> &simplices, std::vector<cp> &points, const double xMin, const double xMax, const double delta);
void subdevide(std::vector<simplex> &simplices, std::vector<cp> &points, const double delta);
void clean(std::vector<simplex> &simplices);
void flow(std::vector<simplex> &simplices, std::vector<cp> &points, const double step, const double mu, const double nu, const double thres, const double delta, const int Niterations);
void flow(std::vector<cp> &points, const double step, const double mu, const double nu, const double thres);
void importPoints(std::vector<simplex> &simplices, std::vector<cp> &points);

std::complex<double> func(const std::complex<double> x, const double mu, const double nu);
std::complex<double> phi(const std::complex<double> x);
double h(const std::complex<double> x, const double mu, const double nu);
double H(const std::complex<double> x, const double mu, const double nu);
std::complex<double> gradient(const std::complex<double> x, const double mu, const double nu, const double epsilon);

std::complex<double> integrate(const std::vector<simplex> &simplices, const double mu, const double nu, const int N);
std::complex<double> integrate(std::vector<simplex> &simplices, std::vector<cp> &points, const double mu, const double nu, const int N);
std::complex<double> romberg(std::complex<double> a, std::complex<double> b, const double mu, const double nu, const int N);

int nearest(const std::vector<double> &list, const double &p);
void writeB(std::vector<std::complex<double> > &psi, std::string fileName);
void writeB(std::complex<double> number, std::ofstream &file);

#include "Exponent.h"
#include "Simplex.h"
#include "Integrate.h"
#include "Utility.h"
