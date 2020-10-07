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
    pointC p;
    bool active;
    
    cp(pointC P)
    {
        p = P;
        active = true;
    }
};

class simplex {
private:
    pointC pnts[4];
    bool active;
    size_t upperleft;
    size_t lowerleft;
    size_t upperright;
    size_t lowerright;
    
public:
    simplex(size_t upperleft, size_t lowerleft, size_t upperright, size_t lowerright);
    simplex(pointC p0, pointC p1, pointC p2, pointC p3);
    void points(std::vector<cp> &points);
    void subdivide(std::vector<simplex> &simplices, std::vector<cp> &points, const double delta);
    std::complex<double> integrate(const pointD mu, const double nu, const int N) const;
    
    bool isActive() const;
    pointC p0() const;
    pointC p1() const;
    pointC p2() const;
    pointC p3() const;
    size_t UL() const;
    size_t LL() const;
    size_t UR() const;
    size_t LR() const;
};

void externalVariables(std::vector<pointD> &mus, const double muMin, const double muMax, int NMu);
void initialize(std::vector<simplex> &simplices, std::vector<cp> &points, const double xMin, const double xMax, const double delta);
void subdevide(std::vector<simplex> &simplices, std::vector<cp> &points, const double delta);
void clean(std::vector<simplex> &simplices);
void flow(std::vector<simplex> &simplices, std::vector<cp> &points, const double step, const pointD mu, const double nu, const double thres, const double delta, const int Niterations);
void flow(std::vector<cp> &points, const double step, const pointD mu, const double nu, const double thres);
void importPoints(std::vector<simplex> &simplices, std::vector<cp> &points);

std::complex<double> func(const pointC p, const pointD mu, const double nu);
std::complex<double> phi(const pointC x);
double h(const std::complex<double> x, const double mu, const double nu);
double H(const std::complex<double> x, const double mu, const double nu);
std::complex<double> gradient(const std::complex<double> x, const double mu, const double nu, const double epsilon);

std::complex<double> integrate(const std::vector<simplex> &simplices, const pointD mu, const double nu, const int N);
std::complex<double> integrate(std::vector<simplex> &simplices, std::vector<cp> &points, const pointD mu, const double nu, const int N);
pointC map(const pointD p,  const pointC pnts[4]);
std::complex<double> jacob(const pointD p, const pointC pnts[4]);
std::complex<double> integrand(const double u, const double v, const pointC pnts[4], const pointD &mu, const double nu);
std::complex<double> romberg(const pointC pnts[4], const pointD &mu, const double nu, const int N);

int nearest(const std::vector<pointD> &list, const pointD &p);
double min(const double a, const double b);
double min(const double a, const double b, const double c);
double max(const double a, const double b);

void writeB(std::vector<std::complex<double> > &psi, std::string fileName);
void writeB(std::complex<double> number, std::ofstream &file);

#include "Exponent.h"
#include "Simplex.h"
#include "Integrate.h"
#include "Utility.h"
