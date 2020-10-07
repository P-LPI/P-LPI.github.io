#include "Header.h"

std::complex<double> func(const std::complex<double> x, const double mu, const double nu)
{
    return I * nu * (pow(x - mu, 2.) / 2. + phi(x));
}

std::complex<double> phi(const std::complex<double> x)
{
    return 2. / (1. + pow(x, 2.));
}

int main(int argc, const char * argv[])
{
    const double nu = 10;

    // Parameters
    const double xMin = -20.;
    const double xMax = +20.;
    const double delta = 0.25;
    
    const double step = 0.05;
    const double thres = -20.;
    
    const int Niterations = 50;
    
    const double muMin = -4.;
    const double muMax = +4.;
    const int NMu = pow(2, 5);
    
    const int N = 3;
    const int NM = pow(2, 9);
    
    {
        // Initial integration domain
        std::vector<simplex> simplicesInit;
        std::vector<cp> pointsInit;
        initialize(simplicesInit, pointsInit, xMin, xMax, delta);
        
        // Compute the Picard-Lefschetz thimbles
        std::vector<double> mus(NMu);
        externalVariables(mus, muMin, muMax, NMu);
        
        std::cout << "Picard-Lefschetz thimble:";
        std::vector<std::vector<simplex>> PL(mus.size());
#pragma omp parallel for
        for(int index = 0; index < mus.size(); index++)
        {
            if((index % (NMu / 16)) == 0){std::cout << "."; std::cout.flush();}
            std::vector<simplex> simplices = simplicesInit;
            std::vector<cp> points = pointsInit;
            flow(simplices, points, step, mus[index], nu, thres, delta, Niterations);
            clean(simplices);
            importPoints(simplices, points);
            PL[index] = simplices;
        }
        std::cout << std::endl;

        // Evaluate the integral
        std::cout << "Integrate along thimbles:";
        std::vector<std::complex<double>> result(NM);
#pragma omp parallel for
        for(int iMu = 0; iMu < NM; iMu++)
        {
            if((iMu % (NM / 16)) == 0){std::cout << "."; std::cout.flush();}
            const double mu = (muMax - muMin) / (NM - 1) * iMu + muMin;
            int index = nearest(mus, mu);
            result[iMu] = integrate(PL[index], mu, nu, N);
        }
        std::cout << std::endl;
        
        // Write result to file
        writeB(result, "result.bin");
    }
    return 0;
}
