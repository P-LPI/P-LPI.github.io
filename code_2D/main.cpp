#include "Header.h"

std::complex<double> func(const pointC p, const pointD mu, const double nu)
{
    return I * nu * (pow(p.x - mu.x, 2.) / 2. + pow(p.y - mu.y, 2.) / 2. + phi(p));
}

std::complex<double> phi(const pointC p)
{
    return 0.65 / (1. + p.x * p.x + 2. * p.y * p.y);
}

int main(int argc, const char * argv[])
{
    const double nu = 100.;
    
    // Parameters
    const double xMin = -5.;
    const double xMax = +5.;
    const double delta = 0.1;
    
    const double step = 0.02;
    const double thres = -20.;
    
    const int Niterations = 20;
    
    const double muMin = -1.;
    const double muMax = +1.;

    const int NMu = pow(2, 3);
    
    const int N = 3;
    const int NM = pow(2, 7);
    
    {
        // Initial integration domain
        std::vector<simplex> simplicesInit;
        std::vector<cp> pointsInit;
        initialize(simplicesInit, pointsInit, xMin, xMax, delta);
        
        // Compute the Picard-Lefschetz thimbles
        std::vector<pointD> mus(NMu * NMu);
        externalVariables(mus, muMin, muMax, NMu);
        
        std::cout << "Picard-Lefschetz thimble:";
        std::vector<std::vector<simplex>> PL(mus.size());
#pragma omp parallel for
        for(int index = 0; index < mus.size(); index++)
        {
            if((index % (NMu * NMu / 16)) == 0){std::cout << "."; std::cout.flush();}
            std::vector<simplex> simplices = simplicesInit;
            std::vector<cp> points = pointsInit;
            flow(simplices, points, step, mus[index], nu, thres, delta, Niterations);
            importPoints(simplices, points);
            PL[index] = simplices;
        }
        std::cout << std::endl;

        // Evaluate the integral
        std::cout << "Integrate along thimbles:";
        std::vector<std::complex<double>> result(NM * NM);
#pragma omp parallel for
        for(int iMu = 0; iMu < NM; iMu++)
        {
            if((iMu % (NM / 16)) == 0){std::cout << "."; std::cout.flush();}
            for(int jMu = 0; jMu < NM; jMu++)
            {
                const pointD mu((muMax - muMin) / (NM - 1) * iMu + muMin,
                                (muMax - muMin) / (NM - 1) * jMu + muMin);
                int index = nearest(mus, mu);
                result[jMu + iMu * NM] = integrate(PL[index], mu, nu, N);
            }
        }
        std::cout << std::endl;

        // Write result to file
        writeB(result, "result.bin");
    }
    return 0;
}
