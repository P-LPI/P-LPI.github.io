std::complex<double> integrate(const std::vector<simplex> &simplices, const pointD mu, const double nu, const int N)
{
    std::complex<double> sum = 0;
    for(int index = 0; index < simplices.size(); index++)
    {
        if(simplices[index].isActive())
        {
            sum = sum + simplices[index].integrate(mu, nu, N);
        }
    }

    return sum;
}

std::complex<double> integrate(std::vector<simplex> &simplices, std::vector<cp> &points, const pointD mu, const double nu, const int N)
{
    importPoints(simplices, points);
    return integrate(simplices, mu, nu, N);
}

pointC map(const pointD p,  const pointC pnts[4])
{
    return (pnts[0] * (1. - p.x) * p.y +
            pnts[1] * (1. - p.x) * (1. - p.y) +
            pnts[2] * p.x * p.y +
            pnts[3] * p.x * (1. - p.y));
}

std::complex<double> jacob(const pointD p, const pointC pnts[4])
{
    std::complex<double> hess[2][2];
    hess[0][0] = (pnts[3].x - pnts[1].x) * (1. - p.y) + (pnts[2].x - pnts[0].x) * p.y;
    hess[0][1] = (pnts[3].y - pnts[1].y) * (1. - p.y) + (pnts[2].y - pnts[0].y) * p.y;
    hess[1][0] = (pnts[0].x - pnts[1].x) * (1. - p.x) + (pnts[2].x - pnts[3].x) * p.x;
    hess[1][1] = (pnts[0].y - pnts[1].y) * (1. - p.x) + (pnts[2].y - pnts[3].y) * p.x;

    return hess[0][0] * hess[1][1] - hess[0][1] * hess[1][0];
}

std::complex<double> integrand(const double u, const double v, const pointC pnts[4], const pointD &mu, const double nu)
{
    return jacob(pointD(u,v), pnts) * std::exp(func(map(pointD(u,v), pnts), mu, nu));
}

std::complex<double> romberg(const pointC pnts[4], const pointD &mu, const double nu, const int N)
{
    double h[N+1];
    std::complex<double> r[N+1][N+1];
    for (int i = 1; i < N + 1; ++i) {
        h[i] = 1. / pow(2, i - 1);
    }

    r[1][1] = h[1] * h[1] * (integrand(0,0, pnts, mu, nu) + integrand(h[1], 0, pnts, mu, nu) + integrand(0, h[1], pnts, mu, nu) + integrand(h[1], h[1], pnts, mu, nu)) / 4.;

    for (int i = 2; i < N + 1; ++i) {
        std::complex<double> coeff = 0;
        for (int k = 1; k <= pow(2, i - 2); ++k) {
            coeff += 2. * (integrand((2 * k - 1) * h[i], 0, pnts, mu, nu) +
                           integrand((2 * k - 1) * h[i], h[1], pnts, mu, nu) +
                           integrand(0, (2 * k - 1) * h[i], pnts, mu, nu) +
                           integrand(h[1], (2 * k - 1) * h[i], pnts, mu, nu));
        }
        for (int k = 1; k <= pow(2, i - 2); ++k) {
            for(int l = 1; l < pow(2, i - 1); ++l) {
                coeff += 4. * integrand((2 * k - 1) * h[i], h[i] * l, pnts, mu, nu);
            }
        }
        for (int k = 1; k <= pow(2, i - 2) - 1; ++k) {
            for (int l = 1; l <= pow(2, i - 2); ++l) {
                coeff += 4. * integrand(2. * h[i] * k, (2 * l - 1) * h[i], pnts, mu, nu);
            }
        }

        r[i][1] = 0.25 * (r[i - 1][1] + h[i] * h[i] * coeff);
    }

    for (int i = 2; i < N + 1; ++i) {
        for (int j = 2; j <= i; ++j) {
            r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
        }
    }
    return r[N][N];
    return 0;
}
