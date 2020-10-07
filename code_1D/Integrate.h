std::complex<double> integrate(const std::vector<simplex> &simplices, const double mu, const double nu, const int N)
{
    std::complex<double> sum = 0;
    for(int index = 0; index < simplices.size(); index++)
    {
        sum = sum + simplices[index].integrate(mu, nu, N);
    }

    return sum;
}

std::complex<double> integrate(std::vector<simplex> &simplices, std::vector<cp> &points, const double mu, const double nu, const int N)
{
    importPoints(simplices, points);
    return integrate(simplices, mu, nu, N);
}

std::complex<double> romberg(std::complex<double> a, std::complex<double> b, const double mu, const double nu, const int N)
{
    std::complex<double> h[N+1], r[N+1][N+1];
    for (int i = 1; i < N + 1; ++i) {
        h[i] = (b - a) / pow(2, i - 1);
    }
    r[1][1] = h[1] / 2. * (std::exp(func(a, mu, nu)) + std::exp(func(b, mu, nu)));
    for (int i = 2; i < N + 1; ++i) {
        std::complex<double> coeff = 0;
        for (int k = 1; k <= pow(2, i - 2); ++k) {
            coeff += std::exp(func(a + (2. * k - 1.) * h[i], mu, nu));
        }
        r[i][1] = 0.5 * (r[i - 1][1] + h[i - 1] * coeff);
    }

    for (int i = 2; i < N + 1; ++i) {
        for (int j = 2; j <= i; ++j) {
            r[i][j] = r[i][j - 1] + (r[i][j - 1] - r[i - 1][j - 1]) / (pow(4, j - 1) - 1);
        }
    }
    return r[N][N];
}
