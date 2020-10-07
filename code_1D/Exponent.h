double h(const std::complex<double> x, const double mu, const double nu)
{
    return real(func(x, mu, nu));
}

double H(const std::complex<double> x, const double mu, const double nu)
{
    return imag(func(x, mu, nu));
}

std::complex<double> gradient(const std::complex<double> x, const double mu, const double nu, const double epsilon)
{
        return ((h(x + epsilon, mu, nu) - h(x - epsilon, mu, nu)) + I * (h(x + I * epsilon, mu, nu) - h(x - I * epsilon, mu, nu))) / (2. * epsilon);
}
