double h(const pointC p, const pointD mu, const double nu)
{
    return real(func(p, mu, nu));
}

double H(const pointC p, const pointD mu, const double nu)
{
    return imag(func(p, mu, nu));
}

pointC gradient(const pointC p, const pointD mu, const double nu, const double epsilon)
{
    const pointC epsilonX(epsilon, 0);
    const pointC epsilonY(0, epsilon);
    
    return pointC((h(p + epsilonX, mu, nu) - h(p - epsilonX, mu, nu)) + I * (h(p + epsilonX * I, mu, nu) - h(p - epsilonX * I, mu, nu)),
                  (h(p + epsilonY, mu, nu) - h(p - epsilonY, mu, nu)) + I * (h(p + epsilonY * I, mu, nu) - h(p - epsilonY * I, mu, nu))) / (2. * epsilon);
}
