// simplex methods
simplex::simplex(std::complex<double> p0, std::complex<double> p1)
{
    pnts[0] = p0;
    pnts[1] = p1;
    active = true;
    
    left = 0;
    right = 0;
}

simplex::simplex(size_t Left, size_t Right)
{
    pnts[0] = 0;
    pnts[1] = 0;
    active = true;
    
    left = Left;
    right = Right;
}

void simplex::points(std::vector<cp> &points)
{
    pnts[0] = points[left].p;
    pnts[1] = points[right].p;
    
    active = points[left].active && points[right].active;
}

void simplex::subdivide(std::vector<simplex> &simplices, std::vector<cp> &points, const double delta)
{
    active = active && points[left].active && points[right].active;
    
    if(active && abs(points[left].p - points[right].p) > delta)
    {
        active = false;
        
        std::complex<double> midpoint = (points[left].p + points[right].p) / 2.;
        points.push_back(midpoint);
        simplex simp0(left, points.size() - 1);
        simplex simp1(points.size() - 1, right);
        
        simplices.push_back(simp0);
        simplices.push_back(simp1);
    }
}

std::complex<double> simplex::integrate(const double mu, const double nu, const int N) const
{
    if(active)
    {
        return romberg(pnts[0], pnts[1], mu, nu, N);
    } else
    {
        return 0.;
    }
}

bool simplex::isActive()
{
    return active;
}

//std::complex<double> simplex::p0()
//{
//    return pnts[0];
//}
//
//std::complex<double> simplex::p1()
//{
//    return pnts[1];
//}
//
//size_t simplex::L()
//{
//    return left;
//}
//
//size_t simplex::R()
//{
//    return right;
//}

void externalVariables(std::vector<double> &mus, const double muMin, const double muMax, int NMu)
{
    for(int iMu = 0; iMu < NMu; iMu++)
    {
        const double mu = (muMax - muMin) / (NMu - 1) * iMu + muMin;
        mus[iMu] = mu;
    }
}

void initialize(std::vector<simplex> &simplices, std::vector<cp> &points, const double xMin, const double xMax, const double delta)
{
    points.push_back(cp(xMin));
    points.push_back(cp(xMax));
    
    simplex sim(0, 1);
    simplices.push_back(sim);
    
    subdevide(simplices, points, delta);
    clean(simplices);
}

void subdevide(std::vector<simplex> &simplices, std::vector<cp> &points, const double delta)
{
    for(int index = 0; index < simplices.size(); index++)
      {
          simplices[index].subdivide(simplices, points, delta);
      }
      clean(simplices);
}

void clean(std::vector<simplex> &simplices)
{
    std::vector<simplex> newSimplices;
    for(int index = 0; index < simplices.size(); index++)
    {
        simplex element = simplices[index];
        if(element.isActive())
        {
            newSimplices.push_back(element);
        }
    }
    simplices = newSimplices;
}

void flow(std::vector<simplex> &simplices, std::vector<cp> &points, const double step, const double mu, const double nu, const double thres, const double delta, const int Niterations)
{
    for(int i = 0; i < Niterations; i++)
    {
        flow(points, step, mu, nu, thres);
        subdevide(simplices, points, delta);
        clean(simplices);
    }
}

void flow(std::vector<cp> &points, const double step, const double mu, const double nu, const double thres)
{
    for(int i = 0; i < points.size(); i++)
    {
        if(points[i].active)
        {
            points[i].active = h(points[i].p, mu, nu) > thres;
            if(points[i].active)
            {
                std::complex<double> g = gradient(points[i].p, mu, nu, 0.00001) / nu;
                
                if(abs(g) > 1){g = g / abs(g);}
                points[i].p = points[i].p - step * g;
            }
        }
    }
}

void importPoints(std::vector<simplex> &simplices, std::vector<cp> &points)
{
    for(int i = 0; i < simplices.size(); i++)
    {
        simplices[i].points(points);
    }
}
