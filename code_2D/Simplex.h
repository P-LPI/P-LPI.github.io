// simplex methods
simplex::simplex(pointC p0, pointC p1, pointC p2, pointC p3)
{
    pnts[0] = p0;
    pnts[1] = p1;
    pnts[2] = p2;
    pnts[3] = p3;
    active = true;
    
    upperleft = 0;
    lowerleft = 0;
    upperright = 0;
    lowerright = 0;
}

simplex::simplex(size_t UpperLeft, size_t LowerLeft, size_t UpperRight, size_t LowerRight)
{
    pnts[0] = pointC(0, 0);
    pnts[1] = pointC(0, 0);
    pnts[2] = pointC(0, 0);
    pnts[3] = pointC(0, 0);
    active = true;
    
    upperleft = UpperLeft;
    lowerleft = LowerLeft;
    upperright = UpperRight;
    lowerright = LowerRight;
}

void simplex::points(std::vector<cp> &points)
{
    pnts[0] = points[upperleft].p;
    pnts[1] = points[lowerleft].p;
    pnts[2] = points[upperright].p;
    pnts[3] = points[lowerright].p;
    
    active = active && points[upperleft].active && points[lowerleft].active && points[upperright].active && points[lowerright].active;
}

void simplex::subdivide(std::vector<simplex> &simplices, std::vector<cp> &points, const double delta)
{
    active = active && points[upperleft].active && points[lowerleft].active && points[upperright].active && points[lowerright].active;
    
    if(active && max(norm(points[upperleft].p  - points[lowerleft].p),
                     norm(points[upperright].p - points[lowerright].p)) > delta)
    {
        active = false;

        const pointC midleft  = (points[upperleft].p  + points[lowerleft].p)  / 2.;
        const pointC midright = (points[upperright].p + points[lowerright].p) / 2.;
        points.push_back(midleft);
        points.push_back(midright);
        simplex simp0(upperleft, points.size() - 2, upperright, points.size() - 1);
        simplex simp1(points.size() - 2, lowerleft, points.size() - 1, lowerright);

        simplices.push_back(simp0);
        simplices.push_back(simp1);
    } else if(active && max(norm(points[upperleft].p  - points[upperright].p),
                            norm(points[lowerleft].p - points[lowerright].p)) > delta)
    {
        active = false;
        
        const pointC midupper  = (points[upperleft].p  + points[upperright].p) / 2.;
        const pointC midlower  = (points[lowerleft].p  + points[lowerright].p) / 2.;
        
        points.push_back(midupper);
        points.push_back(midlower);
        simplex simp0(upperleft, lowerleft, points.size() - 2, points.size() - 1);
        simplex simp1(points.size() - 2, points.size() - 1, upperright, lowerright);
        
        simplices.push_back(simp0);
        simplices.push_back(simp1);
    }
}

std::complex<double> simplex::integrate(const pointD mu, const double nu, const int N) const
{
    if(active)
    {
        return romberg(pnts, mu, nu, N);
    } else
    {
        return 0.;
    }
}

bool simplex::isActive() const
{
    return active;
}

//pointC simplex::p0() const
//{
//    return pnts[0];
//}
//
//pointC simplex::p1() const
//{
//    return pnts[1];
//}
//
//pointC simplex::p2() const
//{
//    return pnts[2];
//}
//
//pointC simplex::p3() const
//{
//    return pnts[3];
//}
//
//size_t simplex::UL() const
//{
//    return upperleft;
//}
//
//size_t simplex::LL() const
//{
//    return lowerleft;
//}
//
//size_t simplex::UR() const
//{
//    return upperright;
//}
//
//size_t simplex::LR() const
//{
//    return lowerright;
//}

void externalVariables(std::vector<pointD> &mus, const double muMin, const double muMax, int NMu)
{
    for(int iMu = 0; iMu < NMu; iMu++)
    {
        const double muX = (muMax - muMin) / (NMu - 1) * iMu + muMin;
        
        for(int jMu = 0; jMu < NMu; jMu++)
        {
            const double muY = (muMax - muMin) / (NMu - 1) * jMu + muMin;
            mus[jMu + NMu * iMu] = pointD(muX, muY);
        }
    }
}

void initialize(std::vector<simplex> &simplices, std::vector<cp> &points, const double xMin, const double xMax, const double delta)
{
    points.push_back(cp(pointC(xMin, xMax)));
    points.push_back(cp(pointC(xMin, xMin)));
    points.push_back(cp(pointC(xMax, xMax)));
    points.push_back(cp(pointC(xMax, xMin)));
    
    simplex sim(0, 1, 2, 3);
    simplices.push_back(sim);
    
    subdevide(simplices, points, delta);
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

void flow(std::vector<simplex> &simplices, std::vector<cp> &points, const double step, const pointD mu, const double nu, const double thres, const double delta, const int Niterations)
{
    for(int i = 0; i < Niterations; i++)
    {
        flow(points, step, mu, nu, thres);
        subdevide(simplices, points, delta);
    }
}

void flow(std::vector<cp> &points, const double step, const pointD mu, const double nu, const double thres)
{
    for(int i = 0; i < points.size(); i++)
    {
        if(points[i].active)
        {
            points[i].active = h(points[i].p, mu, nu) > thres;
            if(points[i].active)
            {
                pointC g = gradient(points[i].p, mu, nu, 0.00001) / nu;
                
                if(norm(g) > 1){g = g / norm(g);}
                points[i].p = points[i].p - g * step;
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



//
