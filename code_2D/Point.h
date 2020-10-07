template<typename coordinate_type>
struct point {
    coordinate_type x;
    coordinate_type y;
    
    point()
    {
        x = 0.;
        y = 0.;
    }
    
    point(coordinate_type X, coordinate_type Y)
    {
        x = X;
        y = Y;
    }
    
    point& operator=(const point& a)
    {
        x = a.x;
        y = a.y;
        return *this;
    }

    point operator+(const point& a) const
    {
        return point(x + a.x, y + a.y);
    }
    
    point operator-(const point& a) const
    {
        return point(x - a.x, y - a.y);
    }
    
    point operator*(const coordinate_type& a) const
    {
        return point(a * x, a * y);
    }
    
    coordinate_type operator*(const point& a) const
    {
        return a.x * x + a.y * y;
    }
    
    point operator/(const double& a) const
    {
        return point(x / a, y / a);
    }
    
    bool operator==(const point& a) const
    {
        return (x == a.x && y == a.y);
    }
};

typedef point<double> pointD;
typedef point<std::complex<double>> pointC;

std::ostream &operator<<(std::ostream &os, pointD const &p) {
    return os << p.x << ", " << p.y;
}

std::ostream &operator<<(std::ostream &os, pointC const &p) {
    return os << p.x << ", " << p.y;
}

double norm(const pointD p)
{
    return sqrt(p.x * p.x + p.y * p.y);
}

double norm(const pointC p)
{
    return real(sqrt(p.x * conj(p.x) + p.y * conj(p.y)));
}
