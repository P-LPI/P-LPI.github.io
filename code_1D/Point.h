struct point {
    double x;
    double y;
    double z;
    
    point()
    {
        x = 0.;
        y = 0.;
        z = 0.;
    }
    
    point(double X, double Y, double Z)
    {
        x = X;
        y = Y;
        z = Z;
    }
    
    point& operator=(const point& a)
    {
        x = a.x;
        y = a.y;
        z = a.z;
        return *this;
    }

    point operator+(const point& a) const
    {
        return point(x + a.x, y + a.y, z + a.z);
    }
    
    point operator-(const point& a) const
    {
        return point(x - a.x, y - a.y, z - a.z);
    }
    
    point operator*(const double& a) const
    {
        return point(a * x, a * y, a * z);
    }
    
    std::complex<double> operator*(const point& a) const
    {
        return a.x * x + a.y * y + a.z * z;
    }
    
    point operator/(const double& a) const
    {
        return point(x / a, y / a, z / a);
    }
    
    bool operator==(const point& a) const
    {
        return (x == a.x && y == a.y && z == a.z);
    }
};

std::ostream &operator<<(std::ostream &os, point const &p) {
    return os << p.x << ", " << p.y << ", " << p.z;
}

double norm(const point p)
{
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}
