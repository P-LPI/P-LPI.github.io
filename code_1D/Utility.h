int nearest(const std::vector<double> &list, const double &p)
{
    int index = 0;
    double min = fabs(list[0] - p);
    for(int i = 1; i < list.size(); i++)
    {
        if(fabs(list[i] - p) < min)
        {
            index = i;
            min = fabs(list[i] - p);
        }
    }
    return index;
}

void writeB(std::vector<std::complex<double> > &psi, std::string fileName)
{
    std::ofstream file; file.open(fileName, std::ios::binary);
    if(file.is_open())
    {
        for(int index = 0; index < psi.size(); index++)
        {
            writeB(psi[index], file);
        }
    } else
    {
        std::cout << "Could not open " << fileName << std::endl;
    }
}

void writeB(std::complex<double> number, std::ofstream &file)
{
    double r = real(number);
    double i = imag(number);
    file.write((char*) &r, sizeof(double));
    file.write((char*) &i, sizeof(double));
}
