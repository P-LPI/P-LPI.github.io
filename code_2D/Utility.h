int nearest(const std::vector<pointD> &list, const pointD &p)
{
    int index = 0;
    double min = norm(list[0] - p);
    for(int i = 1; i < list.size(); i++)
    {
        if(norm(list[i] - p) < min)
        {
            index = i;
            min = norm(list[i] - p);
        }
    }
    return index;
}

double min(const double a, const double b)
{
    if(a < b)
    {
        return a;
    } else
    {
        return b;
    }
}

double min(const double a, const double b, const double c)
{
    return min(a, min(b, c));
}

double max(const double a, const double b)
{
    if(a > b)
    {
        return a;
    } else
    {
        return b;
    }
}

//void writeB(std::vector<std::vector<simplex> > &PL, std::vector<pointD> mus, std::string fileName, std::string fileNameLog, std::string fileNameMus)
//{
//    std::ofstream file; file.open(fileName, std::ios::binary);
//    std::ofstream fileLog; fileLog.open(fileNameLog);
//    std::ofstream fileMus; fileMus.open(fileNameMus);
//    if(file.is_open() && fileLog.is_open() && fileMus.is_open())
//    {
//        for(int i = 0; i < PL.size(); i++)
//        {
//            std::vector<simplex> simplices = PL[i];
//            fileLog << int(simplices.size()) << std::endl;
//            writeB(mus[i], fileMus);
//            writeB(simplices, file);
//        }
//        
//        file.close();
//        fileLog.close();
//        fileMus.close();
//    } else
//    {
//        std::cout << "Could not open " << fileName << " or " << fileNameLog << " or " << fileNameMus << std::endl;
//    }
//}
//
//void writeB(std::vector<simplex> &simplices, std::string fileName)
//{
//    std::ofstream file; file.open(fileName, std::ios::binary);
//    if(file.is_open())
//    {
//        writeB(simplices, file);
//        file.close();
//    } else
//    {
//        std::cout << "Could not open " << fileName << std::endl;
//    }
//}
//
//void writeB(std::vector<simplex> &simplices, std::vector<cp> &points, std::string fileName)
//{
//    importPoints(simplices, points);
//    std::ofstream file; file.open(fileName, std::ios::binary);
//    if(file.is_open())
//    {
//        writeB(simplices, file);
//        file.close();
//    } else
//    {
//        std::cout << "Could not open " << fileName << std::endl;
//    }
//}
//
//void writeB(std::vector<simplex> &simplices, std::ofstream &file)
//{
//    for(int index = 0; index < simplices.size(); index++)
//    {
//        if(simplices[index].isActive())
//        {
//            writeB(simplices[index], file);
//        }
//    }
//}
//
//void writeB(simplex &simp, std::ofstream &file)
//{
//    pointC p0 = simp.p0();
//    pointC p1 = simp.p1();
//    pointC p2 = simp.p2();
//    pointC p3 = simp.p3();
//
//    writeB(p0, file);
//    writeB(p1, file);
//    writeB(p2, file);
//    writeB(p3, file);
//}
//
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
//
//void writeB(const pointC &p, std::ofstream &file)
//{
//    writeB(p.x, file);
//    writeB(p.y, file);
//}
//
//void writeB(const pointD &p, std::ofstream &file)
//{
//    writeB(p.x, file);
//    writeB(p.y, file);
//}

void writeB(std::complex<double> number, std::ofstream &file)
{
    double r = real(number);
    double i = imag(number);
    file.write((char*) &r, sizeof(double));
    file.write((char*) &i, sizeof(double));
}

//// Import
////void read(std::string fileName, std::string fileNameLog, std::string fileNameMus, std::vector<std::vector<simplex> > &PL, std::vector<point> &mus)
////{
////    std::ifstream fileLog; fileLog.open(fileNameLog);
////    std::vector<int> length;
////    if(fileLog.is_open())
////    {
////        int number;
////        while (fileLog >> number) {
////            length.push_back(number);
////        }
////        fileLog.close();
////    } else
////    {
////        std::cout << "Could not open: " + fileNameLog << std::endl;
////        exit (EXIT_FAILURE);
////    }
////
////    std::ifstream fileMus; fileMus.open(fileNameMus);
////    if(fileMus.is_open())
////    {
////        for(int i = 0; i < length.size(); i++)
////        {
////            mus.push_back(readB(fileMus));
////        }
////    } else
////    {
////        std::cout << "Could not open: " + fileNameMus << std::endl;
////        exit (EXIT_FAILURE);
////    }
////
////    std::ifstream file; file.open(fileName, std::ios::binary);
////    if(file.is_open())
////    {
////        for(int i = 0; i < length.size(); i++)
////        {
////            std::vector<simplex> simplices;
////            for(int j = 0; j < length[i]; j++)
////            {
////                double pxr, pxi, pyr, pyi;
////                file.read((char*) &pxr, sizeof(double));
////                file.read((char*) &pxi, sizeof(double));
////                file.read((char*) &pyr, sizeof(double));
////                file.read((char*) &pyi, sizeof(double));
////                simplices.push_back(simplex(pxr + I * pxi, pyr + I * pyi));
////            }
////            PL.push_back(simplices);
////        }
////    } else
////    {
////        std::cout << "Could not open: " + fileName << std::endl;
////        exit (EXIT_FAILURE);
////    }
////}
////
////point readB(std::ifstream &file)
////{
////    double px, py, pz;
////    file.read((char*) &px, sizeof(double));
////    file.read((char*) &py, sizeof(double));
////    file.read((char*) &pz, sizeof(double));
////
////    point p(px, py, pz);
////
////    return p;
////}
//
