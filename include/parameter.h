#ifndef parameter_h
#define parameter_h
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <stdio.h>
#include <sstream>


struct Parameter
{
    std::string pattern;
    double chiAB;
    std::array<int, 2> degree;
    double dsmax;
    double u;
    std::array<int, 3> dim;
    std::array<double, 3> fA1;
    std::array<double, 3> fA2;
    std::array<double, 9> domain;
    std::string phiA;
    std::string phiB;
    std::string filename;

    /* constructor
     * read params from file fname
     */
    Parameter(const char * fname)
    {
        filename = fname;
    }

    void read()
    {
        std::ifstream file(filename.c_str());
        std::string line;
        while(std::getline(file, line))
        {
            auto found = line.find(":");
            if(found != std::string::npos)
            {
                auto name = line.substr(0, found);
                remove_space(name);
                std::cout << name << std::endl;
                auto data = line.substr(found+1);
                if(name == "pattern")
                {
                    pattern = data;
                    remove_space(pattern);
                }
                else if(name == "dim")
                {
                    read_array(data, dim, name);
                }
                else if(name == "dsmax" )
                {
                    read_data(data, dsmax);
                }
 /*               else if(name == "u" )
                {
                    read_data(data, u);
                }
*/
                else if(name == "degree")
                {
                    read_array(data, degree, name);
                }
                else if(name == "chiAB")
                {
                    read_data(data, chiAB);
                }
                else if(name == "fA1")
                {
                    read_array(data, fA1, name);
                }
                else if(name == "fA2")
                {
                    read_array(data, fA2, name);
                }
                else if(name == "domain")
                {
                    read_array(data, domain, name);
                }
                else if(name == "phiA")
                {
                    phiA = data;
                    remove_space(phiA);
                }
                else if(name == "phiB")
                {
                    phiB = data;
                    remove_space(phiB);
                }
            }
        }
    }

    void remove_space(std::string & str)
    {
        str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
        return;
    }

    /* read a array*/
    template<typename Array>
    void read_array(std::string & line, Array & a, std::string & name)
    { 
        std::istringstream iss(line);
        int i = 0;
        while(iss >> a[i])
        {
            i++;
            if(i == a.size())
                break;
        }
        if(i < a.size())
        {
            std::cout << "Warning! " 
                << name 
                << " parameter need " 
                << a.size() 
                << " values, but you just input " 
                << i 
                << " values "
                << std::endl;
            exit(-1);
        }
    }

    /* read one values*/
    template<typename Type>
    void read_data(std::string & line, Type & val)
    {
        std::istringstream iss(line);
        iss >> val;
    }
/*
    "pattern": "HHC",
    "chiBC": 0.30,
    "chiAC": 0.30,
    "chiAB": 0.30,
    "fA": [0.2, 0.4, 100],
    "fB": [0.2, 0.4, 100],
    "domain": [4.8, 4.8, 4.8],
    "phiA": "phiA.[30.00.30.00.30.00].[0.19.0.62.0.19]-[32]-[7.57397744]-[-0.57672435]-[4.88641].8.txt",
    "phiB": "phiB.[30.00.30.00.30.00].[0.19.0.62.0.19]-[32]-[7.57397744]-[-0.57672435]-[4.88641].8.txt",
    "phiC": "phiC.[30.00.30.00.30.00].[0.19.0.62.0.19]-[32]-[7.57397744]-[-0.57672435]-[4.88641].8.txt"
*/

    /* show the struct */
    void print()
    {
        std::cout << " degree : " << degree[0] << " " <<degree[1] << std::endl;
        std::cout << "  dsmax : " << dsmax << std::endl;
//        std::cout << "      u : " << u << std::endl;
        std::cout << "    dim : " << dim[0] << " " << dim[1] << " " << dim[2] << std::endl;
        std::cout << "pattern : " << pattern << std::endl;
        std::cout << "  chiAB : " << chiAB << std::endl;
        std::cout << "     fA1 : " << fA1[0] << " " << fA1[1] << " " << fA1[2] << std::endl;
        std::cout << "     fA2 : " << fA2[0] << " " << fA2[1] << " " << fA2[2] << std::endl;
        std::cout << " domain : " 
            << domain[0] << " "
            << domain[1] << " "
            << domain[2] << " "
            << domain[3] << " "
            << domain[4] << " "
            << domain[5] << " "
            << domain[6] << " "
            << domain[7] << " "
            << domain[8] << std::endl;
        std::cout << "   phiA : " << phiA << std::endl;
        std::cout << "   phiB : " << phiB << std::endl;
        printf("test\n");
    }
};
#endif // end of parameter_h
