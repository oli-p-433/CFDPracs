#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <array>

#include "riem.H"

int main(){

    // set states
    std::array<double,3> sL = {50,1000,1e5};//{1000,0,1e9};
    double gamma1 = 1.4;

    std::array<double,3> sR = {1000,1000,1e5}; // {50,0,1e5};
    double gamma2 = 4.4;

    double time = 250e-6;
    double discPosition = 0.5;

    riemann solution(gamma1,gamma2,sL,sR,0,1,discPosition,time,500,0,0);

    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    solution.dirName = dirname;

    // initialise data locations
    std::array<std::string,6> variables = {"rho","v","p","rhoExact","vExact","pExact"};
    for (size_t i = 0; i<variables.size(); ++i){
        std::string name = dirname + "/" + variables[i];
        std::filesystem::create_directory(name);
    }

    // Exact Riemann solver
    double pStar = solution.exctRiemann()[2];
    solution.exactRiemannSolution(pStar,discPosition);
}