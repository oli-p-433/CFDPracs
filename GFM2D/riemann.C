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
    std::array<double,4> sL = {1,0,0,1};//{1000,0,1e9};
    double gamma1 = 1.4;

    std::array<double,4> sR = {0.125,0,0,0.1}; // {50,0,1e5};
    double gamma2 = 1.67;

    double time = 0.0002;
    double discPosition = 0.5;

    riemann solution(gamma1,gamma2,sL,sR,1,-1,1,0,1,time,0,0); // need p_inf setup

    std::array<double,4> resultL = solution.interfaceRiemann(1);
    std::array<double,4> resultR = solution.interfaceRiemann(0);

    std::cout << "Left interface state:" << std::endl;
    for (size_t i = 0; i<resultL.size(); ++i){
        std::cout << resultL[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Right interface state:" << std::endl; 
    for (size_t i = 0; i<resultR.size(); ++i){
        std::cout << resultR[i] << " ";
    }
    std::cout << std::endl;


    /*
    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    solution.dirName = dirname;

    // initialise data locations
    std::vector<std::string> variables = {"rho","v","p","rhoExact","vExact","pExact"};
    for (size_t i = 0; i<variables.size(); ++i){
        std::string name = dirname + "/" + variables[i];
        std::filesystem::create_directory(name);
    }
        */

    // Exact Riemann solver
}