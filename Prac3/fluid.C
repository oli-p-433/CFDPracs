#include "fluid.H" 
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>

// constructor

fluid::fluid(int nCells, int nGhost){
    u.resize(nCells+2*nGhost);
    fluxes.resize(nCells+1);
    halfSlopes.resize(nCells+3);
    r.resize(nCells+2); // r can only be defined for inner cells
    uBarL.resize(nCells+2);
    uBarR.resize(nCells+2);
    uBarLupd.resize(nCells+2);
    uBarRupd.resize(nCells+2);

    slopeWeight = 1;
};

// Misc

void fluid::print_arr(std::vector< std::array<double,3> > arr, int var){
    for (std::vector<double>::size_type i = 0; i < arr.size(); ++i) {
        std::cout << arr[i][var] << " ";
    }
    std::cout << std::endl;
}

void fluid::print_vect(std::array<double,3> v){
    for (int j = 0; j<3; ++j){
        std::cout << v[j] << " ";
    }
    std::cout << std::endl;
}

