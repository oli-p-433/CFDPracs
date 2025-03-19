#include "EOS.H"
#include <cmath>
#include <iostream>
#include <cassert>


// Sound speeds

double idealGas::calcSoundSpeed(const std::array<double,4>& arr){ // <<<<<< takes primitive
    return std::sqrt(gamma*arr[3]/arr[0]);
}

double stiffenedGas::calcSoundSpeed(const std::array<double,4>& arr){ // <<<<<<
    return std::sqrt(gamma*(arr[3]+p_inf)/arr[0]);
}

// Getters

double EOS::get_gamma(){
    return gamma;
}



// Converters

// Consv: a1, rho1*a1, rho2*a2, momentum, energy
// Prim: a1, rho1, rho2, v, p

std::array<double,4> idealGas::primToConsv(const std::array<double,4>& arr)const{
    std::array<double,4> result;
    //assert(arr[0] > 0);
    result[RHO] = arr[RHO];
    result[XMOM] = arr[RHO]*arr[UX]; //rho*vx
    result[YMOM] = arr[RHO]*arr[UY]; // rho*vy
    result[ENE] = (arr[PRES]/(gamma - 1)) + (arr[RHO]*(arr[UX]*arr[UX]+arr[UY]*arr[UY])) / 2.0; // E modify for total

    return result;
};

std::array<double,4> idealGas::consvToPrim(const std::array<double,4>& arr)const{
    std::array<double,4> result;
    //assert(arr[0] > 0);
    result[RHO] = arr[RHO];
    result[UX] = arr[XMOM]/arr[RHO];
    result[UY] = arr[YMOM]/arr[RHO];
    result[PRES] = (gamma - 1)*(arr[ENE] - (arr[XMOM]*arr[XMOM]+arr[YMOM]*arr[YMOM]) / (2 * arr[RHO]) ); // modify for total
    //print_vect(result);
    return result;
}

std::array<double,4> stiffenedGas::primToConsv(const std::array<double,4>& arr)const{
    std::array<double,4> result;
    //assert(arr[RHO] > 0); assert(arr[PRES] > 0);
    result[RHO] = arr[RHO]; // rho
    result[XMOM] = arr[RHO]*arr[XMOM]; // rho*vx
    result[YMOM] = arr[RHO]*arr[YMOM]; // rho*vy
    result[ENE] = ((arr[PRES]+gamma*p_inf)/(gamma - 1)) + (arr[RHO]*(arr[UX]*arr[UX]+arr[UY]*arr[UY])) / 2.0; // E modify for total
    return result;
}

std::array<double,4> stiffenedGas::consvToPrim(const std::array<double,4>& arr)const{
    std::array<double,4> result;
    //assert(arr[RHO] > 0);
    result[RHO] = arr[RHO];
    result[UX] = arr[XMOM]/arr[RHO];
    result[UY] = arr[YMOM]/arr[RHO];
    result[PRES] = (gamma - 1)*(arr[ENE] - (arr[XMOM]*arr[XMOM]+arr[YMOM]*arr[YMOM]) / (2 * arr[RHO]) ) - gamma*p_inf; // modify for total
    //print_vect(result);
    //assert(result[RHO] > 0);
    return result;
}

// constructors

EOS::EOS(double set_gamma){
    gamma = set_gamma;

}

idealGas::idealGas(double set_gamma) : EOS(set_gamma){
}

stiffenedGas::stiffenedGas(double set_gamma, double set_p_inf, double set_e_inf) : EOS(set_gamma){
    p_inf = set_p_inf;
    e_inf = set_e_inf;

}