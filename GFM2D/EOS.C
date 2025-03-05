#include "EOS.H"
#include<cmath>
#include<iostream>
#include<cassert>


// Sound speeds

double idealGas::calcSoundSpeed(const std::array<double,4> arr){ // takes primitive variables
    assert(arr[RHO] > 0); assert(arr[PRES] > 0);
    return std::sqrt(gamma*arr[PRES]/arr[RHO]);
}

double stiffenedGas::calcSoundSpeed(const std::array<double,4> arr){
    assert(arr[RHO] > 0); assert(arr[PRES] > 0);
    return std::sqrt(gamma*(arr[PRES]+p_inf)/arr[RHO]);
}

// Getters

double EOS::get_gamma(){
    return gamma;
}

// Converters

std::array<double,4> idealGas::primToConsv(std::array<double,4> arr)const{
    std::array<double,4> result;
    assert(arr[RHO] > 0); assert(arr[PRES] > 0);
    result[RHO] = arr[RHO]; // rho
    result[XMOM] = arr[RHO]*arr[XMOM]; // rho*vx
    result[YMOM] = arr[RHO]*arr[YMOM]; // rho*vx
    result[ENE] = (arr[PRES]/(gamma - 1)) + (arr[RHO]*(arr[UX]*arr[UX]+arr[UY]*arr[UY])) / 2.0; // E modify for total
    return result;
}

std::array<double,4> idealGas::consvToPrim(std::array<double,4> arr)const{
    std::array<double,4> result;
    assert(arr[RHO] > 0);
    result[RHO] = arr[RHO];
    result[UX] = arr[XMOM]/arr[RHO];
    result[UY] = arr[YMOM]/arr[RHO];
    result[PRES] = (gamma - 1)*(arr[ENE] - (arr[XMOM]*arr[XMOM]+arr[YMOM]*arr[YMOM]) / (2 * arr[RHO]) ); // modify for total
    //print_vect(result);
    assert(result[PRES] > 0);
    assert(result[RHO] > 0);
    return result;
}

std::array<double,4> stiffenedGas::primToConsv(std::array<double,4> arr)const{
    std::array<double,4> result;
    assert(arr[RHO] > 0); assert(arr[PRES] > 0);
    result[RHO] = arr[RHO]; // rho
    result[XMOM] = arr[RHO]*arr[XMOM]; // rho*vx
    result[YMOM] = arr[RHO]*arr[YMOM]; // rho*vy
    result[ENE] = ((arr[PRES]+gamma*p_inf)/(gamma - 1)) + (arr[RHO]*(arr[UX]*arr[UX]+arr[UY]*arr[UY])) / 2.0; // E modify for total
    return result;
}

std::array<double,4> stiffenedGas::consvToPrim(std::array<double,4> arr)const{
    std::array<double,4> result;
    assert(arr[RHO] > 0);
    result[RHO] = arr[RHO];
    result[UX] = arr[XMOM]/arr[RHO];
    result[UY] = arr[YMOM]/arr[RHO];
    result[PRES] = (gamma - 1)*(arr[ENE] - (arr[XMOM]*arr[XMOM]+arr[YMOM]*arr[YMOM]) / (2 * arr[RHO]) ) - gamma*p_inf; // modify for total
    //print_vect(result);
    assert(result[PRES] > 0);
    assert(result[RHO] > 0);
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