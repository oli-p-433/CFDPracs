#include "EOS.H"
#include<cmath>
#include<iostream>
#include<cassert>


// Sound speeds

double idealGas::calcSoundSpeed(const std::array<double,3> arr){
    return std::sqrt(gamma*arr[2]/arr[0]);
}

double stiffenedGas::calcSoundSpeed(const std::array<double,3> arr){
    return std::sqrt(gamma*(arr[2]+p_inf)/arr[0]);
}

// Converters

std::array<double,3> idealGas::primToConsv(const std::array<double,3> arr)const{
    std::array<double,3> result;
    assert(arr[0] != 0);
    result[0] = arr[0]; // rho
    result[1] = arr[0]*arr[1]; // rho*v
    result[2] = (arr[2]/(gamma - 1)) + (arr[0]*arr[0]*arr[1]*arr[1]) / (2 * arr[0]); // E
    return result;
};

std::array<double,3> idealGas::consvToPrim(const std::array<double,3> arr)const{
    std::array<double,3> result;
    assert(arr[0] != 0);
    result[0] = arr[0];
    result[1] = arr[1]/arr[0];
    result[2] = (gamma - 1)*(arr[2] - (arr[1]*arr[1]) / (2 * arr[0]) );
    //print_vect(result);
    return result;
}

std::array<double,3> stiffenedGas::primToConsv(const std::array<double,3> arr)const{
    std::array<double,3> result;
    assert(arr[0] != 0);
    result[0] = arr[0]; // rho
    result[1] = arr[0]*arr[1]; // rho*v
    result[2] = ((arr[2]+gamma*p_inf)/(gamma - 1)) + (0.5*arr[0]*arr[1]*arr[1]); // E
    return result;
};

std::array<double,3> stiffenedGas::consvToPrim(const std::array<double,3> arr)const{
    std::array<double,3> result;
    assert(arr[0] != 0);
    result[0] = arr[0];
    result[1] = arr[1]/arr[0];
    result[2] = (arr[2]-(arr[1]*arr[1]/(2.0*arr[0])))*(gamma-1)-gamma*p_inf;
    //print_vect(result);
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