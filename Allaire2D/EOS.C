#include "EOS.H"
#include<cmath>
#include<iostream>
#include<cassert>


// Sound speeds

double idealGas::calcSoundSpeed(const std::array<double,5>& arr){ // <<<<<< takes primitive
    double eps = arr[0]*(1.0/(gamma1-1))+(1-arr[0])*(1.0/(gamma2-1));
    double Y1 = arr[1]/(arr[1]+arr[2]);
    double Y2 = arr[2]/(arr[1]+arr[2]);
    double cs1 = (gamma1*arr[4]/(arr[1]/arr[0]));
    double cs2 = (gamma2*arr[4]/(arr[2]/(1-arr[0])));
    return std::sqrt((Y1*(1.0/(gamma1-1))*cs1+Y2*(1.0/(gamma2-1))*cs2)/eps);
}

double stiffenedGas::calcSoundSpeed(const std::array<double,5>& arr){ // <<<<<<
    return std::sqrt(gamma1*(arr[2]+p_inf1)/arr[0]);
}

// Getters

std::array<double,2> EOS::get_gamma(){
    return {gamma1,gamma2};
}

// Converters

// Consv: a1, rho1*a1, rho2*a2, momentum, energy
// Prim: a1, rho1, rho2, v, p

std::array<double,5> idealGas::primToConsv(const std::array<double,5>& arr)const{
    std::array<double,5> result;
    //assert(arr[0] > 0);
    result[0] = arr[0]; // alpha
    result[1] = arr[1]; // alpha1 * rho1
    result[2] = arr[2]; // alpha2 * rho2 
    result[3] = (arr[1]+arr[2])*arr[3]; // (a1r1 +a2r2)*v
    result[4] = arr[4]*((arr[0]/(gamma1-1))+(1-arr[0])/(gamma2-1)) + 0.5*(arr[1]+arr[2])*arr[3]*arr[3]; // allaire eq 31
    return result;
};

std::array<double,5> idealGas::consvToPrim(const std::array<double,5>& arr)const{
    std::array<double,5> result;
    //assert(arr[0] > 0);
    result[0] = arr[0];
    result[1] = arr[1];
    result[2] = arr[2];
    result[3] = arr[3]/(arr[1]+arr[2]);
    result[4] = (arr[4]-0.5*arr[3]*arr[3]/(arr[1]+arr[2]))/((arr[0]/(gamma1-1))+(1-arr[0])/(gamma2-1));
    //print_vect(result);
    return result;
}

std::array<double,5> stiffenedGas::primToConsv(const std::array<double,5>& arr)const{
    std::array<double,5> result;
    assert(arr[0] != 0);
    result[0] = arr[0]; // rho
    result[1] = arr[0]*arr[1]; // rho*v
    result[2] = ((arr[2]+gamma1*p_inf1)/(gamma1 - 1)) + (0.5*arr[0]*arr[1]*arr[1]); // E
    return result;
};

std::array<double,5> stiffenedGas::consvToPrim(const std::array<double,5>& arr)const{
    std::array<double,5> result;
    assert(arr[0] != 0);
    result[0] = arr[0];
    result[1] = arr[1]/arr[0];
    result[2] = (arr[2]-(arr[1]*arr[1]/(2.0*arr[0])))*(gamma1-1)-gamma1*p_inf1;
    //print_vect(result);
    return result;
}

// constructors

EOS::EOS(double set_gamma1, double set_gamma2){
    gamma1 = set_gamma1;
    gamma2 = set_gamma2;
}

idealGas::idealGas(double set_gamma1,double set_gamma2) : EOS(set_gamma1,set_gamma2){
}

stiffenedGas::stiffenedGas(double set_gamma1, double set_gamma2, double set_p_inf1, double set_p_inf2, double set_e_inf1, double set_e_inf2) : EOS(set_gamma1,set_gamma2){
    p_inf1 = set_p_inf1;
    p_inf2 = set_p_inf2;
    e_inf1 = set_e_inf1;
    e_inf2 = set_e_inf2;
}