#include "EOS.H"
#include<cmath>
#include<iostream>
#include<cassert>


// Sound speeds

double idealGas::calcSoundSpeed(const std::array<double,7> arr){
    return std::sqrt(gamma*arr[PRES]/arr[RHO]);
}

double stiffenedGas::calcSoundSpeed(const std::array<double,7> arr){
    return std::sqrt(gamma*(arr[PRES]+p_inf)/arr[RHO]);
}

// Converters

std::array<double,7> idealGas::primToConsv(const std::array<double,7> arr)const{
    std::array<double,7> result;
    assert(arr[RHO] != 0);
    result[RHO] = arr[RHO]; // rho
    result[XMOM] = arr[RHO]*arr[UX]; // rho*v
    result[YMOM] = arr[RHO]*arr[UY]; // rho*v
    result[ZMOM] = arr[RHO]*arr[UZ]; // rho*v
    result[ENE] = (arr[PRES]/(gamma - 1)) + 0.5*arr[RHO]*(arr[UX]*arr[UX]+arr[UY]*arr[UY]+arr[UZ]*arr[UZ])+0.5*(BX*BX+arr[BY]*arr[BY]+arr[BZ]*arr[BZ]); // E
    result[BY] = arr[BY];
    result[BZ] = arr[BZ];

    return result;
};

std::array<double,7> idealGas::consvToPrim(const std::array<double,7> arr)const{
    std::array<double,7> result;
    assert(arr[RHO] != 0);
    result[RHO] = arr[RHO];
    result[UX] = arr[XMOM]/arr[RHO];
    result[UY] = arr[YMOM]/arr[RHO];
    result[UZ] = arr[ZMOM]/arr[RHO];
    result[PRES] = (gamma - 1)*(arr[ENE] - (arr[XMOM]*arr[XMOM]+arr[YMOM]*arr[YMOM]+arr[ZMOM]*arr[ZMOM]) / (2 * arr[RHO]) - 0.5*(BX*BX+arr[BY]*arr[BY]+arr[BZ]*arr[BZ]));
    result[BY] = arr[BY];
    result[BZ] = arr[BZ];

    return result;
}

std::array<double,7> stiffenedGas::primToConsv(const std::array<double,7> arr)const{
    std::array<double,7> result;
    assert(arr[0] != 0);
    result[0] = arr[0]; // rho
    result[1] = arr[0]*arr[1]; // rho*v
    result[2] = ((arr[2]+gamma*p_inf)/(gamma - 1)) + (0.5*arr[0]*arr[1]*arr[1]); // E
    return result;
};

std::array<double,7> stiffenedGas::consvToPrim(const std::array<double,7> arr)const{
    std::array<double,7> result;
    assert(arr[0] != 0);
    result[0] = arr[0];
    result[1] = arr[1]/arr[0];
    result[2] = (arr[2]-(arr[1]*arr[1]/(2.0*arr[0])))*(gamma-1)-gamma*p_inf;
    //print_vect(result);
    return result;
}

// constructors

EOS::EOS(double set_gamma, double set_Bx){
    gamma = set_gamma;
    BX = set_Bx;
}

idealGas::idealGas(double set_gamma, double set_Bx) : EOS(set_gamma, set_Bx){
}

stiffenedGas::stiffenedGas(double set_gamma, double set_p_inf, double set_e_inf, double set_Bx) : EOS(set_gamma, set_Bx){
    p_inf = set_p_inf;
    e_inf = set_e_inf;
}