#include "solver1D.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>

void solver::run(){
    alpha = 0;
    std::cout << "pre-bc" << std::endl;

    transmissiveBC();
    std::cout << "set BCs" << std::endl;

    do{
        this->findBoundary();
        this->reinitPhi();
        this->interfaceRiem();
        this->setGhostFluids();


        this->setDt();
        // dt = 0.1*dx;
        std::cout << "dt is" << dt << std::endl;

        transmissiveBC();


        this->phiUpdate();

        //this->SLIC();
        //this->MUSCL();
        this->HLLCGodunov();
        this->pointsUpdate();
        std::cout << "points updated" << std::endl;

        //this->sourceUpdate();
        this->interface_cells={};
        this->interfaceStates={};

        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData("rho"); writeData("v"); writeData("p");
        }
        
    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData("rho"); writeData("v"); writeData("p");
};


void solver::SLIC(){
    std::cout << "starting SLIC" << std::endl;

    calcHalfSlopes();

    calcr();

    calcUBars();

    updUBars();

    for (size_t i = 0; i < fluxes1.size(); ++i) {
        fluxes1[i] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(u1BarRupd[i],u1BarLupd[i+1],eos[0])),eos[0])+LF(u1BarRupd[i],u1BarLupd[i+1],eos[0])); // fluxes stored in conserved variable form
        fluxes2[i] = 0.5*(flux(eos[1]->consvToPrim(laxFhalf(u2BarRupd[i],u2BarLupd[i+1],eos[1])),eos[1])+LF(u2BarRupd[i],u2BarLupd[i+1],eos[1])); // fluxes stored in conserved variable form

    }

    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
}

void solver::MUSCL(){
    std::cout << "starting MUSCL" << std::endl;

    calcHalfSlopes();

    calcr();

    calcUBars();

    updUBars();

    for (size_t i=0; i < u1BarRupd.size()-1; ++i){
        riemann solution1(eos[0]->get_gamma(),eos[0]->get_gamma(),eos[0]->consvToPrim(u1BarRupd[i]),eos[0]->consvToPrim(u1BarLupd[i+1]),0,1,0,1,2,0,0); // need p_inf setup
        fluxes1[i] = flux(solution1.exctRiemann(),eos[0]);

        riemann solution2(eos[1]->get_gamma(),eos[1]->get_gamma(),eos[1]->consvToPrim(u2BarRupd[i]),eos[1]->consvToPrim(u2BarLupd[i+1]),0,1,0,1,2,0,0); // need p_inf setup
        fluxes2[i] = flux(solution2.exctRiemann(),eos[1]);
    }
}

void solver::exactGodunov(){

    for (size_t i=0; i < fluxes1.size(); ++i){
        riemann solution1(eos[0]->get_gamma(),eos[0]->get_gamma(),eos[0]->consvToPrim(u1[i+nGhost-1]),eos[0]->consvToPrim(u1[i+nGhost]),0,1,0,1,2,0,0); // need p_inf setup
        fluxes1[i] = flux(solution1.exctRiemann(),eos[0]);

        riemann solution2(eos[1]->get_gamma(),eos[1]->get_gamma(),eos[1]->consvToPrim(u2[i+nGhost-1]),eos[1]->consvToPrim(u2[i+nGhost]),0,1,0,1,2,0,0); // need p_inf setup
        fluxes2[i] = flux(solution2.exctRiemann(),eos[1]);
    }
}

void solver::HLLCGodunov(){

    for (size_t i=0; i < fluxes1.size(); ++i){
        riemann solution1(eos[0]->get_gamma(),eos[0]->get_gamma(),eos[0]->consvToPrim(u1[i+nGhost-1]),eos[0]->consvToPrim(u1[i+nGhost]),0,1,0,1,2,0,0); // need p_inf setup
        fluxes1[i] = flux(solution1.exctRiemann(),eos[0]);

        riemann solution2(eos[1]->get_gamma(),eos[1]->get_gamma(),eos[1]->consvToPrim(u2[i+nGhost-1]),eos[1]->consvToPrim(u2[i+nGhost]),0,1,0,1,2,0,0); // need p_inf setup
        fluxes2[i] = flux(solution2.exctRiemann(),eos[1]);
    }
}

std::array<double,3> solver::HLLC(std::array<double,3> left,std::array<double,3> right, EOS* eos){ // takes conserved variables
    std::array<double,3> HLLCFlux;
    double SL = 0.0, SR = 0.0;
    std::array<double,3> LPrim = eos->consvToPrim(left);
    std::array<double,3> RPrim = eos->consvToPrim(right);

    // Pressure - based wavespeed estimation
    /*
    double pEst = 0.5*(LPrim[PRES]+RPrim[PRES])-0.5*(RPrim[dir]-LPrim[dir])*0.25*(LPrim[RHO]+RPrim[RHO])*(eos->calcSoundSpeed(LPrim)+eos->calcSoundSpeed(RPrim));
    double qL,qR;
    
    qL = (pEst < LPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/LPrim[PRES])-1));
    qR = (pEst < RPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/RPrim[PRES])-1));

    SL = LPrim[dir] - eos->calcSoundSpeed(LPrim)*qL;
    SR = RPrim[dir] + eos->calcSoundSpeed(RPrim)*qR;
    */

    // easy wavespeed (doesnt work)
    SL = std::min(LPrim[1] - eos->calcSoundSpeed(LPrim), RPrim[1] - eos->calcSoundSpeed(RPrim));
    SR = std::max(LPrim[1] + eos->calcSoundSpeed(LPrim), RPrim[1] + eos->calcSoundSpeed(RPrim));

    if ((std::isnan(SL) == 1) || (std::isnan(SR) == 1)){
        std::cout << "SL SR " << SL << " " << SR << std::endl;
    }

    // Calculate the numerator of sStar
    double numerator = RPrim[2] - LPrim[2] + LPrim[0] * LPrim[1] * (SL - LPrim[1]) - RPrim[0] * RPrim[1] * (SR - RPrim[1]);
    
    // Calculate the denominator of sStar
    double denominator = LPrim[0] * (SL - LPrim[1]) - RPrim[0] * (SR - RPrim[1]);
    
    // Calculate sStar
    double sStar = numerator / denominator;
                    (LPrim[0]*(SL-LPrim[1]) - RPrim[0]*(SR-RPrim[1]));
    if (std::isnan(sStar) == 1){

        std::cout << "numerator " <<  (RPrim[2]-LPrim[2]+LPrim[0]*LPrim[1]*(SL-LPrim[1])-RPrim[0]*RPrim[1]*(SR-RPrim[1])) << std::endl;
        std::cout << "denominator " << (LPrim[0]*(SL-LPrim[1]) - RPrim[0]*(SR-RPrim[1])) << std::endl; 
        throw std::runtime_error("sStar is nan");
    }

    // calculating u_HLLC
    std::array<double,3> Dstar;
    Dstar = {0, 1, sStar};


    double SdiffL{SL - sStar},SdiffR{SR - sStar};
    if (std::abs(SdiffL) < 1e-8) SdiffL = (SL < 0) ? -1e-8 : 1e-8;
    if (std::abs(SdiffR) < 1e-8) SdiffR = (SR < 0) ? -1e-8 : 1e-8;

    std::array<double,3> FhllcLeft = (sStar*(SL*LPrim-flux(LPrim,eos))+SL*(LPrim[2]+LPrim[0]*(SL-LPrim[1])*(sStar-LPrim[1]))*Dstar)/(SdiffL);
    std::array<double,3> FhllcRight = (sStar*(SR*RPrim-flux(RPrim,eos))+SR*(RPrim[2]+RPrim[0]*(SR-RPrim[1])*(sStar-RPrim[1]))*Dstar)/(SdiffR);

    if (0 <= SL){
        return flux(LPrim,eos);
    } else if (SL < 0 && 0 <= sStar){
        return FhllcLeft;
    } else if (sStar < 0 && 0 <= SR){
        return FhllcRight;
    } else if (SR < 0){
        return flux(RPrim,eos);
    } else {
        throw std::runtime_error("wavespeed condition invalid");
    }
}

void solver::transmissiveBC(){
    if (alpha == 0){
        for (int i = 0; i < nGhost; ++i){ // u BCs
            u1Plus1[i] = u1[nGhost];
            u1[i] = u1[nGhost];
            u1Plus1[nCells+2*nGhost-1-i] = u1[nCells+nGhost-1];
            u1[nCells+2*nGhost-1-i] = u1[nCells+nGhost-1];
        }
        for (int i = 0; i < nGhost; ++i){ // u BCs
            u2Plus1[i] = u2[nGhost];
            u2[i] = u2[nGhost];
            u2Plus1[nCells+2*nGhost-1-i] = u2[nCells+nGhost-1];
            u2[nCells+2*nGhost-1-i] = u2[nCells+nGhost-1];
        }
        for (int i = 0; i < nGhost; ++i){ // phi BCs
            phiPlus1[i] = phi[nGhost];
            phi[i] = phi[nGhost];
            phiPlus1[nCells+2*nGhost-1-i] = phi[nCells+nGhost-1];
            phi[nCells+2*nGhost-1-i] = phi[nCells+nGhost-1];
        }
    } else {
        assert(0 != 0);
        /*
        for (int i = 0; i < nGhost; ++i){
            for (int var = 0; var<3; ++var){
                if (var == 0 || var == 2){
                    // Outside (transmissive)
                    uPlus1[nCells+2*nGhost-1-i][var] = u[nCells+nGhost-1][var];
                    u[nCells+2*nGhost-1-i][var] = u[nCells+nGhost-1][var];
                    // Central (reflective)
                    u[nGhost-1-i][var] = u[nGhost+i][var];
                    uPlus1[nGhost-1-i][var] = u[nGhost+i][var];
                } else {
                    // external (transmissive)
                    uPlus1[nCells+2*nGhost-1-i][var] = u[nCells+nGhost-1][var];
                    u[nCells+2*nGhost-1-i][var] = u[nCells+nGhost-1][var];
                    // central (reflective)
                    u[nGhost-1-i][var] = (-1)*u[nGhost+i][var];
                    uPlus1[nGhost-1-i][var] = (-1)*u[nGhost+i][var];
                }

            }
        } */
    }
}

void solver::findBoundary(){ // checks for change in sign of level set. emplaces index of cell left of the boundary.
    for (size_t i = 0; i < phi.size(); ++i) {  // Loop through all cells (except last)
        if (phi[i] * phi[i+1] < 0) {  // Sign change detected
            //std::cout << phi[i] << " " << phi[i+1] << std::endl;
            interface_cells.push_back(i);
            std::cout << "interface at cell " << i << std::endl; 
        }
    }
    interfaceStates.resize(interface_cells.size());
    
}

void solver::interfaceRiem(){
    for (size_t i=0; i<interface_cells.size(); ++i){
        print_vect(eos[0]->consvToPrim(u1[interface_cells[i]])); print_vect(eos[1]->consvToPrim(u2[interface_cells[i]+1]));
        riemann solution(eos[0]->get_gamma(),eos[1]->get_gamma(),eos[0]->consvToPrim(u1[interface_cells[i]]),eos[1]->consvToPrim(u2[interface_cells[i]+1]),0,1,0.5,0.12,100,0,0); // need p_inf setup
        interfaceStates[i] = solution.exctRiemann();
    }
}

void solver::setGhostFluids(){
    for (int i=interface_cells[0]+1;i<nCells+2*nGhost;++i){
        u1[i]=eos[0]->primToConsv(interfaceStates[0]);
        //std::cout << i << " " << phi[i] << std::endl;
        u2[nCells+2*nGhost-(i+1)]=eos[1]->primToConsv(interfaceStates[0]);
    }
}

void solver::reinitPhi(){
    for (size_t i=0; i<(interface_cells[0]); ++i){ // reinitialisation
        phi[i] = (static_cast<int>(i) - static_cast<int>(interface_cells[0]) - 0.5)*dx;
    }
    for (size_t i=interface_cells[0]+2; i<phi.size(); ++i){ // reinitialisation
        phi[i] = (static_cast<int>(i) - static_cast<int>(interface_cells[0]) - 0.5)*dx;
    }
}




void solver::phiUpdate(){
    transmissiveBC();
    for (std::vector<double>::size_type i = nGhost; i < phi.size()-nGhost; ++i) {
        double uVal{0};
        if (phi[i]<=0){
            uVal = eos[0]->consvToPrim(u1[i])[1];
        } else if (phi[i] > 0) {
            uVal = eos[1]->consvToPrim(u2[i])[1];
        };
        if (uVal >= 0){
            phiPlus1[i] = phi[i] - uVal*(dt/dx)*(phi[i]-phi[i-1]);
        } else {
            phiPlus1[i] = phi[i] - uVal*(dt/dx)*(phi[i+1]-phi[i]);
        }
        
    }

    phi = phiPlus1;

}

void solver::pointsUpdate(){
    transmissiveBC();

    for (std::vector<double>::size_type i = nGhost; i < u1.size() - nGhost; ++i) {
        u1Plus1[i] = u1[i] - (dt / dx) * (fluxes1[i-nGhost+1] - fluxes1[i-nGhost]); // flux[i + 1] and flux[i] for the update
    }

    for (std::vector<double>::size_type i = nGhost; i < u2.size() - nGhost; ++i) {
        u2Plus1[i] = u2[i] - (dt / dx) * (fluxes2[i-nGhost+1] - fluxes2[i-nGhost]); // flux[i + 1] and flux[i] for the update
    }

    u1 = u1Plus1;
    u2 = u2Plus1;
    // print_arr(u,1);
    std::cout << "updated u" << std::endl;
    // }
};

// ---------------- Source Terms ---------------------------- //

void solver::sourceUpdate(){
    transmissiveBC();
    std::vector< std::array<double,3> > sourceArr = sourceTerm(u1,alpha);
    std::cout << "printing source array" << std::endl;
    // print_arr(sourceArr,0);

    for (std::vector<double>::size_type i = nGhost; i < u1.size() - nGhost; ++i) {
        u1Plus1[i] = u1[i] + dt * sourceArr[i-nGhost]; // flux[i + 1] and flux[i] for the update

        //print_vect(uPlus1[i]);
    }
    
    std::cout << std::endl << "-------------------------source terms above------------------------" << std::endl;
    // for (int i = 0; i < u.size(); ++i) {
    u1 = u1Plus1;
    std::cout << "evolved source terms" << std::endl;
    // }
};

std::vector< std::array<double,3> > solver::sourceTerm(const std::vector< std::array<double,3> > arr, int alpha){

    for (size_t i = nGhost; i<(arr.size()-nGhost); ++i){
        double pref = -(alpha/((i-nGhost+0.5)*dx));
        sourceResult[i-nGhost][0] = pref*arr[i][1]; // rho*v
        sourceResult[i-nGhost][1] = pref*(arr[i][1]*arr[i][1])/arr[i][0]; // (rho*v)^2 / rho
        sourceResult[i-nGhost][2] = pref*(eos[0]->consvToPrim(arr[i])[1])*(arr[i][2]+eos[0]->consvToPrim(arr[i])[2]);
    }
    return sourceResult;
}

// -------------------- Flux functions --------------------- //

std::array<double,3> solver::finvectLF(const std::array<double,3> x)const{
    double a=1;
    return a*x;
};

std::array<double,3> solver::fBurgersLF(const std::array<double,3> x)const{
    return 0.5*(x * x);
};

std::array<double,3> solver::fEuler(std::array<double,3> arr, EOS* eos){ // takes primitive variables
    std::array<double,3> result;
    result[0] = arr[0]*arr[1]; // rho*v
    result[1] = arr[0]*arr[1]*arr[1]+arr[2]; // rho*v^2 + p
    result[2] = (eos->primToConsv(arr)[2] + arr[2])*arr[1]; // (E + p)v
    return result;
}

// ------------------------- secondary fluxes --------------------- //

std::array<double,3> solver::LF(std::array<double,3> v1, std::array<double,3> v2, EOS* eos)const{
    return 0.5 * (flux(eos->consvToPrim(v1),eos) + flux(eos->consvToPrim(v2),eos)) + 0.5 * (dx / dt) * (v1 - v2);
};

std::array<double,3> solver::laxFhalf(std::array<double,3> ui, std::array<double,3> uip1, EOS* eos)const{
    return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos->consvToPrim(uip1),eos)-flux(eos->consvToPrim(ui),eos)));
};



// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(){
    for (size_t i = 0; i < halfSlopes1.size(); ++i){
        halfSlopes1[i] = u1[i+nGhost-1]-u1[i+nGhost-2];
    }
    for (size_t i = 0; i < halfSlopes2.size(); ++i){
        halfSlopes2[i] = u2[i+nGhost-1]-u2[i+nGhost-2];
    }
};

void solver::calcr(){
    for (size_t i = 0; i < r1.size(); ++i){
        r1[i] = elementDivide(halfSlopes1[i],halfSlopes1[i+1]);
    }
    for (size_t i = 0; i < r2.size(); ++i){
        r2[i] = elementDivide(halfSlopes2[i],halfSlopes2[i+1]);
    }
    //print_arr(r);
    std::cout << "calcd r" << std::endl;
};


void solver::calcUBars(){ // calculates for all u except leftmost cell
    for (size_t i = 0; i<u1BarL.size(); ++i){
        u1BarL[i] = u1[i+nGhost-1] - 0.5 * ( slopeLim(r1[i]) * calcSlope(slopeWeight, halfSlopes1[i], halfSlopes1[i+1]) );
        u1BarR[i] = u1[i+nGhost-1] + 0.5 * ( slopeLim(r1[i]) * calcSlope(slopeWeight, halfSlopes1[i], halfSlopes1[i+1]) );
    }
    for (size_t i = 0; i<u2BarL.size(); ++i){
        u2BarL[i] = u2[i+nGhost-1] - 0.5 * ( slopeLim(r2[i]) * calcSlope(slopeWeight, halfSlopes2[i], halfSlopes2[i+1]) );
        u2BarR[i] = u2[i+nGhost-1] + 0.5 * ( slopeLim(r2[i]) * calcSlope(slopeWeight, halfSlopes2[i], halfSlopes2[i+1]) );
    }
};


void solver::updUBars(){
    for (size_t i = 0; i<u1BarLupd.size(); ++i){
        u1BarLupd[i] = u1BarL[i]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(u1BarR[i]),eos[0])-flux(eos[0]->consvToPrim(u1BarL[i]),eos[0]));
        u1BarRupd[i] = u1BarR[i]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(u1BarR[i]),eos[0])-flux(eos[0]->consvToPrim(u1BarL[i]),eos[0]));
    }
    for (size_t i = 0; i<u2BarLupd.size(); ++i){
        u2BarLupd[i] = u2BarL[i]-0.5*(dt/dx)*(flux(eos[1]->consvToPrim(u2BarR[i]),eos[1])-flux(eos[1]->consvToPrim(u2BarL[i]),eos[1]));
        u2BarRupd[i] = u2BarR[i]-0.5*(dt/dx)*(flux(eos[1]->consvToPrim(u2BarR[i]),eos[1])-flux(eos[1]->consvToPrim(u2BarL[i]),eos[1]));
    }
};

std::array<double,3> solver::calcSlope(double omega, std::array<double,3> slopeleft, std::array<double,3> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,3> solver::minbee(std::array<double,3> slp){
    std::array<double,3> minbArr{0,0,0};
    double minb;
    for (int k = 0; k<3; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY){
            minb = 0;
        } else if (slp[k] > 1){
            minb = std::min(1.0,2.0/(1+slp[k]));
        } else {
            minb = slp[k];
        }
        minbArr[k] = minb;
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    return {minBval,minBval,minBval};
    //return minbArr;
};

std::array<double,3> solver::superbee(std::array<double,3> slp){
    std::array<double,3> minbArr{0,0,0};
    double minb;
    for (int k = 0; k<3; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY){
            minb = 0;
        } else if (slp[k] <= 0.5){
            minb = 2*slp[k];
        } else if (slp[k] <= 1.0){
            minb = 1.0;
        } else {
            minb = std::min(slp[k],std::min(2.0/(1+slp[k]),2.0));
        }
        minbArr[k] = minb;
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    //return {minBval,minBval,minBval};
    return minbArr;
};



std::array<double,3> solver::vanLeer(std::array<double,3> slp){
    std::array<double,3> minbArr{0,0,0};
    for (int k = 0; k<3; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY){
            minbArr[k] = 0;
        } else {
            minbArr[k] = std::min(2.0*slp[k]/(1.0+slp[k]),2.0/(1+slp[k]));
        }
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    //return minbArr;
    return {minBval,minBval,minBval};
};


// ---------------------- Misc ------------------------------ //

void solver::print_arr(std::vector< std::array<double,3> > arr, int var){
    for (std::vector<double>::size_type i = 0; i < arr.size(); ++i) {
        std::cout << arr[i][var] << " ";
    }
    std::cout << std::endl;
}

void solver::print_vect(std::array<double,3> v){
    for (int j = 0; j<3; ++j){
        std::cout << v[j] << " ";
        //std::cout << "laxFhalf " << laxFhalf(u[i],u[i+1])[j] << std::endl;
    }
    std::cout << std::endl;
}

void solver::setDt(){
    double aMax = -1e14;
    std::cout << "setting dt" << std::endl;
    for (std::vector<double>::size_type i = nGhost; i < u1.size()-nGhost; ++i) {
        aMax = std::max(aMax, std::max(std::fabs(eos[1]->consvToPrim(u2[i])[1])+eos[1]->calcSoundSpeed(eos[1]->consvToPrim(u2[i])), std::fabs(eos[0]->consvToPrim(u1[i])[1])+eos[0]->calcSoundSpeed(eos[0]->consvToPrim(u1[i]))));
    }
    dt = cour * dx / std::fabs(aMax);

    std::cout << "set dt" << std::endl;

    std::cout << "aMax =" << aMax << std::endl;


    double writeTime = writeInterval*(std::floor((time+1e-12) / writeInterval)+1);

    if (endTime - time <= dt && writeTime-time >= dt){
        dt = endTime - time;
        assert(dt != 0);
        time = endTime;
        std::cout << "last step:" << dt << std::endl;
    } else if (writeTime-time <= dt && writeTime-time > 0){ // manually reducing dt to hit the write times
        dt = writeTime - time;
        assert(dt != 0);
        time = writeTime;
    } else {
        time = time + dt;
    }

    if (std::fabs(writeTime - time) < 1e-6){
        checkWrite = 1;
    } else {
        checkWrite = 0;
    }
    // std::cout << "dt is" << dt << std::endl;

};

// ---------------------- Constructor ------------------------ //
solver::solver(double x_0, double x_1, double t0, double t1, int n_Cells, int n_Ghosts, double c){
    assert(t1>t0);
    x0 = x_0; x1 = x_1;
    startTime = t0; endTime = t1;
    timeMulti = 1.0/(endTime-startTime);
    nCells = n_Cells;
    nGhost = n_Ghosts;
    dx = (x1 - x0)/nCells;

    phi.resize(nCells+2*nGhost);
    phiPlus1.resize(nCells+2*nGhost);

    u1Plus1.resize(nCells+2*nGhost);

    u1.resize(nCells+2*nGhost);
    fluxes1.resize(nCells+1);
    halfSlopes1.resize(nCells+3);
    r1.resize(nCells+2); // r can only be defined for inner cells
    u1BarL.resize(nCells+2);
    u1BarR.resize(nCells+2);
    u1BarLupd.resize(nCells+2);
    u1BarRupd.resize(nCells+2);
    // fluid 2
    u2Plus1.resize(nCells+2*nGhost);
    u2.resize(nCells+2*nGhost);
    fluxes2.resize(nCells+1);
    halfSlopes2.resize(nCells+3);
    r2.resize(nCells+2); // r can only be defined for inner cells
    u2BarL.resize(nCells+2);
    u2BarR.resize(nCells+2);
    u2BarLupd.resize(nCells+2);
    u2BarRupd.resize(nCells+2);

    sourceResult.resize(nCells);

    slopeWeight = 1;

    cour = c;

    std::vector<std::string> variables = {"rho","v","p"};

    //varMap["rho"] = 0; varMap["v"] = 1; varMap["p"] = 2;
    for (size_t i = 0; i<variables.size(); ++i){
        varMap[variables[i]] = i;
    }
    
};

// ---------------------- utilities -------------------------- //

int solver::ghosts(){
    return nGhost;
};

double solver::get_dx()const{
    return dx;
};

void solver::init(std::vector< std::array<double,3> > init1,std::vector< std::array<double,3> > init2) {
    assert(init1.size() == u1.size() && init2.size() == u2.size());
    u1 = init1;
    u2 = init2;
};

void solver::phiInit(std::vector< double > init) {
    assert(init.size() == phi.size());
    phi = init;
};

// -------------------- writing functions --------------------- //

void solver::setWriteInterval(double wt){
    writeInterval = wt;
}

void solver::writeData(std::string varName) const{
        // making data file
        std::string filename = dirName + "/" + varName + "/" + std::to_string(time*timeMulti).substr(0,4);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (std::vector<double>::size_type i=nGhost; i<(u1.size()-nGhost);i++){
                double x = x0 + (i-nGhost+0.5)*dx;
                writeFile << x << " " << eos[0]->consvToPrim(u1[i])[varMap.at(varName)] << " " << eos[1]->consvToPrim(u2[i])[varMap.at(varName)] << " " << phi[i] << std::endl;
                //std::cout << u[i] << " ";
            }
            writeFile.close();
            std::cout << filename << "written successfully" << std::endl;
        }else{
            std::cout << "failed to write." << std::endl;
        }

}

// ------------------- Overloading operators ------------------ //

/*
std::array<double, 3> operator-(const std::array<double, 3>& lhs, const std::array<double, 3>& rhs) {
    std::array<double, 3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

std::array<double, 3> operator+(const std::array<double, 3>& lhs, const std::array<double, 3>& rhs) {
    std::array<double, 3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

std::array<double, 3> operator*(const double scalar,const std::array<double, 3>& lhs) {
    std::array<double, 3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        result[i] = scalar*lhs[i];
    }
    return result;
}

std::array<double, 3> operator*(const std::array<double, 3>& rhs,const std::array<double, 3>& lhs) {
    std::array<double, 3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        result[i] = rhs[i]*lhs[i];
    }
    return result;
}

std::array<double, 3> elementDivide(const std::array<double, 3>& lhs, const std::array<double, 3>& rhs) {
    std::array<double, 3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        if (rhs[i] == lhs[i]){
            result[i] = 1;
        }
        else if (rhs[i] == 0){
            result[i] = INFINITY;
        }
        else{
            result[i] = lhs[i] / rhs[i];
        }
    }
    return result;
}

std::array<double, 3> operator/(const std::array<double, 3>& lhs, const double scalar) {
    std::array<double, 3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        result[i] = lhs[i]/scalar;
    }
    return result;
}

*/



