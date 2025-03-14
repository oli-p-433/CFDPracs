#include "diffuseSolver.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>

void solver::run(){
    alpha = 0;
    std::cout << "pre-bc" << std::endl;

    setBCs(f);
    std::cout << "set BCs" << std::endl;

    splitFlip = 0;

    do{

        this->setDt();
        std::cout << "time = " << time << "dt is" << dt << std::endl;

        setBCs(f);

        int nSubcycle = 1;
        double sub = static_cast<double>(nSubcycle);

        splitFlip++;
        //double halfDt = dt/2.0;
        double halfDt = dt;
        if (splitFlip%2 == 0){
            
            f.uPrev = f.u;
            this->fluxMethod(XDIR); this->RK2(XDIR);
            this->fluxMethod(YDIR); this->RK2(YDIR);

            for (int i = 0; i < nSubcycle; ++i){
                if (splitFlip%4 == 0){
                    this->allaireSource(XDIR,halfDt/sub); this->allaireSource(YDIR,halfDt/sub);
                } else {
                    this->allaireSource(YDIR,halfDt/sub); this->allaireSource(XDIR,halfDt/sub);
                }
            }

        } else {

            f.uPrev = f.u;
            this->fluxMethod(YDIR); this->RK2(YDIR);
            this->fluxMethod(XDIR); this->RK2(XDIR);

            for (int i = 0; i < nSubcycle; ++i){
                if (splitFlip%4 == 0){
                    this->allaireSource(YDIR,halfDt/sub); this->allaireSource(XDIR,halfDt/sub);
                } else {
                    this->allaireSource(XDIR,halfDt/sub); this->allaireSource(YDIR,halfDt/sub);
                }
            }

            
        }


        


        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData();
        }
        
    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData();
};


void solver::SLIC(bool direction){
    if (PRIM == 1){
        //std::cout << "starting primitive SLIC" << std::endl;
        for (size_t i = 0; i < f.uPrim.size(); ++i){
            for (size_t j = 0; j < f.uPrim[0].size(); ++j){
                f.uPrim[i][j] = eos[0]->consvToPrim(f.u[i][j]);
            }
        }

    } else {
        //std::cout << "starting conservative SLIC" << std::endl;
    }

    //std::cout << "uPrim calculated, calc halfslopes" << std::endl;
    calcHalfSlopes(PRIM,direction);

    //std::cout << "calcr" << std::endl;
    calcr(direction);

    //std::cout << "calcUBars" << std::endl;
    calcUBars(PRIM,direction);

    //std::cout << "updUBars" << std::endl;
    updUBars(PRIM,direction);

    if (direction == XDIR){
        for (size_t i = 0; i < f.fluxesX.size(); ++i) {
            for (size_t j = 0; j < f.fluxesX[0].size(); ++j){
                f.fluxesX[i][j] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(f.uBarRX[i][j],f.uBarLX[i][j+1],eos[0],XDIR)),eos[0],XDIR)+LF(f.uBarRX[i][j],f.uBarLX[i][j+1],eos[0],XDIR)); // fluxesX stored in conserved variable form
            }
        }
    } else {
        for (size_t j = 0; j < f.fluxesY[0].size(); ++j) {
            for (size_t i = 0; i < f.fluxesY.size(); ++i){
                f.fluxesY[i][j] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(f.uBarRY[i][j],f.uBarLY[i+1][j],eos[0],YDIR)),eos[0],YDIR)+LF(f.uBarRY[i][j],f.uBarLY[i+1][j],eos[0],YDIR)); // fluxesX stored in conserved variable form
            }
        }
    }

    //std::cout << std::endl << "-----------------------fluxesX above-------------------------" << std::endl;
}

void solver::MUSCL(bool direction){
    
    if (PRIM == 1){
        //std::cout << "starting Primitive MUSCL" << std::endl;
        for (size_t i = 0; i < f.uPrim.size(); ++i){
            for (size_t j = 0; j < f.uPrim[0].size(); ++j){
                f.uPrim[i][j] = eos[0]->consvToPrim(f.u[i][j]);
            }
        }
    } else {
        //std::cout << "starting conservative MUSCL" << std::endl;
    }

    calcHalfSlopes(PRIM,direction);

    calcr(direction);

    calcUBars(PRIM,direction);

    updUBars(PRIM,direction);

    if (direction == XDIR){
        for (size_t i=0; i < f.fluxesX.size(); ++i){
            for (size_t j = 0; j < f.fluxesX[0].size(); ++j){
                f.fluxesX[i][j] = HLLC(f.uBarRX[i][j],f.uBarLX[i][j+1],eos[0],i,j,XDIR);
            }
        }
    } else {
        for (size_t j=0; j < f.fluxesY[0].size(); ++j){
            for (size_t i = 0; i < f.fluxesY.size(); ++i){
                f.fluxesY[i][j] = HLLC(f.uBarRY[i][j],f.uBarLY[i+1][j],eos[0],i,j,YDIR);
            }
        }
    }
}

void solver::HLLCGodunov(bool direction){
    if (direction == XDIR){
        for (size_t i=0; i < f.fluxesX.size(); ++i){
            for (size_t j = 0; j < f.fluxesX[0].size(); ++j){
                f.fluxesX[i][j] = HLLC(f.u[i+nG][j+nG-1],f.u[i+nG][j+nG],eos[0],i,j,XDIR);
            }
        }
    } else {
        for (size_t j=0; j < f.fluxesY[0].size(); ++j){
            for (size_t i = 0; i < f.fluxesY.size(); ++i){
                f.fluxesY[i][j] = HLLC(f.u[i+nG-1][j+nG],f.u[i+nG][j+nG],eos[0],i,j,YDIR);
            }
        }
    }
}

std::array<double,6> solver::HLLC(std::array<double,6> left,std::array<double,6> right, EOS* eos, int i, int j, bool direction){ // takes conserved variables
    double SL = 0.0, SR = 0.0;
    std::array<double,6> LPrim = eos->consvToPrim(left);
    std::array<double,6> RPrim = eos->consvToPrim(right);

    //printState(left); printState(right);
    //printState(LPrim); printState(RPrim);

    // Pressure - based wavespeed estimation
    /*
    double pEst = 0.5*(LPrim[PRES]+RPrim[PRES])-0.5*(RPrim[UX]-LPrim[UX])*0.25*(LPrim[RHO]+RPrim[RHO])*(eos->calcSoundSpeed(LPrim)+eos->calcSoundSpeed(RPrim));
    double qL,qR;
    
    qL = (pEst < LPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/LPrim[PRES])-1));
    qR = (pEst < RPrim[PRES]) ? 1 : sqrt(1+((eos->get_gamma()+1)/(2.0*eos->get_gamma()))*((pEst/RPrim[PRES])-1));

    SL = LPrim[UX] - eos->calcSoundSpeed(LPrim)*qL;
    SR = RPrim[UX] + eos->calcSoundSpeed(RPrim)*qR;
    */

    // easy wavespeed
    SL = (direction == XDIR) ? std::min(LPrim[UX] - eos->calcSoundSpeed(LPrim), RPrim[UX] - eos->calcSoundSpeed(RPrim)) : std::min(LPrim[UY] - eos->calcSoundSpeed(LPrim), RPrim[UY] - eos->calcSoundSpeed(RPrim));
    SR = (direction == XDIR) ? std::max(LPrim[UX] + eos->calcSoundSpeed(LPrim), RPrim[UX] + eos->calcSoundSpeed(RPrim)) : std::max(LPrim[UY] + eos->calcSoundSpeed(LPrim), RPrim[UY] + eos->calcSoundSpeed(RPrim));
    //std::cout << SL << " " << SR << std::endl;
    if ((std::isnan(SL) == 1) || (std::isnan(SR) == 1)){
        std::cout << "SL SR " << SL << " " << SR << std::endl;
    }

    double rhoL = left[RHO1]+left[RHO2];
    double rhoR = right[RHO1]+right[RHO2];

    // Calculate the numerator of sStar
    double numerator = (direction == XDIR) ? RPrim[PRES] - LPrim[PRES] + rhoL * LPrim[UX] * (SL - LPrim[UX]) - rhoR * RPrim[UX] * (SR - RPrim[UX])
                                        : RPrim[PRES] - LPrim[PRES] + rhoL * LPrim[UY] * (SL - LPrim[UY]) - rhoR * RPrim[UY] * (SR - RPrim[UY]);
    
    // Calculate the denominator of sStar
    double denominator = (direction == XDIR) ? rhoL * (SL - LPrim[UX]) - rhoR * (SR - RPrim[UX])
                                            : rhoL * (SL - LPrim[UY]) - rhoR * (SR - RPrim[UY]);

    
    // Calculate sStar
    double sStar = 0;
    if (numerator == 0 && denominator == 0){
        sStar = 0;
    } else if (denominator == 0){
        sStar = 1e15;
    } else {
        sStar = numerator / denominator;
    }
    
    
    if (std::isnan(sStar) == 1){
        std::cout << "numerator " <<  numerator << std::endl;
        std::cout << "denominator " << denominator << std::endl; 
        throw std::runtime_error("sStar is nan");
    }

    if (isGodunov == true){
        if (direction == XDIR){
            f.sStarsX[i][j] = sStar;
        } else {
            f.sStarsY[i][j] = sStar;
        }
    }
    
    //
    double prefL = (direction == XDIR) ? (SL-LPrim[UX])/(SL-sStar) : (SL-LPrim[UY])/(SL-sStar);
    double prefR = (direction == XDIR) ? (SR-RPrim[UX])/(SR-sStar) : (SR-RPrim[UY])/(SR-sStar);

    double eL,eR;
    eL = (direction == XDIR) ? (left[ENE]/rhoL)+(sStar-LPrim[UX])*(sStar+LPrim[PRES]/(rhoL*(SL-LPrim[UX])))
                            : (left[ENE]/rhoL)+(sStar-LPrim[UY])*(sStar+LPrim[PRES]/(rhoL*(SL-LPrim[UY])));
    eR = (direction == XDIR) ? (right[ENE]/rhoR)+(sStar-RPrim[UX])*(sStar+RPrim[PRES]/(rhoR*(SR-RPrim[UX])))
                            : (right[ENE]/rhoR)+(sStar-RPrim[UY])*(sStar+RPrim[PRES]/(rhoR*(SR-RPrim[UY])));


    std::array<double,6> uHLLCL,uHLLCR;
    uHLLCL = (direction == XDIR) ? std::array<double,6>{left[FRAC], left[RHO1]*prefL, left[RHO2]*prefL, rhoL*prefL*sStar, rhoL*prefL*LPrim[UY], rhoL*prefL*eL}
                                : std::array<double,6>{left[FRAC], left[RHO1]*prefL, left[RHO2]*prefL, rhoL*prefL*LPrim[UX], rhoL*prefL*sStar, rhoL*prefL*eL}; // ?
    
    uHLLCR = (direction == XDIR) ? std::array<double,6>{right[FRAC], right[RHO1]*prefR, right[RHO2]*prefR, rhoR*prefR*sStar, rhoR*prefR*RPrim[UY], rhoR*prefR*eR}
                                : std::array<double,6>{right[FRAC], right[RHO1]*prefR, right[RHO2]*prefR, rhoR*prefR*RPrim[UX], rhoR*prefR*sStar, rhoR*prefR*eR}; // ?

    std::array<double,6> alpFlux = {0,0,0,0,0,0};
    alpFlux[0] = (direction == XDIR) ? 0.5*(RPrim[UX]*right[FRAC]+LPrim[UX]*left[FRAC]) - 0.5*std::abs(0.5*(RPrim[UX]+LPrim[UX]))*(right[FRAC]-left[FRAC])
                                : 0.5*(RPrim[UY]*right[FRAC]+LPrim[UY]*left[FRAC]) - 0.5*std::abs(0.5*(RPrim[UY]+LPrim[UY]))*(right[FRAC]-left[FRAC]);
    
    if (0 <= SL){
        return flux(LPrim,eos,direction) + alpFlux;
    } else if (SL < 0 && 0 <= sStar){
        return flux(LPrim,eos,direction)+SL*(uHLLCL-left) + alpFlux;
    } else if (sStar < 0 && 0 <= SR){
        return flux(RPrim,eos,direction)+SR*(uHLLCR-right) + alpFlux;
    } else if (SR < 0){
        return flux(RPrim,eos,direction) + alpFlux;
    } else {
        throw std::runtime_error("wavespeed condition invalid");
    }

}

/*
void solver::setBCs(){
        for (int i = 0; i < nG; ++i){ // u BCs
            uPlus1[i] = u[nG];
            u[i] = u[nG];
            uPlus1[nCellsX+2*nG-1-i] = u[nCellsX+nG-1];
            u[nCellsX+2*nG-1-i] = u[nCellsX+nG-1];
        }
}
*/








void solver::pointsUpdate(bool direction){
    if (direction == XDIR){
        for (size_t i=nG; i < f.u.size()-nG; ++i){
            for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
                //std::cout << (dt / dx) * (fluxesX[i-nG+1] - fluxesX[i-nG])[0] << std::endl;
                f.uPlus1[i][j] = f.u[i][j] - (dt / dx) * (f.fluxesX[i-nG][j-nG+1] - f.fluxesX[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
            }
        }
    } else {
        for (size_t j=nG; j < f.u[0].size()-nG; ++j){
            for (size_t i = nG; i < f.u.size() - nG; ++i) {
                //std::cout << (dt / dx) * (fluxesX[i-nG+1] - fluxesX[i-nG])[0] << std::endl;
                f.uPlus1[i][j] = f.u[i][j] - (dt / dy) * (f.fluxesY[i-nG+1][j-nG] - f.fluxesY[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
            }
        }
    }
    f.u = f.uPlus1;
    setBCs(f);
    // print_arr(u,1);
    //std::cout << "updated u" << std::endl;
    // }
};

void solver::RK2(bool direction){
    std::vector< std::vector< std::array<double,6> > > uStored = f.u;
    if (direction == XDIR){
        this->pointsUpdate(XDIR); // gives updated u
        f.uPlus1 = f.u; // stores new u as uPlus1
        this->fluxMethod(XDIR); // calculates fluxes for new u

        for (size_t i = nG; i < f.u.size() - nG; ++i){
            for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
                f.uPlus1[i][j] = 0.5*(uStored[i][j]+f.uPlus1[i][j]) - 0.5 * (dt / dx) * (f.fluxesX[i-nG][j-nG+1] - f.fluxesX[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
            }
        }
    } else {
        this->pointsUpdate(YDIR); // gives updated u
        f.uPlus1 = f.u; // stores new u as uPlus1
        this->fluxMethod(YDIR); // calculates fluxes for new u

        for (size_t j = nG; j < f.u[0].size() - nG; ++j){
            for (size_t i = nG; i < f.u.size() - nG; ++i) {
                f.uPlus1[i][j] = 0.5*(uStored[i][j]+f.uPlus1[i][j]) - 0.5 * (dt / dx) * (f.fluxesY[i-nG+1][j-nG] - f.fluxesY[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
            }
        }

    }
    f.u = f.uPlus1;
    setBCs(f);

}

// ---------------- Source Terms ---------------------------- //

/*
void solver::sourceUpdate(){
    setBCs();
    std::vector< std::array<double,6> > sourceArr = sourceTerm(u,alpha);
    //std::cout << "printing source array" << std::endl;
    // print_arr(sourceArr,0);

    for (std::vector<double>::size_type i = nG; i < u.size() - nG; ++i) {
        uPlus1[i] = u[i] + dt * sourceArr[i-nG]; // flux[i + 1] and flux[i] for the update

        //print_vect(uPlus1[i]);
    }
    
    //std::cout << std::endl << "-------------------------source terms above------------------------" << std::endl;
    // for (int i = 0; i < u.size(); ++i) {
    u = uPlus1;
    std::cout << "evolved source terms" << std::endl;
    // }
};
*/

std::vector< std::array<double,6> > solver::sourceTerm(const std::vector< std::array<double,6> > arr, int alpha){

    for (size_t i = nG; i<(arr.size()-nG); ++i){
        double pref = -(alpha/((i-nG+0.5)*dx));
        sourceResult[i-nG][0] = pref*arr[i][1]; // rho*v
        sourceResult[i-nG][1] = pref*(arr[i][1]*arr[i][1])/arr[i][0]; // (rho*v)^2 / rho
        sourceResult[i-nG][2] = pref*(eos[0]->consvToPrim(arr[i])[1])*(arr[i][2]+eos[0]->consvToPrim(arr[i])[2]);
    }
    return sourceResult;
}

void solver::allaireSource(bool direction, double timestep){
    if (direction == XDIR){
        for (size_t i = nG; i < f.u.size() - nG; ++i) {
            for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
                //double sDiff = (u[i][XMOM]<0) ? (sStarsX[i-nG+1]-sStarsX[i-nG]) : (sStarsX[i-nG+2]-sStarsX[i-nG+1]);
                f.uPlus1[i][j][FRAC] = f.u[i][j][FRAC] + (timestep/dx) * f.uPrev[i][j][FRAC] * (f.sStarsX[i-nG][j-nG+1]-f.sStarsX[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
                //std::cout << u[i][FRAC] << " " << uPrev[i][FRAC] << " " << sStarsX[i-nG+1] << " " << sStarsX[i-nG] << " " << (sStarsX[i-nG+1]-sStarsX[i-nG]) << " " << (dt/dx) * uPrev[i][FRAC] * (sStarsX[i-nG+1]-sStarsX[i-nG]) << std::endl;

            }
        }
    } else {
        for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
            for (size_t i = nG; i < f.u.size() - nG; ++i) {
                //double sDiff = (u[i][XMOM]<0) ? (sStarsX[i-nG+1]-sStarsX[i-nG]) : (sStarsX[i-nG+2]-sStarsX[i-nG+1]);
                f.uPlus1[i][j][FRAC] = f.u[i][j][FRAC] + (timestep/dy) * f.uPrev[i][j][FRAC] * (f.sStarsY[i-nG+1][j-nG]-f.sStarsY[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
                //std::cout << u[i][FRAC] << " " << uPrev[i][FRAC] << " " << sStarsX[i-nG+1] << " " << sStarsX[i-nG] << " " << (sStarsX[i-nG+1]-sStarsX[i-nG]) << " " << (dt/dx) * uPrev[i][FRAC] * (sStarsX[i-nG+1]-sStarsX[i-nG]) << std::endl;

            }
        }
    }
    f.u = f.uPlus1;
    setBCs(f);
}

// -------------------- Flux functions --------------------- //

std::array<double,6> solver::finvectLF(const std::array<double,6> x)const{
    double a=1;
    return a*x;
};

std::array<double,6> solver::fBurgersLF(const std::array<double,6> x)const{
    return 0.5*(x * x);
};

std::array<double,6> solver::fEuler(std::array<double,6> arr, EOS* eos, bool direction){ // takes primitive variables
    std::array<double,6> result;
    if (direction == XDIR){
        result[FRAC] = 0; // vol fraction 
        result[RHO1] = arr[RHO1]*arr[UX]; // rho1*v
        result[RHO2] = arr[RHO2]*arr[UX]; // rho2*v
        result[XMOM] = (arr[RHO1]+arr[RHO2])*arr[UX]*arr[UX] + arr[PRES]; // rho*v^2 + p
        result[YMOM] = (arr[RHO1]+arr[RHO2])*arr[UX]*arr[UY]; // rho*vx*vy
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UX]; // (E + p)v
        return result;
    } else {
        result[FRAC] = 0; // vol fraction 
        result[RHO1] = arr[RHO1]*arr[UY]; // rho1*v
        result[RHO2] = arr[RHO2]*arr[UY]; // rho2*v
        result[XMOM] = (arr[RHO1]+arr[RHO2])*arr[UX]*arr[UY]; // rho*vx*vy
        result[YMOM] = (arr[RHO1]+arr[RHO2])*arr[UY]*arr[UY] + arr[PRES]; // rho*v^2 + p
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UY]; // (E + p)v
        return result;
    }

}

// ------------------------- secondary fluxesX --------------------- //

std::array<double,6> solver::LF(std::array<double,6> v1, std::array<double,6> v2, EOS* eos,bool direction)const{
    if (direction == XDIR){
        return 0.5 * (flux(eos->consvToPrim(v1),eos,XDIR) + flux(eos->consvToPrim(v2),eos,XDIR)) + 0.5 * (dx / dt) * (v1 - v2);
    } else {
        return 0.5 * (flux(eos->consvToPrim(v1),eos,YDIR) + flux(eos->consvToPrim(v2),eos,YDIR)) + 0.5 * (dy / dt) * (v1 - v2);
    }
    
};

std::array<double,6> solver::laxFhalf(std::array<double,6> ui, std::array<double,6> uip1, EOS* eos, bool direction)const{
    if (direction == XDIR){
        return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos->consvToPrim(uip1),eos,XDIR)-flux(eos->consvToPrim(ui),eos,XDIR)));
    } else {
        return 0.5*(ui+uip1)-(0.5*(dt/dy)*(flux(eos->consvToPrim(uip1),eos,YDIR)-flux(eos->consvToPrim(ui),eos,YDIR)));
    }
};



// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(bool prim, bool direction){
    if (direction == XDIR){
        for (size_t i=0; i < f.halfSlopesX.size(); ++i){
            for (size_t j = 0; j < f.halfSlopesX[0].size(); ++j){ // halfslopes[i] = halfslopes i-1/2
                f.halfSlopesX[i][j] = (prim == 0) ? f.u[i+nG][j+1]-f.u[i+nG][j] : f.uPrim[i+nG][j+1]-f.uPrim[i+nG][j];
            }
        }
    } else {
        for (size_t j=0; j < f.halfSlopesY[0].size(); ++j){
            for (size_t i = 0; i < f.halfSlopesY.size(); ++i){ // halfslopes[i] = halfslopes i-1/2
                f.halfSlopesY[i][j] = (prim == 0) ? f.u[i+1][j+nG]-f.u[i][j+nG] : f.uPrim[i+1][j+nG]-f.uPrim[i][j+nG];
            }
        }
    }

};

void solver::calcr(bool direction){
    if (direction == XDIR){
        for (size_t i=0; i < f.rX.size(); ++i){
            for (size_t j = 0; j < f.rX[0].size(); ++j){
                f.rX[i][j] = elementDivide(f.halfSlopesX[i][j],f.halfSlopesX[i][j+1]);
            }
        }
    } else {
        for (size_t j=0; j < f.rY[0].size(); ++j){
            for (size_t i = 0; i < f.rY.size(); ++i){
                f.rY[i][j] = elementDivide(f.halfSlopesY[i][j],f.halfSlopesY[i+1][j]);
            }
        }
    }
    //print_arr(r);
    //std::cout << "calcd r" << std::endl;
};


/*

+---------------+--------------+--------------+
|   u[0][0]     |   u[0][1]    |    u[0][2]   |  <- Top row (Row 0)
|               |              |  (u[0][nG])  |
+---------------+--------------+--------------+
|               |              |              |  <- Middle row (Row 1)
|               |              |              |
+---------------+--------------+--------------+
|               |  uBar[0][0]  |   u[nG][nG]  |  <- Bottom row (Row 2)
|        halfSlopesX[0][0]  f[0][0](&sStarsX)   |
|               |   r[0][0]    |              |
+---------------+--------------+--------------+

*/


void solver::calcUBars(bool prim, bool direction){ // calculates for all u except leftmost cell
    if (direction == XDIR){
        for (size_t i = 0; i < f.uBarLX.size(); ++i){
            for (size_t j = 0; j < f.uBarLX[0].size(); ++j){
                f.uBarLX[i][j] = (prim == 0) ? f.u[i+nG][j+nG-1] - 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1])) : f.uPrim[i+nG][j+nG-1] - 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1]));
                f.uBarRX[i][j] = (prim == 0) ? f.u[i+nG][j+nG-1] + 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1])) : f.uPrim[i+nG][j+nG-1] + 0.5 * slopeLim(f.rX[i][j]) * (0.5 * (f.halfSlopesX[i][j]+f.halfSlopesX[i][j+1]));
            }
        }
        for (size_t i = 0; i < f.sStarsX.size(); ++i){
            for (size_t j = 0; j < f.sStarsX[0].size(); ++j){
                f.sStarsX[i][j] = (prim == 0) ? 0.5*(eos[0]->consvToPrim(f.uBarRX[i][j])[UX]+eos[0]->consvToPrim(f.uBarLX[i][j+1])[UX]) : 0.5*(f.uBarRX[i][j][UX]+f.uBarLX[i][j+1][UX]);
            }
        }
    } else {
        for (size_t j = 0; j < f.uBarLY[0].size(); ++j){
            for (size_t i = 0; i < f.uBarLY.size(); ++i){
                f.uBarLY[i][j] = (prim == 0) ? f.u[i+nG-1][j+nG] - 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j])) : f.uPrim[i+nG-1][j+nG] - 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j]));
                f.uBarRY[i][j] = (prim == 0) ? f.u[i+nG-1][j+nG] + 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j])) : f.uPrim[i+nG-1][j+nG] + 0.5 * slopeLim(f.rY[i][j]) * (0.5 * (f.halfSlopesY[i][j]+f.halfSlopesY[i+1][j]));
            }
        }
        for (size_t j = 0; j < f.sStarsY[0].size(); ++j){
            for (size_t i = 0; i < f.sStarsY.size(); ++i){
                f.sStarsY[i][j] = (prim == 0) ? 0.5*(eos[0]->consvToPrim(f.uBarRY[i][j])[UY]+eos[0]->consvToPrim(f.uBarLY[i+1][j])[UY]) : 0.5*(f.uBarRY[i][j][UY]+f.uBarLY[i+1][j][UY]);
            }
        }
    }

};

std::array<double,6> multiplyMatrixVector(
    const std::vector<std::vector<double>>& B,
    const std::array<double,6>& x) {
    
    std::array<double,6> result = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            result[i] += B[i][j] * x[j];
        }
    }
    
    return result;
}

void solver::updUBars(bool prim,bool direction){
    if (direction == XDIR){
        if (prim == false){
            for (size_t i = 0; i<f.uBarLupdX.size(); ++i){
                for (size_t j = 0; j<f.uBarLupdX[0].size(); ++j){
                    f.uBarLupdX[i][j] = f.uBarLX[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(f.uBarRX[i][j]),eos[0],XDIR)-flux(eos[0]->consvToPrim(f.uBarLX[i][j]),eos[0],XDIR));
                    f.uBarRupdX[i][j] = f.uBarRX[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(f.uBarRX[i][j]),eos[0],XDIR)-flux(eos[0]->consvToPrim(f.uBarLX[i][j]),eos[0],XDIR));
                }
            }
        } else {
            std::vector< std::vector <double> > B;
            for (size_t i = 0; i<f.uBarLupdX.size(); ++i){
                for (size_t j = 0; j<f.uBarLupdX[0].size(); ++j){
                    double v = f.uPrim[i+nG][j+nG-1][UX]; double rho = f.uPrim[i+nG][j+nG-1][RHO1]+f.uPrim[i+nG][j+nG-1][RHO2];
                    B = {   {v, 0,  0,  0,                           0,                              0},
                            {0, v,  0,  f.uPrim[i+nG][j+nG-1][RHO1], 0,                              0},
                            {0, 0,  v,  f.uPrim[i+nG][j+nG-1][RHO2], 0,                              0},
                            {0, 0,  0,  v,                           0,                        1.0/rho},
                            {0, 0,  0,  0,                           v,                              0},
                            {0, 0,  0,  rho*pow(eos[0]->calcSoundSpeed(f.uPrim[i+nG][j+nG-1]),2), 0, v}     };

                    std::array<double,6> deltaU = f.uBarRX[i][j]-f.uBarLX[i][j];
                    f.uBarLupdX[i][j] = eos[0]->primToConsv(f.uBarLX[i][j]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
                    f.uBarRupdX[i][j] = eos[0]->primToConsv(f.uBarRX[i][j]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
                }
            }

        }
        f.uBarLX = f.uBarLupdX;
        f.uBarRX = f.uBarRupdX;
    } else {
        if (prim == false){
            for (size_t j = 0; j<f.uBarLupdY[0].size(); ++j){
                for (size_t i = 0; i<f.uBarLupdY.size(); ++i){
                    f.uBarLupdY[i][j] = f.uBarLY[i][j]-0.5*(dt/dy)*(flux(eos[0]->consvToPrim(f.uBarRY[i][j]),eos[0],YDIR)-flux(eos[0]->consvToPrim(f.uBarLY[i][j]),eos[0],YDIR));
                    f.uBarRupdY[i][j] = f.uBarRY[i][j]-0.5*(dt/dy)*(flux(eos[0]->consvToPrim(f.uBarRY[i][j]),eos[0],YDIR)-flux(eos[0]->consvToPrim(f.uBarLY[i][j]),eos[0],YDIR));
                }
            }
        } else {
            std::vector< std::vector <double> > B;
            for (size_t j = 0; j<f.uBarLupdY[0].size(); ++j){
                for (size_t i = 0; i<f.uBarLupdY.size(); ++i){
                    double v = f.uPrim[i+nG-1][j+nG][UY]; double rho = f.uPrim[i+nG-1][j+nG][RHO1]+f.uPrim[i+nG-1][j+nG][RHO2];
                    B = {   {v, 0,  0,  0,                           0,                              0},
                            {0, v,  0,  0, f.uPrim[i+nG-1][j+nG][RHO1],                              0},
                            {0, 0,  v,  0, f.uPrim[i+nG-1][j+nG][RHO2],                              0},
                            {0, 0,  0,  v,                           0,                              0},
                            {0, 0,  0,  0,                           v,                        1.0/rho},
                            {0, 0,  0,  0, rho*pow(eos[0]->calcSoundSpeed(f.uPrim[i+nG-1][j+nG]),2), v}     };

                    std::array<double,6> deltaU = f.uBarRY[i][j]-f.uBarLY[i][j];
                    f.uBarLupdY[i][j] = eos[0]->primToConsv(f.uBarLY[i][j]-0.5*(dt/dy)*multiplyMatrixVector(B,deltaU));
                    f.uBarRupdY[i][j] = eos[0]->primToConsv(f.uBarRY[i][j]-0.5*(dt/dy)*multiplyMatrixVector(B,deltaU));
                }
            }

        }
        f.uBarLY = f.uBarLupdY;
        f.uBarRY = f.uBarRupdY;
    }
};

std::array<double,6> solver::calcSlope(double omega, std::array<double,6> slopeleft, std::array<double,6> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,6> solver::minbee(std::array<double,6> slp){
    std::array<double,6> minbArr{0,0,0,0,0,0};
    double minb;
    for (int k = 0; k<6; ++k){
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
    //return {minBval,minBval,minBval,minBval};
    return minbArr;
};

std::array<double,6> solver::superbee(std::array<double,6> slp){
    std::array<double,6> minbArr{0,0,0,0,0,0};
    double minb;
    for (int k = 0; k<6; ++k){
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

std::array<double,6> solver::vanLeer(std::array<double,6> slp){
    std::array<double,6> minbArr{0,0,0,0,0,0};
    for (int k = 0; k<6; ++k){
        if (slp[k] <= 0 || slp[k] == INFINITY){
            minbArr[k] = 0;
        } else {
            minbArr[k] = std::min(2.0*slp[k]/(1.0+slp[k]),2.0/(1+slp[k]));
        }
    };
    auto minminB = std::min_element(minbArr.begin(),minbArr.end());
    double minBval = *minminB;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    return minbArr;
    //return {minBval,minBval,minBval,minBval};
};


// ---------------------- Misc ------------------------------ //

void solver::print_arr(std::vector< std::array<double,6> > arr, int var){
    for (std::vector<double>::size_type i = 0; i < arr.size(); ++i) {
        std::cout << arr[i][var] << " ";
    }
    std::cout << std::endl;
}

void solver::print_vect(std::array<double,6> v){
    for (int j = 0; j<4; ++j){
        std::cout << v[j] << " ";
        //std::cout << "laxFhalf " << laxFhalf(u[i],u[i+1])[j] << std::endl;
    }
    std::cout << std::endl;
}

void solver::setDt(){
    double aMax = -1e14;
    std::cout << "setting dt" << std::endl;
    for (size_t i = nG; i < f.u.size()-nG; ++i) {
        for (size_t j = nG; j < f.u[0].size()-nG; ++j) {
            //std::cout << "sound speed " << eos[0]->calcSoundSpeed(u[i]) << std::endl;
            aMax = std::max(aMax, std::sqrt(pow(eos[0]->consvToPrim(f.u[i][j])[UX],2)+pow(eos[0]->consvToPrim(f.u[i][j])[UY],2)) + eos[0]->calcSoundSpeed(eos[0]->consvToPrim(f.u[i][j])));
        }
    }

    dtReducer = (splitFlip < 20) ? 0.5 : 1.0;

    dt = dtReducer * cour * std::min(dx,dy) / std::fabs(aMax);

    std::cout << "set dt" << std::endl;

    std::cout << "aMax =" << aMax << std::endl;


    double writeTime = writeInterval*(std::floor((time+1e-9) / writeInterval)+1);

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

    if (std::fabs(writeTime - time) < 1e-9){
        checkWrite = 1;
    } else {
        checkWrite = 0;
    }
    // std::cout << "dt is" << dt << std::endl;

};

// ---------------------- Constructor ------------------------ //
solver::solver(double x_0, double x_1, double y_0, double y_1, double t0, double t1, int n_CellsX, int n_CellsY, int n_Ghosts, double c)
    : f(fluid(n_CellsX, n_CellsY, n_Ghosts)), bound(boundary(n_CellsX,n_CellsY,n_Ghosts)),
        nCellsX(n_CellsX), nCellsY(n_CellsY), nG(n_Ghosts), 
        x0(x_0), x1(x_1),
        y0(y_0), y1(y_1),
        startTime(t0), endTime(t1),
        
        cour(c){
    
    assert(t1>t0);
    timeMulti = 1000; // 1.0/(endTime-startTime);

    dx = (x1 - x0)/nCellsX;
    dy = (y1 - y0)/nCellsY;
    

    // fluid 2

    sourceResult.resize(nCellsX);

    slopeWeight = 1;

    cour = c;

    std::vector<std::string> variables = {"frac","rho1","rho2","v","p"};

    //varMap["rho"] = 0; varMap["v"] = 1; varMap["p"] = 2;
    for (size_t i = 0; i<variables.size(); ++i){
        varMap[variables[i]] = i;
    }
    
};

// ---------------------- utilities -------------------------- //

int solver::ghosts(){
    return nG;
};

std::array<double,2> solver::get_dxdy()const{
    std::array<double,2> result = {dx,dy};
    return result;
};


void solver::init(std::vector < std::vector< std::array<double,6> > > init) {
    assert(init.size() == f.u.size());
    assert(init[0].size() == f.u[0].size());
    f.u = init;
};

// -------------------- writing functions --------------------- //

void solver::setWriteInterval(double wt){
    writeInterval = wt;
}

void solver::writeData() const{
        // making data file
        std::string filename = dirName + "/" + std::to_string(time*timeMulti).substr(0,5);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        writeFile << "x y alpha rho1 rho2 vx vy p" << std::endl;
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (size_t i=0; i < f.u.size();i++){
                double y = y0 + (i-static_cast<double>(nG)+0.5)*dy;
                for (size_t j=0; j < f.u[0].size();j++){
                    double x = x0 + (j-static_cast<double>(nG)+0.5)*dx;
                    std::array<double,6> prim = eos[0]->consvToPrim(f.u[i][j]);
                    writeFile << x << " " << y << " " << prim[FRAC] << " " << prim[RHO1] << " " << prim[RHO2] << " " << prim[UX] << " " << prim[UY] << " " << prim[PRES] << std::endl;
                    //std::cout << u[i] << " ";
                }
                writeFile << "\n";
            }
            writeFile.close();
            std::cout << filename << "written successfully" << std::endl;
        }else{
            std::cout << "failed to write." << std::endl;
        }

}



