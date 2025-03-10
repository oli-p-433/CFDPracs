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

    do{

        this->setDt();
        // dt = 0.1*dx;
        std::cout << "dt is" << dt << std::endl;

        setBCs(f);

        //this->SLIC();
        //this->MUSCL();
        f.uPrev = f.u;
        this->fluxMethod();
        //std::cout << "Fluxes" << std::endl;
        //print_arr(fluxes,0);
        //std::cout << "u" << std::endl;

        //print_arr(u,0);
        this->pointsUpdate();
        //std::cout << "points updated" << std::endl;
        //print_arr(u,0);
        this->allaireSource();
        //std::cout << "allaire source updated" << std::endl;
        //print_arr(u,0);

        


        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData();
        }
        
    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData();
};


void solver::SLIC(){
    if (PRIM == 1){
        std::cout << "starting primitive SLIC" << std::endl;
        for (size_t i = 0; i < f.uPrim.size(); ++i){
            for (size_t j = 0; j < f.uPrim[0].size(); ++j){
                f.uPrim[i][j] = eos[0]->consvToPrim(f.u[i][j]);
            }
        }

    } else {
        std::cout << "starting conservative SLIC" << std::endl;
    }

    std::cout << "uPrim calculated, calc halfslopes" << std::endl;
    calcHalfSlopes(PRIM);

    std::cout << "calcr" << std::endl;
    calcr();

    std::cout << "calcUBars" << std::endl;
    calcUBars(PRIM);

    std::cout << "updUBars" << std::endl;
    updUBars(PRIM);

    for (size_t i = 0; i < f.fluxes.size(); ++i) {
        for (size_t j = 0; j < f.fluxes[0].size(); ++j){
            f.fluxes[i][j] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(f.uBarR[i][j],f.uBarL[i][j+1],eos[0])),eos[0])+LF(f.uBarR[i][j],f.uBarL[i][j+1],eos[0])); // fluxes stored in conserved variable form
        }
    }

    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
}

void solver::MUSCL(){
    
    if (PRIM == 1){
        std::cout << "starting Primitive MUSCL" << std::endl;
        for (size_t i = 0; i < f.uPrim.size(); ++i){
            for (size_t j = 0; j < f.uPrim[0].size(); ++j){
                f.uPrim[i][j] = eos[0]->consvToPrim(f.u[i][j]);
            }
        }
    } else {
        std::cout << "starting conservative MUSCL" << std::endl;
    }

    calcHalfSlopes(PRIM);

    calcr();

    calcUBars(PRIM);

    updUBars(PRIM);

    for (size_t i=0; i < f.fluxes.size(); ++i){
        for (size_t j = 0; j < f.fluxes[0].size(); ++j){
            f.fluxes[i][j] = HLLC(f.uBarR[i][j],f.uBarL[i][j+1],eos[0],i,j);
        }
    }
}

void solver::HLLCGodunov(){

    for (size_t i=0; i < f.fluxes.size(); ++i){
        for (size_t j = 0; j < f.fluxes[0].size(); ++j){
            f.fluxes[i][j] = HLLC(f.u[i+nG][j+nG-1],f.u[i+nG][j+nG],eos[0],i,j);
        }
    }
}

std::array<double,5> solver::HLLC(std::array<double,5> left,std::array<double,5> right, EOS* eos, int i, int j){ // takes conserved variables
    double SL = 0.0, SR = 0.0;
    std::array<double,5> LPrim = eos->consvToPrim(left);
    std::array<double,5> RPrim = eos->consvToPrim(right);

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
    SL = std::min(LPrim[UX] - eos->calcSoundSpeed(LPrim), RPrim[UX] - eos->calcSoundSpeed(RPrim));
    SR = std::max(LPrim[UX] + eos->calcSoundSpeed(LPrim), RPrim[UX] + eos->calcSoundSpeed(RPrim));
    //std::cout << SL << " " << SR << std::endl;
    if ((std::isnan(SL) == 1) || (std::isnan(SR) == 1)){
        std::cout << "SL SR " << SL << " " << SR << std::endl;
    }

    double rhoL = left[RHO1]+left[RHO2];
    double rhoR = right[RHO1]+right[RHO2];

    // Calculate the numerator of sStar
    double numerator = RPrim[PRES] - LPrim[PRES] + rhoL * LPrim[UX] * (SL - LPrim[UX]) - rhoR * RPrim[UX] * (SR - RPrim[UX]);
    
    // Calculate the denominator of sStar
    double denominator = rhoL * (SL - LPrim[UX]) - rhoR * (SR - RPrim[UX]);

    
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
        std::cout << "numerator " <<  (RPrim[PRES]-LPrim[PRES]+rhoL*LPrim[UX]*(SL-LPrim[UX])-rhoR*RPrim[UX]*(SR-RPrim[UX])) << std::endl;
        std::cout << "denominator " << (rhoL*(SL-LPrim[UX]) - rhoR*(SR-RPrim[UX])) << std::endl; 
        throw std::runtime_error("sStar is nan");
    }

    f.sStars[i][j] = sStar;
    //
    double prefL = (SL-LPrim[UX])/(SL-sStar);
    double prefR = (SR-RPrim[UX])/(SR-sStar);

    double eL,eR;
    eL = (left[ENE]/rhoL)+(sStar-LPrim[UX])*(sStar+LPrim[PRES]/(rhoL*(SL-LPrim[UX])));
    eR = (right[ENE]/rhoR)+(sStar-RPrim[UX])*(sStar+RPrim[PRES]/(rhoR*(SR-RPrim[UX])));


    std::array<double,5> uHLLCL,uHLLCR;
    uHLLCL = {left[FRAC], left[RHO1]*prefL, left[RHO2]*prefL, rhoL*prefL*sStar, rhoL*prefL*eL}; // ?
    uHLLCR = {right[FRAC], right[RHO1]*prefR, right[RHO2]*prefR, rhoR*prefR*sStar, rhoR*prefR*eR}; // ?

    
    if (0 <= SL){
        return flux(LPrim,eos);
    } else if (SL < 0 && 0 <= sStar){
        return flux(LPrim,eos)+SL*(uHLLCL-left);
    } else if (sStar < 0 && 0 <= SR){
        return flux(RPrim,eos)+SR*(uHLLCR-right);
    } else if (SR < 0){
        return flux(RPrim,eos);
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








void solver::pointsUpdate(){
    setBCs(f);
    for (size_t i=nG; i < f.u.size()-nG; ++i){
        for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
            //std::cout << (dt / dx) * (fluxes[i-nG+1] - fluxes[i-nG])[0] << std::endl;
            f.uPlus1[i][j] = f.u[i][j] - (dt / dx) * (f.fluxes[i-nG][j-nG+1] - f.fluxes[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
        }
    }
    f.u = f.uPlus1;
    // print_arr(u,1);
    //std::cout << "updated u" << std::endl;
    // }
};

// ---------------- Source Terms ---------------------------- //

/*
void solver::sourceUpdate(){
    setBCs();
    std::vector< std::array<double,5> > sourceArr = sourceTerm(u,alpha);
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

std::vector< std::array<double,5> > solver::sourceTerm(const std::vector< std::array<double,5> > arr, int alpha){

    for (size_t i = nG; i<(arr.size()-nG); ++i){
        double pref = -(alpha/((i-nG+0.5)*dx));
        sourceResult[i-nG][0] = pref*arr[i][1]; // rho*v
        sourceResult[i-nG][1] = pref*(arr[i][1]*arr[i][1])/arr[i][0]; // (rho*v)^2 / rho
        sourceResult[i-nG][2] = pref*(eos[0]->consvToPrim(arr[i])[1])*(arr[i][2]+eos[0]->consvToPrim(arr[i])[2]);
    }
    return sourceResult;
}

void solver::allaireSource(){
    setBCs(f);
    for (size_t i = nG; i < f.u.size() - nG; ++i) {
        for (size_t j = nG; j < f.u[0].size() - nG; ++j) {
            //double sDiff = (u[i][XMOM]<0) ? (sStars[i-nG+1]-sStars[i-nG]) : (sStars[i-nG+2]-sStars[i-nG+1]);
            f.uPlus1[i][j][FRAC] = f.u[i][j][FRAC] + (dt/dx) * f.uPrev[i][j][FRAC] * (f.sStars[i-nG][j-nG+1]-f.sStars[i-nG][j-nG]); // flux[i + 1] and flux[i] for the update
            //std::cout << u[i][FRAC] << " " << uPrev[i][FRAC] << " " << sStars[i-nG+1] << " " << sStars[i-nG] << " " << (sStars[i-nG+1]-sStars[i-nG]) << " " << (dt/dx) * uPrev[i][FRAC] * (sStars[i-nG+1]-sStars[i-nG]) << std::endl;

        }
    }
    
    f.u = f.uPlus1;
}

// -------------------- Flux functions --------------------- //

std::array<double,5> solver::finvectLF(const std::array<double,5> x)const{
    double a=1;
    return a*x;
};

std::array<double,5> solver::fBurgersLF(const std::array<double,5> x)const{
    return 0.5*(x * x);
};

std::array<double,5> solver::fEuler(std::array<double,5> arr, EOS* eos){ // takes primitive variables
    std::array<double,5> result;
    result[FRAC] = arr[FRAC]*arr[UX]; // vol fraction 
    result[RHO1] = arr[RHO1]*arr[UX]; // rho1*v
    result[RHO2] = arr[RHO2]*arr[UX]; // rho2*v
    result[XMOM] = (arr[RHO1]+arr[RHO2])*arr[UX]*arr[UX] + arr[PRES]; // rho*v^2 + p
    result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UX]; // (E + p)v
    return result;
}

// ------------------------- secondary fluxes --------------------- //

std::array<double,5> solver::LF(std::array<double,5> v1, std::array<double,5> v2, EOS* eos)const{
    return 0.5 * (flux(eos->consvToPrim(v1),eos) + flux(eos->consvToPrim(v2),eos)) + 0.5 * (dx / dt) * (v1 - v2);
};

std::array<double,5> solver::laxFhalf(std::array<double,5> ui, std::array<double,5> uip1, EOS* eos)const{
    return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos->consvToPrim(uip1),eos)-flux(eos->consvToPrim(ui),eos)));
};



// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(bool prim){
    for (size_t i=0; i < f.halfSlopes.size(); ++i){
        for (size_t j = 0; j < f.halfSlopes[0].size(); ++j){ // halfslopes[i] = halfslopes i-1/2
            f.halfSlopes[i][j] = (prim == 0) ? f.u[i+nG][j+1]-f.u[i+nG][j] : f.uPrim[i+nG][j+1]-f.uPrim[i+nG][j];
        }
    }
};

void solver::calcr(){
    for (size_t i=0; i < f.r.size(); ++i){
        for (size_t j = 0; j < f.r[0].size(); ++j){
            f.r[i][j] = elementDivide(f.halfSlopes[i][j],f.halfSlopes[i][j+1]);
        }
    }
    //print_arr(r);
    std::cout << "calcd r" << std::endl;
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
|        halfSlopes[0][0]  f[0][0](&sStars)   |
|               |   r[0][0]    |              |
+---------------+--------------+--------------+

*/


void solver::calcUBars(bool prim){ // calculates for all u except leftmost cell
    for (size_t i = 0; i < f.uBarL.size(); ++i){
        for (size_t j = 0; j < f.uBarL[0].size(); ++j){
            f.uBarL[i][j] = (prim == 0) ? f.u[i+nG][j+nG-1] - 0.5 * slopeLim(f.r[i][j]) * (0.5 * (f.halfSlopes[i][j]+f.halfSlopes[i][j+1])) : f.uPrim[i+nG][j+nG-1] - 0.5 * slopeLim(f.r[i][j]) * (0.5 * (f.halfSlopes[i][j]+f.halfSlopes[i][j+1]));
            f.uBarR[i][j] = (prim == 0) ? f.u[i+nG][j+nG-1] + 0.5 * slopeLim(f.r[i][j]) * (0.5 * (f.halfSlopes[i][j]+f.halfSlopes[i][j+1])) : f.uPrim[i+nG][j+nG-1] + 0.5 * slopeLim(f.r[i][j]) * (0.5 * (f.halfSlopes[i][j]+f.halfSlopes[i][j+1]));
        }
    }
    for (size_t i = 0; i < f.sStars.size(); ++i){
        for (size_t j = 0; j < f.sStars[0].size(); ++j){
            f.sStars[i][j] = (prim == 0) ? 0.5*(eos[0]->consvToPrim(f.uBarR[i][j])[UX]+eos[0]->consvToPrim(f.uBarL[i][j+1])[UX]) : 0.5*(f.uBarR[i][j][UX]+f.uBarL[i][j+1][UX]);
        }
    }

};

std::array<double, 5> multiplyMatrixVector(
    const std::vector<std::vector<double>>& B,
    const std::array<double, 5>& x) {
    
    std::array<double, 5> result = {0.0, 0.0, 0.0, 0.0, 0.0};
    
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            result[i] += B[i][j] * x[j];
        }
    }
    
    return result;
}

void solver::updUBars(bool prim){
    if (prim == 0){
        for (size_t i = 0; i<f.uBarLupd.size(); ++i){
            for (size_t j = 0; j<f.uBarLupd[0].size(); ++j){
                f.uBarLupd[i][j] = f.uBarL[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(f.uBarR[i][j]),eos[0])-flux(eos[0]->consvToPrim(f.uBarL[i][j]),eos[0]));
                f.uBarRupd[i][j] = f.uBarR[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(f.uBarR[i][j]),eos[0])-flux(eos[0]->consvToPrim(f.uBarL[i][j]),eos[0]));
            }
        }
    } else {
        std::vector< std::vector <double> > B;
        for (size_t i = 0; i<f.uBarLupd.size(); ++i){
            for (size_t j = 0; j<f.uBarLupd[0].size(); ++j){
                double v = f.uPrim[i+nG][j+nG-1][UX]; double rho = f.uPrim[i+nG][j+nG-1][RHO1]+f.uPrim[i+nG][j+nG-1][RHO2];
                B = {   {v, 0,  0,  0,                                            0},
                        {0, v,  0,  f.uPrim[i+nG][j+nG-1][RHO1],                               0},
                        {0, 0,  v,  f.uPrim[i+nG][j+nG-1][RHO2],                               0},
                        {0, 0,  0,  v,                                      1.0/rho},
                        {0, 0,  0,  rho*pow(eos[0]->calcSoundSpeed(f.uPrim[i+nG][j+nG-1]),2)  ,v}  };

                std::array<double,5> deltaU = f.uBarR[i][j]-f.uBarL[i][j];
                f.uBarLupd[i][j] = eos[0]->primToConsv(f.uBarL[i][j]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
                f.uBarRupd[i][j] = eos[0]->primToConsv(f.uBarR[i][j]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
            }
        }

    }
    f.uBarL = f.uBarLupd;
    f.uBarR = f.uBarRupd;
};

std::array<double,5> solver::calcSlope(double omega, std::array<double,5> slopeleft, std::array<double,5> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,5> solver::minbee(std::array<double,5> slp){
    std::array<double,5> minbArr{0,0,0};
    double minb;
    for (int k = 0; k<4; ++k){
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
    return {minBval,minBval,minBval,minBval};
    //return minbArr;
};

std::array<double,5> solver::superbee(std::array<double,5> slp){
    std::array<double,5> minbArr{0,0,0};
    double minb;
    for (int k = 0; k<4; ++k){
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



std::array<double,5> solver::vanLeer(std::array<double,5> slp){
    std::array<double,5> minbArr{0,0,0};
    for (int k = 0; k<4; ++k){
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
    return {minBval,minBval,minBval,minBval};
};


// ---------------------- Misc ------------------------------ //

void solver::print_arr(std::vector< std::array<double,5> > arr, int var){
    for (std::vector<double>::size_type i = 0; i < arr.size(); ++i) {
        std::cout << arr[i][var] << " ";
    }
    std::cout << std::endl;
}

void solver::print_vect(std::array<double,5> v){
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
            aMax = std::max(aMax, std::abs(eos[0]->consvToPrim(f.u[i][j])[UX]) + eos[0]->calcSoundSpeed(f.u[i][j]));
        }
    }
    dt = cour * dx / std::fabs(aMax);

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

    if (std::fabs(writeTime - time) < 1e-6){
        checkWrite = 1;
    } else {
        checkWrite = 0;
    }
    // std::cout << "dt is" << dt << std::endl;

};

// ---------------------- Constructor ------------------------ //
solver::solver(double x_0, double x_1, double y_0, double y_1, double t0, double t1, int n_CellsX, int n_CellsY, int n_Ghosts, double c)
    : f(fluid(n_CellsX, n_CellsY, n_Ghosts)), bound(boundary(n_CellsX,n_CellsY,n_Ghosts)),
        x0(x_0), x1(x_1),
        y0(y_0), y1(y_1),
        startTime(t0), endTime(t1),
        nCellsX(n_CellsX), nCellsY(n_CellsY), nG(n_Ghosts), 
        cour(c){
    
    assert(t1>t0);
    timeMulti = 1; // 1.0/(endTime-startTime);

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


void solver::init(std::vector < std::vector< std::array<double,5> > > init) {
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
        writeFile << "x y alpha rho1 rho2 v p" << std::endl;
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (size_t i=0; i < f.u.size();i++){
                double y = y0 + (i-static_cast<double>(nG)+0.5)*dy;
                for (size_t j=0; j < f.u[0].size();j++){
                    double x = x0 + (j-static_cast<double>(nG)+0.5)*dx;
                    std::array<double,5> prim = eos[0]->consvToPrim(f.u[i][j]);
                    writeFile << x << " " << y << " " << prim[0] << " " << prim[1] << " " << prim[2] << " " << prim[3] << " " << prim[4] << std::endl;
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



