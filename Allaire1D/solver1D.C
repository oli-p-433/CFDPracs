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

    int nSubcycle = 1;
    double sub = static_cast<double>(nSubcycle);

    writeData();

    std::cout << "first update done, starting loop" << std::endl;

    do{

        this->setDt();

        // dt = 0.1*dx;
        std::cout << "dt is" << dt << std::endl;

        uPrev = u;
        this->fluxMethod();
        //this->pointsUpdate();
        this->RK2();

        for (int i = 0; i < nSubcycle; ++i){
            this->allaireSource(dt/sub);
        }
        uPrev = u;

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
        for (size_t i = 0; i < uPrim.size(); ++i){
            uPrim[i] = eos[0]->consvToPrim(u[i]);
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

    for (size_t i = 0; i < fluxes.size(); ++i) {
        fluxes[i] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(uBarR[i],uBarL[i+1],eos[0])),eos[0])+LF(uBarR[i],uBarL[i+1],eos[0])); // fluxes stored in conserved variable form
    }

    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
}

void solver::MUSCL(){
    
    if (PRIM == 1){
        std::cout << "starting Primitive MUSCL" << std::endl;
        for (size_t i = 0; i < uPrim.size(); ++i){
            uPrim[i] = eos[0]->consvToPrim(u[i]);
        }
    } else {
        std::cout << "starting conservative MUSCL" << std::endl;
    }

    calcHalfSlopes(PRIM);

    calcr();

    calcUBars(PRIM);

    updUBars(PRIM);

    for (size_t i=0; i < fluxes.size(); ++i){
        fluxes[i] = HLLC(uBarR[i],uBarL[i+1],eos[0],i);
    }
}


void solver::HLLCGodunov(){

    for (size_t i=0; i < fluxes.size(); ++i){
        fluxes[i] = HLLC(u[i+nGhost-1],u[i+nGhost],eos[0],i);
    }
}

std::array<double,5> solver::HLLC(std::array<double,5> left,std::array<double,5> right, EOS* eos, int i){ // takes conserved variables
    double SL = 0.0, SR = 0.0;
    std::array<double,5> LPrim = eos->consvToPrim(left);
    std::array<double,5> RPrim = eos->consvToPrim(right);

    double rhoL = left[RHO1]+left[RHO2]; double rhoR = right[RHO1]+right[RHO2];

    // Pressure - based wavespeed estimation
    /*
    double xiL = left[FRAC]*(1/(eos->get_gamma()[1]-1)) + (1-left[FRAC])*(1/(eos->get_gamma()[0]-1));
    double xiR = right[FRAC]*(1/(eos->get_gamma()[1]-1)) + (1-right[FRAC])*(1/(eos->get_gamma()[0]-1));

    double gammaL = 1 + 1/xiL; double gammaR = 1 + 1/xiR;
    
    double pEst = 0.5*(LPrim[PRES]+RPrim[PRES])-0.5*(RPrim[UX]-LPrim[UX])*0.25*(rhoL+rhoR)*(eos->calcSoundSpeed(LPrim)+eos->calcSoundSpeed(RPrim));
    double qL,qR;
    
    qL = (pEst < LPrim[PRES]) ? 1 : sqrt(1+((gammaL+1)/(2.0*gammaL))*((pEst/LPrim[PRES])-1));
    qR = (pEst < RPrim[PRES]) ? 1 : sqrt(1+((gammaR+1)/(2.0*gammaR))*((pEst/RPrim[PRES])-1));

    SL = LPrim[UX] - eos->calcSoundSpeed(LPrim)*qL;
    SR = RPrim[UX] + eos->calcSoundSpeed(RPrim)*qR;
    */
    

    // easy wavespeed
    SL = std::min(LPrim[UX] - eos->calcSoundSpeed(LPrim), RPrim[UX] - eos->calcSoundSpeed(RPrim));
    SR = std::max(LPrim[UX] + eos->calcSoundSpeed(LPrim), RPrim[UX] + eos->calcSoundSpeed(RPrim));


    if ((std::isnan(SL) == 1) || (std::isnan(SR) == 1)){
        std::cout << "SL SR " << SL << " " << SR << std::endl;
    }

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

    sStars[i] = sStar;
    //
    double prefL = (SL-LPrim[UX])/(SL-sStar);
    double prefR = (SR-RPrim[UX])/(SR-sStar);

    double eL,eR;
    eL = (left[ENE]/rhoL)+(sStar-LPrim[UX])*(sStar+LPrim[PRES]/(rhoL*(SL-LPrim[UX])));
    eR = (right[ENE]/rhoR)+(sStar-RPrim[UX])*(sStar+RPrim[PRES]/(rhoR*(SR-RPrim[UX])));


    std::array<double,5> uHLLCL,uHLLCR;
    uHLLCL = {left[FRAC], left[RHO1]*prefL, left[RHO2]*prefL, rhoL*prefL*sStar, rhoL*prefL*eL}; // ?
    uHLLCR = {right[FRAC], right[RHO1]*prefR, right[RHO2]*prefR, rhoR*prefR*sStar, rhoR*prefR*eR}; // ?

    std::array<double,5> alpFlux = {0,0,0,0,0};
    alpFlux[0] = 0.5*(RPrim[UX]*right[FRAC]+LPrim[UX]*left[FRAC]) - 0.5*std::abs(0.5*(RPrim[UX]+LPrim[UX]))*(right[FRAC]-left[FRAC]);
    //alpFlux[0] = 0.5*(RPrim[UX]*right[FRAC]+LPrim[UX]*left[FRAC]);

    if (0 <= SL){
        return flux(LPrim,eos) + alpFlux;
    } else if (SL < 0 && 0 <= sStar){
        return flux(LPrim,eos)+SL*(uHLLCL-left) + alpFlux;
    } else if (sStar < 0 && 0 <= SR){
        return flux(RPrim,eos)+SR*(uHLLCR-right) + alpFlux;
    } else if (SR < 0){
        return flux(RPrim,eos) + alpFlux;
    } else {
        throw std::runtime_error("wavespeed condition invalid");
    }

}

void solver::transmissiveBC(){
        for (int i = 0; i < nGhost; ++i){ // u BCs
            uPlus1[i] = u[nGhost];
            u[i] = u[nGhost];
            uPlus1[nCells+2*nGhost-1-i] = u[nCells+nGhost-1];
            u[nCells+2*nGhost-1-i] = u[nCells+nGhost-1];
        }
}








void solver::pointsUpdate(){

    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        // 1st order
        uPlus1[i] = u[i] - (dt / dx) * (fluxes[i-nGhost+1] - fluxes[i-nGhost]); // flux[i + 1] and flux[i] for the update
    }
    u = uPlus1;
    transmissiveBC();
    // print_arr(u,1);
    //std::cout << "updated u" << std::endl;
    // }
};

void solver::RK2(){
    std::vector<std::array<double,5>> uStored = u;
    this->pointsUpdate(); // gives updated u
    uPlus1 = u; // stores new u as uPlus1
    this->fluxMethod(); // calculates fluxes for new u

    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        uPlus1[i] = 0.5*(uStored[i]+uPlus1[i]) - 0.5 * (dt / dx) * (fluxes[i-nGhost+1] - fluxes[i-nGhost]); // flux[i + 1] and flux[i] for the update
    }
    u = uPlus1;
    transmissiveBC();

}


// ---------------- Source Terms ---------------------------- //

void solver::sourceUpdate(){
    transmissiveBC();
    std::vector< std::array<double,5> > sourceArr = sourceTerm(u,alpha);
    //std::cout << "printing source array" << std::endl;
    // print_arr(sourceArr,0);

    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        uPlus1[i] = u[i] + dt * sourceArr[i-nGhost]; // flux[i + 1] and flux[i] for the update

        //print_vect(uPlus1[i]);
    }
    
    //std::cout << std::endl << "-------------------------source terms above------------------------" << std::endl;
    // for (int i = 0; i < u.size(); ++i) {
    u = uPlus1;
    std::cout << "evolved source terms" << std::endl;
    // }
};

std::vector< std::array<double,5> > solver::sourceTerm(const std::vector< std::array<double,5> > arr, int alpha){

    for (size_t i = nGhost; i<(arr.size()-nGhost); ++i){
        double pref = -(alpha/((i-nGhost+0.5)*dx));
        sourceResult[i-nGhost][0] = pref*arr[i][1]; // rho*v
        sourceResult[i-nGhost][1] = pref*(arr[i][1]*arr[i][1])/arr[i][0]; // (rho*v)^2 / rho
        sourceResult[i-nGhost][2] = pref*(eos[0]->consvToPrim(arr[i])[1])*(arr[i][2]+eos[0]->consvToPrim(arr[i])[2]);
    }
    return sourceResult;
}

void solver::allaireSource(double timestep){

    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        
        // Murrone & guillard source term (compaction)
        
        std::array<double,5> uPrevPrim = eos[0]->consvToPrim(uPrev[i]);
        double r1 = uPrevPrim[1]/uPrevPrim[0];
        double r2 = uPrevPrim[2]/(1-uPrevPrim[0]);

        auto sgPtr = dynamic_cast<stiffenedGas*>(eos[0]);
        double cs1 = (eos[0]->get_gamma()[0]*(uPrevPrim[4]+sgPtr->get_pInf()[0])/(uPrevPrim[1]/uPrevPrim[0]));
        double cs2 = (eos[0]->get_gamma()[1]*(uPrevPrim[4]+sgPtr->get_pInf()[1])/(uPrevPrim[2]/(1-uPrevPrim[0])));

        double pref = uPrevPrim[0]*(1-uPrevPrim[0])*(r2*cs2-r1*cs1)/(uPrevPrim[1]*cs1+uPrevPrim[2]*cs2);

        uPlus1[i][FRAC] = u[i][FRAC] + (timestep/dx) * pref * (sStars[i-nGhost+1]-sStars[i-nGhost]); // flux[i + 1] and flux[i] for the update
        

        // Original allaire source term
        //uPlus1[i][FRAC] = u[i][FRAC] + (timestep/dx) * uPrev[i][FRAC] * (sStars[i-nGhost+1]-sStars[i-nGhost]); // flux[i + 1] and flux[i] for the update

        // Heun's method with Allaire source term
        //double K1 = (timestep/dx)* uPrev[i][FRAC] * (sStars[i-nGhost+1]-sStars[i-nGhost]);
        //double K2 = (timestep/dx)* (uPrev[i][FRAC] + K1) * (sStars[i-nGhost+1]-sStars[i-nGhost]);
        //uPlus1[i][FRAC] = u[i][FRAC] + 0.5*(K1+K2);

    }
    
    u = uPlus1;
    transmissiveBC();
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
    //result[FRAC] = arr[FRAC]*arr[UX]; // vol fraction 
    result[FRAC] = 0;
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
    for (size_t i = 0; i < halfSlopes.size(); ++i){ // halfslopes[i] = halfslopes i-1/2
        halfSlopes[i] = (prim == 0) ? u[i+1]-u[i] : uPrim[i+1]-uPrim[i];
    }
};

void solver::calcr(){
    for (size_t i = 0; i < r.size(); ++i){
        r[i] = elementDivide(halfSlopes[i],halfSlopes[i+1]);
    }
    //print_arr(r);
    std::cout << "calcd r" << std::endl;
};


void solver::calcUBars(bool prim){ // calculates for all u except leftmost cell
    for (size_t i = 0; i<uBarL.size(); ++i){
        uBarL[i] = (prim == 0) ? u[i+nGhost-1] - 0.5 * slopeLim(r[i]) * (0.5 * (halfSlopes[i]+halfSlopes[i+1])) : uPrim[i+nGhost-1] - 0.5 * slopeLim(r[i]) * (0.5 * (halfSlopes[i]+halfSlopes[i+1]));
        uBarR[i] = (prim == 0) ? u[i+nGhost-1] + 0.5 * slopeLim(r[i]) * (0.5 * (halfSlopes[i]+halfSlopes[i+1])) : uPrim[i+nGhost-1] + 0.5 * slopeLim(r[i]) * (0.5 * (halfSlopes[i]+halfSlopes[i+1]));
    }
    for (size_t i = 0; i < sStars.size(); ++i){
        sStars[i] = (prim == 0) ? 0.5*(eos[0]->consvToPrim(uBarR[i])[UX]+eos[0]->consvToPrim(uBarL[i+1])[UX]) : 0.5*(uBarR[i][UX]+uBarL[i+1][UX]);
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
        for (size_t i = 0; i<uBarLupd.size(); ++i){
            uBarLupd[i] = uBarL[i]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(uBarR[i]),eos[0])-flux(eos[0]->consvToPrim(uBarL[i]),eos[0]));
            uBarRupd[i] = uBarR[i]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(uBarR[i]),eos[0])-flux(eos[0]->consvToPrim(uBarL[i]),eos[0]));
        }
    } else {
        std::vector< std::vector <double> > B;
        for (size_t i = 0; i<uBarLupd.size(); ++i){
            double v = uPrim[i+nGhost-1][UX]; double rho = uPrim[i+nGhost-1][RHO1]+uPrim[i+nGhost-1][RHO2];
            B = {   {v, 0,  0,  0,                                            0},
                    {0, v,  0,  uPrim[i+nGhost-1][RHO1],                               0},
                    {0, 0,  v,  uPrim[i+nGhost-1][RHO2],                               0},
                    {0, 0,  0,  v,                                      1.0/rho},
                    {0, 0,  0,  rho*pow(eos[0]->calcSoundSpeed(uPrim[i+nGhost-1]),2)  ,v}  };

            std::array<double,5> deltaU = uBarR[i]-uBarL[i];
            uBarLupd[i] = eos[0]->primToConsv(uBarL[i]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
            uBarRupd[i] = eos[0]->primToConsv(uBarR[i]-0.5*(dt/dx)*multiplyMatrixVector(B,deltaU));
        }

    }
    uBarL = uBarLupd;
    uBarR = uBarRupd;
};

std::array<double,5> solver::calcSlope(double omega, std::array<double,5> slopeleft, std::array<double,5> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,5> solver::minbee(std::array<double,5> slp){
    std::array<double,5> minbArr{0,0,0,0,0};
    double minb;
    for (int k = 0; k<5; ++k){
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

std::array<double,5> solver::superbee(std::array<double,5> slp){
    std::array<double,5> minbArr{0,0,0,0,0};
    double minb;
    for (int k = 0; k<5; ++k){
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
    std::array<double,5> minbArr{0,0,0,0,0};
    for (int k = 0; k<5; ++k){
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
    for (std::vector<double>::size_type i = nGhost; i < u.size()-nGhost; ++i) {
        aMax = std::max(aMax, std::abs(eos[0]->consvToPrim(u[i])[UX]) + eos[0]->calcSoundSpeed(eos[0]->consvToPrim(u[i])));
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
solver::solver(double x_0, double x_1, double t0, double t1, int n_Cells, int n_Ghosts, double c){
    assert(t1>t0);
    x0 = x_0; x1 = x_1;
    startTime = t0; endTime = t1;
    timeMulti = (endTime-startTime <= 1e-3) ? 1000 : 1; // 1.0/(endTime-startTime);
    nCells = n_Cells;
    nGhost = n_Ghosts;
    dx = (x1 - x0)/nCells;

    phi.resize(nCells+2*nGhost);
    phiPlus1.resize(nCells+2*nGhost);

    uPlus1.resize(nCells+2*nGhost);

    u.resize(nCells+2*nGhost);
    uPrim.resize(nCells+2*nGhost);
    fluxes.resize(nCells+1);
    sStars.resize(nCells+1);
    halfSlopes.resize(nCells+3);
    r.resize(nCells+2); // r can only be defined for inner cells
    uBarL.resize(nCells+2);
    uBarR.resize(nCells+2);
    uBarLupd.resize(nCells+2);
    uBarRupd.resize(nCells+2);
    // fluid 2

    sourceResult.resize(nCells);

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
    return nGhost;
};

double solver::get_dx()const{
    return dx;
};


void solver::init(std::vector< std::array<double,5> > init) {
    assert(init.size() == u.size());
    u = init;
};

// -------------------- writing functions --------------------- //

void solver::setWriteInterval(double wt){
    writeInterval = wt;
}

void solver::writeData() const{
        // making data file
        std::string filename = dirName + "/" + std::to_string(time*timeMulti).substr(0,4);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        writeFile << "x alpha rho1 rho2 v p" << std::endl;
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (std::vector<double>::size_type i=nGhost; i<(u.size()-nGhost);i++){
                double x = x0 + (i-nGhost+0.5)*dx;
                std::array<double,5> prim = eos[0]->consvToPrim(u[i]);
                writeFile << x << " " << prim[0] << " " << prim[1] << " " << prim[2] << " " << prim[3] << " " << prim[4] << std::endl;
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



