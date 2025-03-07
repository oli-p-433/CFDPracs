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

        this->setDt();
        // dt = 0.1*dx;
        std::cout << "dt is" << dt << std::endl;

        transmissiveBC();

        //this->SLIC();
        //this->MUSCL();
        this->HLLCGodunov();
        uPrev = u;
        this->pointsUpdate();
        this->allaireSource();

        std::cout << "points updated" << std::endl;


        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData();
        }
        
    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData();
};


void solver::SLIC(){
    std::cout << "starting SLIC" << std::endl;

    calcHalfSlopes();

    calcr();

    calcUBars();

    updUBars();

    for (size_t i = 0; i < fluxes.size(); ++i) {
        fluxes[i] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(uBarR[i],uBarL[i+1],eos[0])),eos[0])+LF(uBarR[i],uBarL[i+1],eos[0])); // fluxes stored in conserved variable form
    }

    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
}

void solver::MUSCL(){
    std::cout << "starting MUSCL" << std::endl;

    calcHalfSlopes();

    calcr();

    calcUBars();

    updUBars();

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

    sStars[i] = sStar;
    //
    double prefL = (SL-LPrim[UX])/(SL-sStar);
    double prefR = (SR-RPrim[UX])/(SR-sStar);

    double eL,eR;
    eL = (left[ENE]/rhoL)+(sStar-LPrim[UX])*(sStar+LPrim[PRES]/(rhoL*(SL-LPrim[UX])));
    eR = (right[ENE]/rhoR)+(sStar-RPrim[UX])*(sStar+RPrim[PRES]/(rhoR*(SR-RPrim[UX])));


    std::array<double,5> uHLLCL,uHLLCR;
    uHLLCL = {1, left[RHO1]*prefL, left[RHO2]*prefL, rhoL*prefL*sStar, rhoL*prefL*eL}; // ?
    uHLLCR = {1, right[RHO1]*prefR, right[RHO2]*prefR, rhoR*prefR*sStar, rhoR*prefR*eR}; // ?

    
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

void solver::transmissiveBC(){
        for (int i = 0; i < nGhost; ++i){ // u BCs
            uPlus1[i] = u[nGhost];
            u[i] = u[nGhost];
            uPlus1[nCells+2*nGhost-1-i] = u[nCells+nGhost-1];
            u[nCells+2*nGhost-1-i] = u[nCells+nGhost-1];
        }
}








void solver::pointsUpdate(){
    transmissiveBC();

    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        uPlus1[i] = u[i] - (dt / dx) * (fluxes[i-nGhost+1] - fluxes[i-nGhost]); // flux[i + 1] and flux[i] for the update
    }
    u = uPlus1;
    // print_arr(u,1);
    std::cout << "updated u" << std::endl;
    // }
};

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

void solver::allaireSource(){
    transmissiveBC();
    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        uPlus1[i][FRAC] = u[i][FRAC] + (dt/dx) * uPrev[i][FRAC] * (sStars[i-nGhost+1]-sStars[i-nGhost]); // flux[i + 1] and flux[i] for the update
        std::cout << uPrev[i][FRAC] << " " << (sStars[i-nGhost+1]-sStars[i-nGhost]) << std::endl;
    }
    
    u = uPlus1;
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

void solver::calcHalfSlopes(){
    for (size_t i = 0; i < halfSlopes.size(); ++i){ // halfslopes[i] = halfslopes i-1/2
        halfSlopes[i] = u[i+1]-u[i];
    }
};

void solver::calcr(){
    for (size_t i = 0; i < r.size(); ++i){
        r[i] = elementDivide((u[i+1]-u[i]),(u[i+2]-u[i+1]));
    }
    //print_arr(r);
    std::cout << "calcd r" << std::endl;
};


void solver::calcUBars(){ // calculates for all u except leftmost cell
    for (size_t i = 0; i<uBarL.size(); ++i){
        uBarL[i] = u[i+nGhost-1] - 0.5 * slopeLim(r[i]) * (0.5 * (halfSlopes[i]+halfSlopes[i+1]));
        uBarR[i] = u[i+nGhost-1] + 0.5 * slopeLim(r[i]) * (0.5 * (halfSlopes[i]+halfSlopes[i+1]));
    }
};


void solver::updUBars(){
    for (size_t i = 0; i<uBarLupd.size(); ++i){
        uBarL[i] = uBarL[i]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(uBarR[i]),eos[0])-flux(eos[0]->consvToPrim(uBarL[i]),eos[0]));
        uBarR[i] = uBarR[i]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(uBarR[i]),eos[0])-flux(eos[0]->consvToPrim(uBarL[i]),eos[0]));
    }
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
    for (std::vector<double>::size_type i = nGhost; i < u.size()-nGhost; ++i) {
        std::cout << "sound speed " << eos[0]->calcSoundSpeed(u[i]) << std::endl;
        aMax = std::max(aMax, std::abs(eos[0]->consvToPrim(u[i])[UX]) + eos[0]->calcSoundSpeed(u[i]));
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
    timeMulti = 1; // 1.0/(endTime-startTime);
    nCells = n_Cells;
    nGhost = n_Ghosts;
    dx = (x1 - x0)/nCells;

    phi.resize(nCells+2*nGhost);
    phiPlus1.resize(nCells+2*nGhost);

    uPlus1.resize(nCells+2*nGhost);

    u.resize(nCells+2*nGhost);
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



