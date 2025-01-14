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
        this->SLIC();
        this->pointsUpdate();
        std::cout << "points updated" << std::endl;

        this->sourceUpdate();

        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData("rho"); writeData("v"); writeData("p");
        }
        
    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData("rho"); writeData("v"); writeData("p");
};

void solver::laxFriedrichs(){
    for (int i = 0; i < nCells + 1; ++i) {

        fluxes[i] = 0.5 * (flux(eos->consvToPrim(u[i])) + flux(eos->consvToPrim(u[i + 1]))) + 0.5 * (dx / dt) * (u[i] - u[i+1]);

    }
    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
};

void solver::Richt(){
    for (int i = 0; i < nCells + 1; ++i) {
        fluxes[i] = flux(eos->consvToPrim(laxFhalf(u[i],u[i+1])));
    }
    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
};

void solver::FORCE(){ // problem is probably with ghost cells / indexing

    for (size_t i = 0; i < fluxes.size(); ++i) {
        //std::cout << "finvect " << finvectLF(u[i]) << " " << finvectLF(u[i]) << std::endl;
        fluxes[i] = 0.5*(flux(eos->consvToPrim(laxFhalf(u[i+nGhost-1],u[i+nGhost])))+LF(u[i+nGhost-1],u[i+nGhost])); // fluxes stored in conserved variable form
    }
    // print_arr(fluxes,0);
    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
};

void solver::SLIC(){
    std::cout << "starting SLIC" << std::endl;

    calcHalfSlopes();

    calcr();

    calcUBars();

    updUBars();

    for (size_t i = 0; i < fluxes.size(); ++i) {
        fluxes[i] = 0.5*(flux(eos->consvToPrim(laxFhalf(uBarRupd[i],uBarLupd[i+1])))+LF(uBarRupd[i],uBarLupd[i+1])); // fluxes stored in conserved variable form
    }

    std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
}

void solver::transmissiveBC(){
    if (alpha == 0){
        for (int i = 0; i < nGhost; ++i){
            uPlus1[i] = u[nGhost];
            u[i] = u[nGhost];
            uPlus1[nCells+2*nGhost-1-i] = u[nCells+nGhost-1];
            u[nCells+2*nGhost-1-i] = u[nCells+nGhost-1];
        }
    } else {
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
        }
    }
}

void solver::pointsUpdate(){
    transmissiveBC();

    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        uPlus1[i] = u[i] - (dt / dx) * (fluxes[i-nGhost+1] - fluxes[i-nGhost]); // flux[i + 1] and flux[i] for the update
    }
    
    std::cout << std::endl << "-------------------------uplus1 above------------------------" << std::endl;
    // for (int i = 0; i < u.size(); ++i) {
    u = uPlus1;
    // print_arr(u,1);
    std::cout << "updated u" << std::endl;
    // }
};

// ---------------- Source Terms ---------------------------- //

void solver::sourceUpdate(){
    transmissiveBC();
    std::vector< std::array<double,3> > sourceArr = sourceTerm(u,alpha);
    std::cout << "printing source array" << std::endl;
    // print_arr(sourceArr,0);

    for (std::vector<double>::size_type i = nGhost; i < u.size() - nGhost; ++i) {
        uPlus1[i] = u[i] + dt * sourceArr[i-nGhost]; // flux[i + 1] and flux[i] for the update

        //print_vect(uPlus1[i]);
    }
    
    std::cout << std::endl << "-------------------------source terms above------------------------" << std::endl;
    // for (int i = 0; i < u.size(); ++i) {
    u = uPlus1;
    std::cout << "evolved source terms" << std::endl;
    // }
};

std::vector< std::array<double,3> > solver::sourceTerm(const std::vector< std::array<double,3> > arr, int alpha){

    for (size_t i = nGhost; i<(arr.size()-nGhost); ++i){
        double pref = -(alpha/((i-nGhost+0.5)*dx));
        sourceResult[i-nGhost][0] = pref*arr[i][1]; // rho*v
        sourceResult[i-nGhost][1] = pref*(arr[i][1]*arr[i][1])/arr[i][0]; // (rho*v)^2 / rho
        sourceResult[i-nGhost][2] = pref*(eos->consvToPrim(arr[i])[1])*(arr[i][2]+eos->consvToPrim(arr[i])[2]);
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

std::array<double,3> solver::fEuler(std::array<double,3> arr){ // takes primitive variables
    std::array<double,3> result;
    result[0] = arr[0]*arr[1]; // rho*v
    result[1] = arr[0]*arr[1]*arr[1]+arr[2]; // rho*v^2 + p
    result[2] = (eos->primToConsv(arr)[2] + arr[2])*arr[1]; // (E + p)v
    return result;
}

// ------------------------- secondary fluxes --------------------- //

std::array<double,3> solver::LF(std::array<double,3> v1, std::array<double,3> v2)const{
    return 0.5 * (flux(eos->consvToPrim(v1)) + flux(eos->consvToPrim(v2))) + 0.5 * (dx / dt) * (v1 - v2);
};

std::array<double,3> solver::laxFhalf(std::array<double,3> ui, std::array<double,3> uip1)const{
    return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos->consvToPrim(uip1))-flux(eos->consvToPrim(ui))));
};



// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(){
    for (size_t i = 0; i < halfSlopes.size(); ++i){
        halfSlopes[i] = u[i+nGhost-1]-u[i+nGhost-2];
    }
};

void solver::calcr(){
    for (size_t i = 0; i < r.size(); ++i){
        r[i] = elementDivide(halfSlopes[i],halfSlopes[i+1]);
    }
    //print_arr(r);
    std::cout << "calcd r" << std::endl;
};


void solver::calcUBars(){ // calculates for all u except leftmost cell
    for (size_t i = 0; i<uBarL.size(); ++i){
        uBarL[i] = u[i+nGhost-1] - 0.5 * ( slopeLim(r[i]) * calcSlope(slopeWeight, halfSlopes[i], halfSlopes[i+1]) );
        uBarR[i] = u[i+nGhost-1] + 0.5 * ( slopeLim(r[i]) * calcSlope(slopeWeight, halfSlopes[i], halfSlopes[i+1]) );
    }
};


void solver::updUBars(){
    for (size_t i = 0; i<uBarLupd.size(); ++i){
        uBarLupd[i] = uBarL[i]-0.5*(dt/dx)*(flux(eos->consvToPrim(uBarR[i]))-flux(eos->consvToPrim(uBarL[i])));
        uBarRupd[i] = uBarR[i]-0.5*(dt/dx)*(flux(eos->consvToPrim(uBarR[i]))-flux(eos->consvToPrim(uBarL[i])));
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
    for (std::vector<double>::size_type i = nGhost; i < u.size()-nGhost; ++i) {
        aMax = std::max(aMax, std::fabs(eos->consvToPrim(u[i])[1])+eos->calcSoundSpeed(eos->consvToPrim(u[i])));
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
    u.resize(nCells+2*nGhost);
    uPlus1.resize(nCells+2*nGhost);
    fluxes.resize(nCells+1);
    halfSlopes.resize(nCells+3);
    r.resize(nCells+2); // r can only be defined for inner cells
    uBarL.resize(nCells+2);
    uBarR.resize(nCells+2);
    uBarLupd.resize(nCells+2);
    uBarRupd.resize(nCells+2);

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

void solver::init(std::vector< std::array<double,3> > init) {
    assert(init.size() == u.size());
    u = init;
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
            for (std::vector<double>::size_type i=nGhost; i<(u.size()-nGhost);i++){
                double x = x0 + (i-nGhost+0.5)*dx;
                writeFile << x << " " << eos->consvToPrim(u[i])[varMap.at(varName)] << std::endl;
                //std::cout << u[i] << " ";
            }
            writeFile.close();
            std::cout << filename << "written successfully" << std::endl;
        }else{
            std::cout << "failed to write." << std::endl;
        }

}

// ------------------- Overloading operators ------------------ //

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



