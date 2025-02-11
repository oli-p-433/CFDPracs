#include "solverGFM.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>

// constuctor

solverGFM::solverGFM(double x_0, double x_1, double t0, double t1, int n_Cells, int n_Ghosts, double c){
    assert(t1>t0);
    x0 = x_0; x1 = x_1;
    startTime = t0; endTime = t1;
    timeMulti = 1.0/(endTime-startTime);
    nCells = n_Cells;
    nGhost = n_Ghosts;
    dx = (x1 - x0)/nCells;

    phi.resize(nCells+2*nGhost);
    phiPlus1.resize(nCells+2*nGhost);

    temp_u.resize(nCells+2*nGhost);

    cour = c;

    std::vector<std::string> variables = {"rho","v","p"};

    //varMap["rho"] = 0; varMap["v"] = 1; varMap["p"] = 2;
    for (size_t i = 0; i<variables.size(); ++i){
        varMap[variables[i]] = i;
    }
    
};

void solverGFM::run(){
    std::cout << "pre-bc" << std::endl;

    transmissiveBC();
    std::cout << "set BCs" << std::endl;

    do{


        this->setDt(); //std::cout << "dt is" << dt << std::endl;

        this->levelSetUpdate();

        this->SLIC();
        
        this->pointsUpdate();
        std::cout << "points updated" << std::endl;

        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData("rho"); writeData("v"); writeData("p");
        }
        
    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData("rho"); writeData("v"); writeData("p");
};

// Level set evolution

void solverGFM::levelSetUpdate(){
    transmissiveBC();
    for (std::vector<double>::size_type i = nGhost; i < phi.size()-nGhost; ++i) {
        double uVal=0;
        if (phi[i]<=0){
            uVal = fluids[0]->eos->consvToPrim(fluids[0]->u[i])[1];
        } else if (phi[i] > 0) {
            uVal = fluids[1]->eos->consvToPrim(fluids[1]->u[i])[1];
        } else {
            throw std::runtime_error("level set undefined");
        }

        if (uVal >= 0){
            phiPlus1[i] = phi[i] - uVal*(dt/dx)*(phi[i]-phi[i-1]);
        } else {
            phiPlus1[i] = phi[i] - uVal*(dt/dx)*(phi[i+1]-phi[i]);
        }
        
    }

    phi = phiPlus1;
    transmissiveBC();

}

// State evolution

void solverGFM::pointsUpdate(){
    transmissiveBC();

    for (fluid* f : fluids){
        for (int i = nGhost; i < nCells+nGhost; ++i) {
            temp_u[i] = f->u[i] - (dt / dx) * (f->fluxes[i-nGhost+1] - f->fluxes[i-nGhost]); // flux[i + 1] and flux[i] for the update
        }
        f->u = temp_u;
    }

    transmissiveBC();
    // print_arr(u,1);
    std::cout << "updated u " << std::endl;
    // }
};

// Boundary conditions

void solverGFM::transmissiveBC(){
    for (fluid* f : fluids){
        for (int i = 0; i < nGhost; ++i){ // u BCs
            f->u[i] = f->u[nGhost];
            f->u[nCells+2*nGhost-1-i] = f->u[nCells+nGhost-1];
        }
    }
    for (int i = 0; i < nGhost; ++i){ // phi BCs
        phi[i] = phi[nGhost];
        phi[nCells+2*nGhost-1-i] = phi[nCells+nGhost-1];
    }
}

// Flux functions

std::array<double,3> solverGFM::fEuler(const EOS* eos, std::array<double,3> arr){ // takes primitive variables
    std::array<double,3> result;
    result[0] = arr[0]*arr[1]; // rho*v
    result[1] = arr[0]*arr[1]*arr[1]+arr[2]; // rho*v^2 + p
    result[2] = (eos->primToConsv(arr)[2] + arr[2])*arr[1]; // (E + p)v
    return result;
}

// ------------------------- secondary fluxes --------------------- //

std::array<double,3> solverGFM::LF(const EOS* eos, std::array<double,3> v1, std::array<double,3> v2)const{
    return 0.5 * (flux(eos,eos->consvToPrim(v1)) + flux(eos,eos->consvToPrim(v2))) + 0.5 * (dx / dt) * (v1 - v2);
};

std::array<double,3> solverGFM::laxFhalf(const EOS* eos, std::array<double,3> ui, std::array<double,3> uip1)const{
    return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos,eos->consvToPrim(uip1))-flux(eos,eos->consvToPrim(ui))));
};

// Slope limited scheme

void solverGFM::SLIC(){
    for (fluid* f : fluids){
        std::cout << "starting SLIC" << std::endl;

        calcHalfSlopes(f);

        calcr(f);

        calcUBars(f);

        updUBars(f);

        for (size_t i = 0; i < f->fluxes.size(); ++i) {
            f->fluxes[i] = 0.5*(flux(f->eos,f->eos->consvToPrim(laxFhalf(f->eos,f->uBarRupd[i],f->uBarLupd[i+1])))+LF(f->eos,f->uBarRupd[i],f->uBarLupd[i+1])); // fluxes stored in conserved variable form
            //f->fluxes[i] = LF(f->eos,f->u[i+nGhost],f->u[i+nGhost+1]);
            /*
            if (i > 40 && i < 60){
                std::cout << f->fluxes[i][0] << std::endl;
            }
            */
        }

        std::cout << std::endl << "-----------------------fluxes above-------------------------" << std::endl;
    }
}

void solverGFM::calcHalfSlopes(fluid* f){
    for (size_t i = 0; i < f->halfSlopes.size(); ++i){
        f->halfSlopes[i] = f->u[i+nGhost-1]-f->u[i+nGhost-2];
    }
};

void solverGFM::calcr(fluid* f){
    for (size_t i = 0; i < f->r.size(); ++i){
        f->r[i] = elementDivide(f->halfSlopes[i],f->halfSlopes[i+1]);
    };
    //print_arr(r);
    std::cout << "calcd r" << std::endl;
};

void solverGFM::calcUBars(fluid* f){ // calculates for all u except leftmost cell
    for (size_t i = 0; i<f->uBarL.size(); ++i){
        f->uBarL[i] = f->u[i+nGhost-1] - 0.5 * ( slopeLim(f->r[i]) * calcSlope(f->get_slopeWeight(), f->halfSlopes[i], f->halfSlopes[i+1]) );
        f->uBarR[i] = f->u[i+nGhost-1] + 0.5 * ( slopeLim(f->r[i]) * calcSlope(f->get_slopeWeight(), f->halfSlopes[i], f->halfSlopes[i+1]) );
    }
};

void solverGFM::updUBars(fluid* f){
    for (size_t i = 0; i<f->uBarLupd.size(); ++i){
        f->uBarLupd[i] = f->uBarL[i]-0.5*(dt/dx)*(flux(f->eos,f->eos->consvToPrim(f->uBarR[i]))-flux(f->eos,f->eos->consvToPrim(f->uBarL[i])));
        f->uBarRupd[i] = f->uBarR[i]-0.5*(dt/dx)*(flux(f->eos,f->eos->consvToPrim(f->uBarR[i]))-flux(f->eos,f->eos->consvToPrim(f->uBarL[i])));
    }
};

std::array<double,3> solverGFM::calcSlope(double omega, std::array<double,3> slopeleft, std::array<double,3> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,3> solverGFM::minbee(std::array<double,3> slp){
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
    //std::cout << "Minbee value is" << minBval << std::endl;
    //return {minbArr[3],minbArr[3],minbArr[3],minbArr[3]};
    return {minBval,minBval,minBval};
    //return minbArr;
};

void solverGFM::setDt(){ // only works for two fluids with phi < 0 and phi > 0 corresponding to the physical fluids, respectively
    double aMax = -1e14;
    std::cout << "setting dt" << std::endl;

    for (int i = nGhost; i < nCells+nGhost; ++i){
        if (phi[i]<0){
            aMax = std::max(aMax, std::fabs(fluids[0]->eos->consvToPrim(fluids[0]->u[i])[1]));
        } else if (phi[i]>=0) {
            aMax = std::max(aMax, std::fabs(fluids[1]->eos->consvToPrim(fluids[1]->u[i])[1]));
        } else {
            writeData("rho"); writeData("v"); writeData("p");
            throw std::runtime_error("level set value undefined");

        }
    }

    dt = cour * dx / std::fabs(aMax);

    std::cout << "set dt = " << dt << std::endl;

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

// utilities

int solverGFM::ghosts(){
    return nGhost;
};

double solverGFM::get_dx()const{
    return dx;
};

void solverGFM::init(const std::vector< std::array<double,3> > &init, fluid* f) {
    if (!f) {
        throw std::runtime_error("Null pointer passed to init!");
    }
    f->u = init;
};

void solverGFM::phiInit(std::vector< double > init) {
    assert(init.size() == phi.size());
    phi = init;
};

// -------------------- writing functions --------------------- //

void solverGFM::setWriteInterval(double wt){
    writeInterval = wt;
}

void solverGFM::writeData(std::string varName) const{
        // making data file
        std::string filename = dirName + "/" + varName + "/" + std::to_string(time).substr(0,4);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (int i=nGhost; i<nCells+nGhost;i++){
                double x = x0 + (i-nGhost+0.5)*dx;
                writeFile << x << " " << phi[i] << " ";
                for (const fluid* f : fluids){
                    writeFile << f->eos->consvToPrim(f->u[i])[varMap.at(varName)] << " ";
                }
                writeFile << std::endl;
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



