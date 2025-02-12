#include "solver.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>


void solver::run(){
    std::cout << "pre-bc" << std::endl;

    transmissiveBC(fluid1);
    transmissiveBC(fluid2);

    int splitFlip = 0;
    do{
        this->setDt();

        std::cout << "dt is" << dt << std::endl;
        std::cout << "time is " << time << std::endl;

        transmissiveBC(fluid1);
        transmissiveBC(fluid2);
        std::cout << "setting phiBC" << std::endl;
        phiBC();
        std::cout << "updating phi" << std::endl;

        phiUpdate();

        splitFlip++;
        if (splitFlip%2 == 0){
            direction = XDIR;
            this->MUSCL(fluid1);
            this->pointsUpdate(fluid1);
            this->MUSCL(fluid2);
            this->pointsUpdate(fluid2);
            direction = YDIR;
            this->MUSCL(fluid1);
            this->pointsUpdate(fluid1);
            this->MUSCL(fluid2);
            this->pointsUpdate(fluid2);
        } else {
            direction = YDIR;
            this->MUSCL(fluid1);
            this->pointsUpdate(fluid1);
            this->MUSCL(fluid2);
            this->pointsUpdate(fluid2);
            direction = XDIR;
            this->MUSCL(fluid1);
            this->pointsUpdate(fluid1);
            this->MUSCL(fluid2);
            this->pointsUpdate(fluid2);
        };
        std::cout << "points updated" << std::endl;

        transmissiveBC(fluid1);
        transmissiveBC(fluid2);

        //this->sourceUpdate();

        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        }

    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
};

void solver::SLIC(fluid& fluid){

    calcHalfSlopes(fluid); // -------------------------------------------- WIP

    calcr(fluid);

    calcUBars(fluid);

    updUBars(fluid); //-----------!!!!!!!!!!!!!!!!

    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1])),eos[0])+LF(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1])); // fluxes stored in conserved variable form
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j])),eos[0])+LF(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j])); // fluxes stored in conserved variable form
            }
        }
    }
}

std::array<double,4> solver::riemannSolver(std::array<double,4> left,std::array<double,4> right){
    riemann solution(eos[0]->get_gamma(),eos[0]->get_gamma(),eos[0]->consvToPrim(left),eos[0]->consvToPrim(right),direction,0,1,0.5,1,2,0,0);
    return solution.exctRiemann();
}

void solver::MUSCL(fluid& fluid){

    calcHalfSlopes(fluid); 

    calcr(fluid);

    calcUBars(fluid);

    updUBars(fluid);

    if (direction == XDIR){
        for (size_t i = 0; i < fluid.fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluid.fluxesX[0].size(); ++j) {
                fluid.fluxesX[i][j] = flux(riemannSolver(fluid.uBarRX[i][j],fluid.uBarLX[i][j+1]),eos[0]);
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluid.fluxesY.size(); ++i) { // why am i getting a Y flux of vx? 
            for (size_t j = 0; j < fluid.fluxesY[0].size(); ++j) {
                fluid.fluxesY[i][j] = flux(riemannSolver(fluid.uBarRY[i][j],fluid.uBarLY[i+1][j]),eos[0]);
            }
        }
    }
}



void solver::transmissiveBC(fluid& fluid){
    for (int i = 0; i < nGhost; ++i){
        for (size_t j=0; j<fluid.uPlus1.size(); ++j){ // sets rows (i)
            fluid.uPlus1[i][j] = fluid.u[nGhost][j];
            fluid.u[i][j] = fluid.u[nGhost][j];
            fluid.uPlus1[nCells+2*nGhost-1-i][j] = fluid.u[nCells+nGhost-1][j];
            fluid.u[nCells+2*nGhost-1-i][j] = fluid.u[nCells+nGhost-1][j];
        }
    }
    for (size_t i=0; i<fluid.uPlus1.size(); ++i){
        for (int j = 0; j < nGhost; ++j){
            fluid.uPlus1[i][j] = fluid.u[i][nGhost];
            fluid.u[i][j] = fluid.u[i][nGhost];
            fluid.uPlus1[i][nCells+2*nGhost-1-j] = fluid.u[i][nCells+nGhost-1];
            fluid.u[i][nCells+2*nGhost-1-j] = fluid.u[i][nCells+nGhost-1];
        }
    }
}

void solver::cylTransmissiveBC(fluid& fluid){
    for (int var =0; var<4;++var){
        if (var == 0 || var == 2 || var == 3){
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<fluid.uPlus1.size(); ++j){ // sets rows (ie z)
                    fluid.uPlus1[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.u[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.uPlus1[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                    fluid.u[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<fluid.uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    fluid.uPlus1[i][nGhost-1-j][var] = fluid.u[i][nGhost+j][var];
                    fluid.u[i][nGhost-1-j][var] = fluid.u[i][nGhost+j][var];
                    fluid.uPlus1[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                    fluid.u[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                }
            }
        } else {
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<fluid.uPlus1.size(); ++j){ // sets rows (ie z)
                    fluid.uPlus1[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.u[i][j][var] = fluid.u[nGhost][j][var];
                    fluid.uPlus1[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                    fluid.u[nCells+2*nGhost-1-i][j][var] = fluid.u[nCells+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<fluid.uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    fluid.uPlus1[i][nGhost-1-j][var] = (-1)*fluid.u[i][nGhost+j][var];
                    fluid.u[i][nGhost-1-j][var] = (-1)*fluid.u[i][nGhost+j][var];
                    fluid.uPlus1[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                    fluid.u[i][nCells+2*nGhost-1-j][var] = fluid.u[i][nCells+nGhost-1][var];
                }
            }
        }
    }
}

void solver::phiBC(){
    for (int i = 0; i < nGhost; ++i){
        for (size_t j=0; j<phi.size(); ++j){ // sets rows (i)
            phiPlus1[i][j] = phi[nGhost][j];
            phi[i][j] = phi[nGhost][j];
            phiPlus1[nCells+2*nGhost-1-i][j] = phi[nCells+nGhost-1][j];
            phi[nCells+2*nGhost-1-i][j] = phi[nCells+nGhost-1][j];
        }
    }
    for (size_t i=0; i<phi.size(); ++i){
        for (int j = 0; j < nGhost; ++j){
            phiPlus1[i][j] = phi[i][nGhost];
            phi[i][j] = phi[i][nGhost];
            phiPlus1[i][nCells+2*nGhost-1-j] = phi[i][nCells+nGhost-1];
            phi[i][nCells+2*nGhost-1-j] = phi[i][nCells+nGhost-1];
        }
    }
}

void solver::phiUpdate(){
    phiBC();

    double uVal;

    for (std::vector<double>::size_type i = nGhost-1; i < phi.size()-nGhost+1; ++i){
        for (std::vector<double>::size_type j = nGhost; j < phi[0].size() - nGhost; ++j) {
            if (phi[i][j] >= 0){
                uVal = eos[0]->consvToPrim(fluid1.u[i][j])[UX];
            } else {
                uVal = eos[1]->consvToPrim(fluid2.u[i][j])[UX];
            }

            if (uVal >= 0){
                phiPlus1[i][j] = phi[i][j] - uVal * (dt / dx) * (phi[i][j] - phi[i][j-1]); // flux[i + 1] and flux[i] for the update
            } else {
                phiPlus1[i][j] = phi[i][j] - uVal * (dt / dx) * (phi[i][j+1] - phi[i][j]); // flux[i + 1] and flux[i] for the update
            }
        }
    }
    phi = phiPlus1;
    //print_arr(u,RHO);
    //std::cout << "phiX updated" << std::endl;
    for (std::vector<double>::size_type i = nGhost; i < phi.size()-nGhost; ++i){ // this is the part that causes vx to grow
        for (std::vector<double>::size_type j = nGhost-1; j < phi[0].size() - nGhost+1; ++j) {
            if (phi[i][j] >= 0){
                uVal = eos[0]->consvToPrim(fluid1.u[i][j])[UY];
            } else {
                uVal = eos[1]->consvToPrim(fluid2.u[i][j])[UY];
            }
            
            if (uVal >= 0){
                phiPlus1[i][j] = phi[i][j] - uVal * (dt / dy) * (phi[i][j] - phi[i-1][j]); // flux[i + 1] and flux[i] for the update
            } else {
                phiPlus1[i][j] = phi[i][j] - uVal * (dt / dy) * (phi[i+1][j] - phi[i][j]); // flux[i + 1] and flux[i] for the update
            }
        }
    }
    phi = phiPlus1;
    //print_arr(u,RHO);
    //std::cout << "Y updated" << std::endl;
}

void solver::pointsUpdate(fluid& fluid){ //second index = x, first index = y 
    //cylTransmissiveBC();
    if (direction == XDIR){
        for (std::vector<double>::size_type i = nGhost-1; i < fluid.u.size()-nGhost+1; ++i){
            for (std::vector<double>::size_type j = nGhost; j < fluid.u[0].size() - nGhost; ++j) {
                fluid.uPlus1[i][j] = fluid.u[i][j] - (dt / dx) * (fluid.fluxesX[i-nGhost+1][j-nGhost+1] - fluid.fluxesX[i-nGhost+1][j-nGhost]); // flux[i + 1] and flux[i] for the update

            }
        }
        fluid.u = fluid.uPlus1;
        //print_arr(u,RHO);
        std::cout << "X updated" << std::endl;
    } else if (direction == YDIR){
        for (std::vector<double>::size_type i = nGhost; i < fluid.u.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost-1; j < fluid.u[0].size() - nGhost+1; ++j) {
                fluid.uPlus1[i][j] = fluid.u[i][j] - (dt / dy) * (fluid.fluxesY[i-nGhost+1][j-nGhost+1] - fluid.fluxesY[i-nGhost][j-nGhost+1]); // flux[i + 1] and flux[i] for the update
            }
        }
        fluid.u = fluid.uPlus1;
        //print_arr(u,RHO);
        std::cout << "Y updated" << std::endl;
    }
    
};

// -------------------- Flux functions --------------------- //

std::array<double,4> solver::finvectLF(const std::array<double,4> x)const{
    double a=1;
    return a*x;
};

std::array<double,4> solver::fBurgersLF(const std::array<double,4> x)const{
    return 0.5*(x * x);
};

std::array<double,4> solver::fEuler(std::array<double,4> arr, EOS* eos){ // takes primitive variables
    std::array<double,4> result;
    if (direction == XDIR){
        result[RHO] = arr[RHO]*arr[UX]; // rho*v
        result[XMOM] = arr[RHO]*arr[UX]*arr[UX]+arr[PRES]; // rho*v^2 + p
        result[YMOM] = arr[RHO]*arr[UX]*arr[UY]; // rho*v^2 + p
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UX]; // (E + p)v
    } else if (direction == YDIR){
        result[RHO] = arr[RHO]*arr[UY]; // rho*v
        result[XMOM] = arr[RHO]*arr[UX]*arr[UY]; // rho*v^2 + p
        result[YMOM] = arr[RHO]*arr[UY]*arr[UY]+arr[PRES]; // rho*v^2 + p
        result[ENE] = (eos->primToConsv(arr)[ENE] + arr[PRES])*arr[UY]; // (E + p)v
    }

    return result;
}

// ------------------------- secondary fluxes --------------------- //

std::array<double,4> solver::LF(std::array<double,4> v1, std::array<double,4> v2)const{
    if (direction == XDIR){
        return 0.5 * (flux(eos[0]->consvToPrim(v1),eos[0]) + flux(eos[0]->consvToPrim(v2),eos[0])) + 0.5 * (dx / dt) * (v1 - v2);
    } else {
        return 0.5 * (flux(eos[0]->consvToPrim(v1),eos[0]) + flux(eos[0]->consvToPrim(v2),eos[0])) + 0.5 * (dy / dt) * (v1 - v2);
    }
};

std::array<double,4> solver::laxFhalf(std::array<double,4> ui, std::array<double,4> uip1)const{
    if (direction == XDIR){
        return 0.5*(ui+uip1)-(0.5*(dt/dx)*(flux(eos[0]->consvToPrim(uip1),eos[0])-flux(eos[0]->consvToPrim(ui),eos[0])));
    } else if (direction == YDIR){
        return 0.5*(ui+uip1)-(0.5*(dt/dy)*(flux(eos[0]->consvToPrim(uip1),eos[0])-flux(eos[0]->consvToPrim(ui),eos[0])));
    } else {
        throw std::runtime_error("laxfhalf doesnt know which dimension");
    }
};



// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(fluid& fluid){ // only using y-gradient - ie correspondng to y-faces; need extra column for x direction //DONE
    //std::cout << "calculating halfslopes n printing u" << std::endl;
    ////print_arr(u,RHO);
    if (direction == XDIR){
        for (size_t i = 0; i < fluid.halfSlopesX.size(); ++i){
            for (size_t j = 0; j < fluid.halfSlopesX[0].size(); ++j){
                fluid.halfSlopesX[i][j] = fluid.u[i+nGhost-1][j+nGhost-1]-fluid.u[i+nGhost-1][j+nGhost-2];
            }
        }
    } else {
        for (size_t i = 0; i < fluid.halfSlopesY.size(); ++i){
            for (size_t j = 0; j < fluid.halfSlopesY[0].size(); ++j){
                fluid.halfSlopesY[i][j] = fluid.u[i+nGhost-1][j+nGhost-1]-fluid.u[i+nGhost-2][j+nGhost-1];
            }
        }
    };
};

void solver::calcr(fluid& fluid){ // DONE
    if (direction == XDIR){
        for (size_t i = 0; i < fluid.rX.size(); ++i){
            for (size_t j = 0; j < fluid.rX[0].size(); ++j){
                fluid.rX[i][j] = elementDivide(fluid.halfSlopesX[i][j],fluid.halfSlopesX[i][j+1]);
            }
        }
    } else {
        for (size_t i = 0; i < fluid.rY.size(); ++i){
            for (size_t j = 0; j < fluid.rY[0].size(); ++j){
                fluid.rY[i][j] = elementDivide(fluid.halfSlopesY[i][j],fluid.halfSlopesY[i+1][j]);
            }
        }
    };
    ////print_arr(r);
    //std::cout << "calcd r" << std::endl;
};


void solver::calcUBars(fluid& fluid){ // calculates for all u except leftmost cell
    if (direction == XDIR){
        for (size_t i = 0; i<fluid.uBarLX.size(); ++i){
            for (size_t j = 0; j<fluid.uBarLX[0].size(); ++j){
            
                //print_vect(slopeLim(r[i-1]));
                //std::cout << std::endl << "slope ";
                //print_vect(calcSlope(slopeWeight, halfSlopes[i-1], halfSlopes[i]));
                //std::cout << std::endl;
                fluid.uBarLX[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] - 0.5 * ( slopeLim(fluid.rX[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesX[i][j], fluid.halfSlopesX[i][j+1]) );
                fluid.uBarRX[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] + 0.5 * ( slopeLim(fluid.rX[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesX[i][j], fluid.halfSlopesX[i][j+1]) );
                //std::cout << slopeLim(rX[i][j])[RHO] << " ";
            }
            //std::cout << std::endl;
        }
    } else {
        for (size_t i = 0; i<(fluid.uBarLY.size()); ++i){
            for (size_t j = 0; j<(fluid.uBarLY[0].size()); ++j){
            
                //print_vect(slopeLim(r[i-1]));
                //std::cout << std::endl << "slope ";
                //print_vect(calcSlope(slopeWeight, halfSlopes[i-1], halfSlopes[i]));
                //std::cout << std::endl;
                fluid.uBarLY[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] - 0.5 * ( slopeLim(fluid.rY[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesY[i][j], fluid.halfSlopesY[i+1][j]) );
                fluid.uBarRY[i][j] = fluid.u[i+nGhost-1][j+nGhost-1] + 0.5 * ( slopeLim(fluid.rY[i][j]) * calcSlope(slopeWeight, fluid.halfSlopesY[i][j], fluid.halfSlopesY[i+1][j]) );
                //std::cout << slopeLim(rY[i][j])[RHO] << " ";
            }
            //std::cout << std::endl;
        }
    };
};


void solver::updUBars(fluid& fluid){
    if (direction == XDIR){
        for (size_t i = 0; i<(fluid.uBarLupdX.size()); ++i){
            for (size_t j = 0; j<(fluid.uBarLupdX[0].size()); ++j){
                fluid.uBarLupdX[i][j] = fluid.uBarLX[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(fluid.uBarRX[i][j]),eos[0])-flux(eos[0]->consvToPrim(fluid.uBarLX[i][j]),eos[0]));
                fluid.uBarRupdX[i][j] = fluid.uBarRX[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(fluid.uBarRX[i][j]),eos[0])-flux(eos[0]->consvToPrim(fluid.uBarLX[i][j]),eos[0]));
            }
        }
        fluid.uBarLX = fluid.uBarLupdX;
        fluid.uBarRX = fluid.uBarRupdX;
    } else if (direction == YDIR){
        for (size_t i = 0; i<(fluid.uBarLupdY.size()); ++i){
            for (size_t j = 0; j<(fluid.uBarLupdY[0].size()); ++j){
                fluid.uBarLupdY[i][j] = fluid.uBarLY[i][j]-0.5*(dt/dy)*(flux(eos[0]->consvToPrim(fluid.uBarRY[i][j]),eos[0])-flux(eos[0]->consvToPrim(fluid.uBarLY[i][j]),eos[0]));
                fluid.uBarRupdY[i][j] = fluid.uBarRY[i][j]-0.5*(dt/dy)*(flux(eos[0]->consvToPrim(fluid.uBarRY[i][j]),eos[0])-flux(eos[0]->consvToPrim(fluid.uBarLY[i][j]),eos[0]));
            }
        }
        fluid.uBarLY = fluid.uBarLupdY;
        fluid.uBarRY = fluid.uBarRupdY;
    }
};

std::array<double,4> solver::calcSlope(double omega, std::array<double,4> slopeleft, std::array<double,4> sloperight){
        return 0.5*(1+omega)*slopeleft+0.5*(1-omega)*sloperight;
};

std::array<double,4> solver::minbee(std::array<double,4> slp){
    std::array<double,4> minbArr{0,0,0,0};
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


std::array<double,4> solver::vanLeer(std::array<double,4> slp){
    std::array<double,4> minbArr{0,0,0,0};
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

void solver::print_arr(std::vector < std::vector< std::array<double,4> > > arr, int var) { 
    std::string row;
    for (size_t j = 0; j < arr.size(); ++j) {
        for (size_t i = 0; i < arr[j].size(); ++i) {
            row += " ";

            // Use std::ostringstream to format the number to 6 decimal places
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << arr[j][i][var];

            // Append the formatted string to the row
            row += oss.str();
        }
        std::cout << row << std::endl;
        row.clear(); // Clear the row for the next line
    }
};

void solver::print_vect(std::vector < std::array<double,4> > v){ // just prints a slice (0 by default)
    for (int j = 0; j<4; ++j){
        std::cout << v[0][j] << " ";
        //std::cout << "laxFhalf " << laxFhalf(u[i],u[i+1])[j] << std::endl;
    }
    std::cout << std::endl;
}

void solver::print_state(std::array<double,4> v){ // just prints a slice (0 by default)
    for (int j = 0; j<4; ++j){
        std::cout << v[j] << " ";
        //std::cout << "laxFhalf " << laxFhalf(u[i],u[i+1])[j] << std::endl;
    }
    std::cout << std::endl;
}

void solver::resize2D(int nx, int ny, std::vector< std::vector< double > >& arr){
    std::cout << "Resizing vector to rows: " << nx << ", cols: " << ny << std::endl;

    arr.resize(ny);
    for (size_t i = 0; i < arr.size(); ++i){
        arr[i].resize(nx);
    }

}


// ---------------------- Time Stepping ---------------------- //

void solver::setDt(){
    double aMax = -1e14;
    std::cout << "setting dt" << std::endl;
    for (std::vector<double>::size_type i = 0; i < fluid1.u.size(); ++i) {
        for (std::vector<double>::size_type j = 0; j < fluid1.u.size(); ++j) {
            //aMax = std::max(aMax, std::fabs(std::sqrt((eos[0]->consvToPrim(fluid1.u[i][j])[UX])*(eos[0]->consvToPrim(fluid1.u[i][j])[UX])+(eos[0]->consvToPrim(fluid1.u[i][j])[UY])*(eos[0]->consvToPrim(fluid1.u[i][j])[UY])))+eos[0]->calcSoundSpeed(eos[0]->consvToPrim(fluid1.u[i][j])));
            aMax = std::max(aMax, std::max(std::sqrt(pow(eos[1]->consvToPrim(fluid2.u[i][j])[UX],2)+pow(eos[1]->consvToPrim(fluid2.u[i][j])[UY],2))+eos[1]->calcSoundSpeed(eos[1]->consvToPrim(fluid2.u[i][j])), std::sqrt(pow(eos[0]->consvToPrim(fluid1.u[i][j])[UX],2)+pow(eos[0]->consvToPrim(fluid1.u[i][j])[UY],2))+eos[0]->calcSoundSpeed(eos[0]->consvToPrim(fluid1.u[i][j]))));

        }
    }
    std::cout << "aMax = " << aMax << std::endl;
    dt = cour * std::min(dx,dy) / std::fabs(aMax);

    std::cout << "set dt = " << dt << std::endl;


    double writeTime = writeInterval*(std::floor((time+1e-9) / writeInterval)+1); //returing 0.29?
    std::cout << "writeTime = " << writeTime << std::endl;

    if (endTime - time <= dt && writeTime-time >= dt){
        std::cout << "option 0 " << std::endl;
        dt = endTime - time;
        assert(dt != 0);
        time = endTime;
        std::cout << "last step:" << dt << std::endl;
    } else if (writeTime-time <= dt && writeTime-time > 0){ // manually reducing dt to hit the write times
        std::cout << "option 1 " << std::endl;

        dt = writeTime - time;
        assert(dt != 0);
        time = writeTime;
    } else {
        std::cout << "option 2 " << std::endl;
        time = time + dt;
    }

    if (std::fabs(writeTime - time) < 1e-10){
        checkWrite = 1;
    } else {
        checkWrite = 0;
    }
    std::cout << "dt is" << dt << std::endl;

};

// ---------------------- Constructor ------------------------ //
solver::solver(double x_0, double x_1, double y_0, double y_1, double t0, double t1, int n_Cells, int n_Ghosts, double c, double g)
    : fluid1(n_Cells, n_Ghosts), fluid2(n_Cells,n_Ghosts), gamma(g) {
    assert(t1>t0);
    assert(g > 0);

    x0 = x_0; x1 = x_1;
    y0 = y_0; y1 = y_1;
    startTime = t0; endTime = t1;
    nCells = n_Cells;
    nGhost = n_Ghosts;

    resize2D(nCells+2*nGhost,nCells+2*nGhost,phi);
    resize2D(nCells+2*nGhost,nCells+2*nGhost,phiPlus1);

    dx = (x1 - x0)/nCells;
    dy = (y1 - y0)/nCells;


    slopeWeight = 1;

    cour = c;

    variables = {"rho","vx","vy","p"};

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

double solver::get_dy()const{
    return dy;
};

void solver::init(std::vector< std::vector< std::array<double,4> > > init, fluid& fluid) {
    assert(init.size() == fluid.u.size());
    fluid.u = init;
};

void solver::phiInit(std::vector< std::vector < double > > init) {
    assert(init.size() == phi.size());
    phi = init;
};

fluid& solver::get_fluid(bool f){
    if (f==0){
        return fluid1;
    } else if (f==1){
        return fluid2;
    } else {
        throw std::invalid_argument("Invalid fluid index");
    }
}

// -------------------- writing functions --------------------- //

void solver::setWriteInterval(double wt){
    writeInterval = wt;
}

void solver::writeData(std::string varName) const{ // modify to write both fluids
        // making data file
        std::string filename = dirName + "/" + varName + "/" + std::to_string(time).substr(0,4);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (std::vector<double>::size_type i=nGhost; i<(fluid1.u.size()-nGhost);i++){
                double y = y0 + (i-nGhost+0.5)*dy;
                for (std::vector<double>::size_type j=nGhost; j<(fluid1.u[0].size()-nGhost);j++){
                    double x = x0 + (j-nGhost+0.5)*dx;
                    writeFile << x << " " << y << " " << eos[0]->consvToPrim(fluid1.u[i][j])[varMap.at(varName)] << " " << eos[1]->consvToPrim(fluid2.u[i][j])[varMap.at(varName)] << " " << phi[i][j] << std::endl;
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

std::array<double,4> solver::set_vals(double rho, double vx, double vy, double pres){
        std::array<double,4> result = {rho,vx,vy,pres};
        return result;
};


// ------------------- Overloading operators ------------------ //

std::array<double, 4> operator-(const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

std::array<double, 4> operator+(const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

std::array<double, 4> operator*(const double scalar,const std::array<double, 4>& lhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = scalar*lhs[i];
    }
    return result;
}

std::array<double, 4> operator*(const std::array<double, 4>& rhs,const std::array<double, 4>& lhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = rhs[i]*lhs[i];
    }
    return result;
}

std::array<double, 4> elementDivide(const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
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

std::array<double, 4> operator/(const std::array<double, 4>& lhs, const double scalar) {
    std::array<double, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = lhs[i]/scalar;
    }
    return result;
}



