#include "solver.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>


void solver::run(){
    std::cout << "pre-bc" << std::endl;

    transmissiveBC();

    int splitFlip = 0;
    do{
        std::cout << "time is " << time << std::endl;
        this->setDt();

        std::cout << "dt is" << dt << std::endl;
        std::cout << "time is " << time << std::endl;

        transmissiveBC();

        splitFlip++;
        if (splitFlip%2 == 0){
            direction = XDIR;
            this->MUSCL();
            this->pointsUpdate();
            direction = YDIR;
            this->MUSCL();
            this->pointsUpdate();
        } else {
            direction = YDIR;
            this->MUSCL();
            this->pointsUpdate();
            direction = XDIR;
            this->MUSCL();
            this->pointsUpdate();
        };
        std::cout << "points updated" << std::endl;

        transmissiveBC();

        //this->sourceUpdate();

        if (checkWrite==1){
            std::cout << "writing " << time << std::endl;
            writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
        }

    } while (time < endTime);
    std::cout << "simulation finished" << std::endl;
    writeData("rho"); writeData("vx"); writeData("vy"); writeData("p");
};

void solver::SLIC(){

    calcHalfSlopes(); // -------------------------------------------- WIP

    calcr();

    calcUBars();

    updUBars(); //-----------!!!!!!!!!!!!!!!!

    if (direction == XDIR){
        for (size_t i = 0; i < fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluxesX[0].size(); ++j) {
                fluxesX[i][j] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(uBarRX[i][j],uBarLX[i][j+1])),eos[0])+LF(uBarRX[i][j],uBarLX[i][j+1])); // fluxes stored in conserved variable form
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluxesY.size(); ++i) { // why am i getting a Y flux of vx? 
            for (size_t j = 0; j < fluxesY[0].size(); ++j) {
                fluxesY[i][j] = 0.5*(flux(eos[0]->consvToPrim(laxFhalf(uBarRY[i][j],uBarLY[i+1][j])),eos[0])+LF(uBarRY[i][j],uBarLY[i+1][j])); // fluxes stored in conserved variable form
            }
        }
    }
}

std::array<double,4> solver::riemannSolver(std::array<double,4> left,std::array<double,4> right){
    riemann solution(eos[0]->get_gamma(),eos[0]->get_gamma(),eos[0]->consvToPrim(left),eos[0]->consvToPrim(right),direction,0,1,0.5,1,2,0,0);
    return solution.exctRiemann();
}

void solver::MUSCL(){

    calcHalfSlopes(); // -------------------------------------------- WIP

    calcr();

    calcUBars();

    updUBars();

    if (direction == XDIR){
        for (size_t i = 0; i < fluxesX.size(); ++i) {
            for (size_t j = 0; j < fluxesX[0].size(); ++j) {
                fluxesX[i][j] = flux(riemannSolver(uBarRX[i][j],uBarLX[i][j+1]),eos[0]);
            }
        }
    } else if (direction == YDIR){
        for (size_t i = 0; i < fluxesY.size(); ++i) { // why am i getting a Y flux of vx? 
            for (size_t j = 0; j < fluxesY[0].size(); ++j) {
                fluxesY[i][j] = flux(riemannSolver(uBarRY[i][j],uBarLY[i+1][j]),eos[0]);
            }
        }
    }
}



void solver::transmissiveBC(){
    for (int i = 0; i < nGhost; ++i){
        for (size_t j=0; j<uPlus1.size(); ++j){ // sets rows (i)
            uPlus1[i][j] = u[nGhost][j];
            u[i][j] = u[nGhost][j];
            uPlus1[nCells+2*nGhost-1-i][j] = u[nCells+nGhost-1][j];
            u[nCells+2*nGhost-1-i][j] = u[nCells+nGhost-1][j];
        }
    }
    for (size_t i=0; i<uPlus1.size(); ++i){
        for (int j = 0; j < nGhost; ++j){
            uPlus1[i][j] = u[i][nGhost];
            u[i][j] = u[i][nGhost];
            uPlus1[i][nCells+2*nGhost-1-j] = u[i][nCells+nGhost-1];
            u[i][nCells+2*nGhost-1-j] = u[i][nCells+nGhost-1];
        }
    }
}

void solver::cylTransmissiveBC(){
    for (int var =0; var<4;++var){
        if (var == 0 || var == 2 || var == 3){
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<uPlus1.size(); ++j){ // sets rows (ie z)
                    uPlus1[i][j][var] = u[nGhost][j][var];
                    u[i][j][var] = u[nGhost][j][var];
                    uPlus1[nCells+2*nGhost-1-i][j][var] = u[nCells+nGhost-1][j][var];
                    u[nCells+2*nGhost-1-i][j][var] = u[nCells+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    uPlus1[i][nGhost-1-j][var] = u[i][nGhost+j][var];
                    u[i][nGhost-1-j][var] = u[i][nGhost+j][var];
                    uPlus1[i][nCells+2*nGhost-1-j][var] = u[i][nCells+nGhost-1][var];
                    u[i][nCells+2*nGhost-1-j][var] = u[i][nCells+nGhost-1][var];
                }
            }
        } else {
            for (int i = 0; i < nGhost; ++i){
                for (size_t j=0; j<uPlus1.size(); ++j){ // sets rows (ie z)
                    uPlus1[i][j][var] = u[nGhost][j][var];
                    u[i][j][var] = u[nGhost][j][var];
                    uPlus1[nCells+2*nGhost-1-i][j][var] = u[nCells+nGhost-1][j][var];
                    u[nCells+2*nGhost-1-i][j][var] = u[nCells+nGhost-1][j][var];
                }
            }
            for (size_t i=0; i<uPlus1.size(); ++i){ // sets columns (ie r)
                for (int j = 0; j < nGhost; ++j){
                    uPlus1[i][nGhost-1-j][var] = (-1)*u[i][nGhost+j][var];
                    u[i][nGhost-1-j][var] = (-1)*u[i][nGhost+j][var];
                    uPlus1[i][nCells+2*nGhost-1-j][var] = u[i][nCells+nGhost-1][var];
                    u[i][nCells+2*nGhost-1-j][var] = u[i][nCells+nGhost-1][var];
                }
            }
        }
    }
}

void solver::pointsUpdate(){ //second index = x, first index = y // DONE
    //cylTransmissiveBC();
    if (direction == XDIR){
        for (std::vector<double>::size_type i = nGhost-1; i < u.size()-nGhost+1; ++i){
            for (std::vector<double>::size_type j = nGhost; j < u[0].size() - nGhost; ++j) {
                uPlus1[i][j] = u[i][j] - (dt / dx) * (fluxesX[i-nGhost+1][j-nGhost+1] - fluxesX[i-nGhost+1][j-nGhost]); // flux[i + 1] and flux[i] for the update

            }
        }
        u = uPlus1;
        //print_arr(u,RHO);
        std::cout << "X updated" << std::endl;
    } else if (direction == YDIR){
        for (std::vector<double>::size_type i = nGhost; i < u.size()-nGhost; ++i){ // this is the part that causes vx to grow
            for (std::vector<double>::size_type j = nGhost-1; j < u[0].size() - nGhost+1; ++j) {
                uPlus1[i][j] = u[i][j] - (dt / dy) * (fluxesY[i-nGhost+1][j-nGhost+1] - fluxesY[i-nGhost][j-nGhost+1]); // flux[i + 1] and flux[i] for the update
            }
        }
        u = uPlus1;
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
        throw std::runtime_error("laxfhalf doesnt know which dimension!");
    }
};



// -------------- Slope limiting ---------------------------- //

void solver::calcHalfSlopes(){ // only using y-gradient - ie correspondng to y-faces; need extra column for x direction //DONE
    //std::cout << "calculating halfslopes n printing u" << std::endl;
    ////print_arr(u,RHO);
    if (direction == XDIR){
        for (size_t i = 0; i < halfSlopesX.size(); ++i){
            for (size_t j = 0; j < halfSlopesX[0].size(); ++j){
                halfSlopesX[i][j] = u[i+nGhost-1][j+nGhost-1]-u[i+nGhost-1][j+nGhost-2];
            }
        }
    } else {
        for (size_t i = 0; i < halfSlopesY.size(); ++i){
            for (size_t j = 0; j < halfSlopesY[0].size(); ++j){
                halfSlopesY[i][j] = u[i+nGhost-1][j+nGhost-1]-u[i+nGhost-2][j+nGhost-1];
            }
        }
    };
};

void solver::calcr(){ // DONE
    if (direction == XDIR){
        for (size_t i = 0; i < rX.size(); ++i){
            for (size_t j = 0; j < rX[0].size(); ++j){
                rX[i][j] = elementDivide(halfSlopesX[i][j],halfSlopesX[i][j+1]);
            }
        }
    } else {
        for (size_t i = 0; i < rY.size(); ++i){
            for (size_t j = 0; j < rY[0].size(); ++j){
                rY[i][j] = elementDivide(halfSlopesY[i][j],halfSlopesY[i+1][j]);
            }
        }
    };
    ////print_arr(r);
    //std::cout << "calcd r" << std::endl;
};


void solver::calcUBars(){ // calculates for all u except leftmost cell
    if (direction == XDIR){
        for (size_t i = 0; i<uBarLX.size(); ++i){
            for (size_t j = 0; j<uBarLX[0].size(); ++j){
            
                //print_vect(slopeLim(r[i-1]));
                //std::cout << std::endl << "slope ";
                //print_vect(calcSlope(slopeWeight, halfSlopes[i-1], halfSlopes[i]));
                //std::cout << std::endl;
                uBarLX[i][j] = u[i+nGhost-1][j+nGhost-1] - 0.5 * ( slopeLim(rX[i][j]) * calcSlope(slopeWeight, halfSlopesX[i][j], halfSlopesX[i][j+1]) );
                uBarRX[i][j] = u[i+nGhost-1][j+nGhost-1] + 0.5 * ( slopeLim(rX[i][j]) * calcSlope(slopeWeight, halfSlopesX[i][j], halfSlopesX[i][j+1]) );
                //std::cout << slopeLim(rX[i][j])[RHO] << " ";
            }
            //std::cout << std::endl;
        }
    } else {
        for (size_t i = 0; i<(uBarLY.size()); ++i){
            for (size_t j = 0; j<(uBarLY[0].size()); ++j){
            
                //print_vect(slopeLim(r[i-1]));
                //std::cout << std::endl << "slope ";
                //print_vect(calcSlope(slopeWeight, halfSlopes[i-1], halfSlopes[i]));
                //std::cout << std::endl;
                uBarLY[i][j] = u[i+nGhost-1][j+nGhost-1] - 0.5 * ( slopeLim(rY[i][j]) * calcSlope(slopeWeight, halfSlopesY[i][j], halfSlopesY[i+1][j]) );
                uBarRY[i][j] = u[i+nGhost-1][j+nGhost-1] + 0.5 * ( slopeLim(rY[i][j]) * calcSlope(slopeWeight, halfSlopesY[i][j], halfSlopesY[i+1][j]) );
                //std::cout << slopeLim(rY[i][j])[RHO] << " ";
            }
            //std::cout << std::endl;
        }
    };
};


void solver::updUBars(){
    if (direction == XDIR){
        for (size_t i = 0; i<(uBarLupdX.size()); ++i){
            for (size_t j = 0; j<(uBarLupdX[0].size()); ++j){
                uBarLupdX[i][j] = uBarLX[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(uBarRX[i][j]),eos[0])-flux(eos[0]->consvToPrim(uBarLX[i][j]),eos[0]));
                uBarRupdX[i][j] = uBarRX[i][j]-0.5*(dt/dx)*(flux(eos[0]->consvToPrim(uBarRX[i][j]),eos[0])-flux(eos[0]->consvToPrim(uBarLX[i][j]),eos[0]));
            }
        }
        uBarLX = uBarLupdX;
        uBarRX = uBarRupdX;
    } else if (direction == YDIR){
        for (size_t i = 0; i<(uBarLupdY.size()); ++i){
            for (size_t j = 0; j<(uBarLupdY[0].size()); ++j){
                uBarLupdY[i][j] = uBarLY[i][j]-0.5*(dt/dy)*(flux(eos[0]->consvToPrim(uBarRY[i][j]),eos[0])-flux(eos[0]->consvToPrim(uBarLY[i][j]),eos[0]));
                uBarRupdY[i][j] = uBarRY[i][j]-0.5*(dt/dy)*(flux(eos[0]->consvToPrim(uBarRY[i][j]),eos[0])-flux(eos[0]->consvToPrim(uBarLY[i][j]),eos[0]));
            }
        }
        uBarLY = uBarLupdY;
        uBarRY = uBarRupdY;
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

void solver::setDt(){
    double aMax = -1e14;
    std::cout << "setting dt" << std::endl;
    for (std::vector<double>::size_type i = 0; i < u.size(); ++i) {
        for (std::vector<double>::size_type j = 0; j < u.size(); ++j) {
            aMax = std::max(aMax, std::fabs(std::sqrt((eos[0]->consvToPrim(u[i][j])[UX])*(eos[0]->consvToPrim(u[i][j])[UX])+(eos[0]->consvToPrim(u[i][j])[UY])*(eos[0]->consvToPrim(u[i][j])[UY])))+eos[0]->calcSoundSpeed(eos[0]->consvToPrim(u[i][j])));
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
    dx = (x1 - x0)/nCells;
    dy = (y1 - y0)/nCells;
    resize2D(nCells+2*nGhost,nCells+2*nGhost,u);
    resize2D(nCells+2*nGhost,nCells+2*nGhost,uPlus1);
    resize2D(nCells+2,nCells+1,fluxesX);
    resize2D(nCells+1,nCells+2,fluxesY);
    resize2D(nCells+2,nCells+3,halfSlopesX);
    resize2D(nCells+3,nCells+2,halfSlopesY);

    resize2D(nCells+2,nCells+2,rX); // r can only be defined for inner cells
    resize2D(nCells+2,nCells+2,rY); // r can only be defined for inner cells

    resize2D(nCells+2,nCells+2,uBarLX);
    resize2D(nCells+2,nCells+2,uBarLY); //"L" ie left in the y direction corresponds to upper; R to lower.
    //std::cout << uBarL.size();
    resize2D(nCells+2,nCells+2,uBarRX);
    resize2D(nCells+2,nCells+2,uBarRY);

    resize2D(nCells+2,nCells+2,uBarLupdX);
    resize2D(nCells+2,nCells+2,uBarLupdY);
    resize2D(nCells+2,nCells+2,uBarRupdX);
    resize2D(nCells+2,nCells+2,uBarRupdY);

    resize2D(nCells,nCells,sourceResult);

    /*
    resize2D(nCells+2*nGhost-3,nCells+2*nGhost-2,uBarLupdX);
    resize2D(nCells+2*nGhost-2,nCells+2*nGhost-3,uBarLupdY);
    resize2D(nCells+2*nGhost-3,nCells+2*nGhost-2,uBarRupdX);
    resize2D(nCells+2*nGhost-2,nCells+2*nGhost-3,uBarRupdY);
    */

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

void solver::init(std::vector< std::vector< std::array<double,4> > > init) {
    assert(init.size() == u.size());
    u = init;
};

// -------------------- writing functions --------------------- //

void solver::setWriteInterval(double wt){
    writeInterval = wt;
}

void solver::writeData(std::string varName) const{
        // making data file
        std::string filename = dirName + "/" + varName + "/" + std::to_string(time).substr(0,4);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (std::vector<double>::size_type i=nGhost; i<(u.size()-nGhost);i++){
                double y = y0 + (i-nGhost+0.5)*dy;
                for (std::vector<double>::size_type j=nGhost; j<(u[0].size()-nGhost);j++){
                    double x = x0 + (j-nGhost+0.5)*dx;
                    writeFile << x << " " << y << " " << eos[0]->consvToPrim(u[i][j])[varMap.at(varName)] << std::endl;
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



