#include "riem.H"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <unordered_map>
#include <functional>
#include "EOS.H"
#include "operators.H"

// -------------------- exact Riemann solver ------------------------ //

double riemann::gamma(bool left){
    if (left == 0){
        return gammaR;
    } else {
        return gammaL;
    }
};

double riemann::pInf(bool left){
    if (left == 0){
        return pInfRight;
    } else {
        return pInfLeft;
    }
};


double riemann::fShock(double pStar, std::array<double,3> s, bool left){
    //std::cout << "Shock" << std::endl;
    double A = 2.0/((gamma(left)+1)*s[0]);
    double B = (s[2]+pInf(left))*(gamma(left)-1)/(gamma(left)+1);
    //correct --> return sqrt((pStar-s[2])/(s[0])*(1-((gamma(left)-1)*(pStar+s[2])+2.0*(s[2]+gamma(left)*pInf(left)))/((gamma(left)-1)*(pStar+s[2])+2.0*(pStar+gamma(left)*pInf(left)))));
    return (pStar-s[2])*sqrt(A/(B+pStar+pInf(left)));
}

double riemann::fRare(double pStar, std::array<double,3> s, bool left){
    //std::cout << "Rarefaction" << std::endl;
    //original --> return (2*sqrt(gamma(left)*(s[2]+pInf(left))/s[0])/(gamma(left)-1))*(pow((pStar+pInf(left))/(s[2]+pInf(left)),(gamma(left)-1)/(2*gamma(left)))-1);
    return (-1)*(2.0*calcSoundSpeed(s,left)/(gamma(left)-1))*(1-pow(((pStar+pInf(left))/(s[2]+pInf(left))),(gamma(left)-1)/(2.0*gamma(left))));
}

double riemann::fPrm(double pStar, std::array<double,3> s, bool left){
    double A = 2.0/((gamma(left)+1)*s[0]);
    double B = (s[2]+pInf(left))*(gamma(left)-1)/(gamma(left)+1);
    if (pStar >= s[2]){ // shock case
        return sqrt(A/(B+pStar+pInf(left)))*(1-0.5*(pStar-s[2])*(1.0/(B+pStar+pInf(left)))); // wrong
        // correct --> return ((pStar - s[2]) * ((gamma(left) + 1) * pStar + (3 * s[2] + 4 * pInf(left)) * gamma(left) - s[2])) / (s[0] * pow(((gamma(left) + 1) * pStar + (s[2] + 2 * pInf(left)) * gamma(left) - s[2]),2) * sqrt(((pStar - s[2]) * (1 - ((gamma(left) - 1) * (pStar + s[2]) + 2 * (pInf(left) * gamma(left) + s[2])) / (2 * (pStar + pInf(left) * gamma(left)) + (gamma(left) - 1) * (pStar + s[2])))) / s[0]));
    } else {
        // original --> return (sqrt((s[2]+pInf(left))/(gamma(left)*s[0])))*pow((pStar+pInf(left)),(-1)*((gamma(left)+1)/(2*gamma(left))))*pow(1/(s[2]+pInf(left)),(gamma(left)-1)/2*gamma(left));
        return (2*calcSoundSpeed(s,left)/(gamma(left)-1))*(pow((1.0/(s[2]+pInf(left))),(gamma(left)-1)/(2.0*gamma(left))))*(pow((pStar+pInf(left)),(-1)*(gamma(left)+1)/(2.0*gamma(left))))*((gamma(left)-1)/(2*gamma(left)));
    }
}

double riemann::fRiemann(double pStar){
    if (pStar >= state1[2]){
        if (pStar >= state2[2]){
            return fShock(pStar,state1,1)+fShock(pStar,state2,0)+(state2[1]-state1[1]);
        }
        else {
            return fShock(pStar,state1,1)+fRare(pStar,state2,0)+(state2[1]-state1[1]);
        }
    }
    else {
        if (pStar >= state2[2]){
            return fRare(pStar,state1,1)+fShock(pStar,state2,0)+(state2[1]-state1[1]);
        }
        else {
            return fRare(pStar,state1,1)+fRare(pStar,state2,0)+(state2[1]-state1[1]);
        }
    }
}

std::array<double,4> riemann::solveRiem(double pStar, std::array<double,3> s, bool left){
    std::array<double,4> result = {0,0,pStar,0};
    double pStarTilde = pStar + pInf(left);
    double pTilde = s[2] + pInf(left);

    if (pStar > s[2]){
        double F = sqrt(s[0]*((pStarTilde+pTilde*((gamma(left)-1)/(gamma(left)+1)))/(2.0/(gamma(left)+1))));
        result[0] = s[0]*((pStarTilde/pTilde)+(gamma(left)-1)/(gamma(left)+1))/(((gamma(left)-1)/(gamma(left)+1))*(pStarTilde/pTilde)+1);
        if (left == 0){
            //result[3] = s[1] + sqrt((gamma(left)*(s[2]+pInf(left))/s[0])*(((gamma(left)+1)/(2*gamma(left)))*(pStar/s[2])+(gamma(left)-1)/(2*gamma(left))));
            result[3] = s[1] + (1.0/s[0])*F;
            result[1] = s[1] + (pStarTilde-pTilde)/F;
        } else if (left == 1) {
            result[3] = s[1] - (1.0/s[0])*F;
            result[1] = s[1] - (pStarTilde-pTilde)/F;
        }
    } else {
        result[0] = s[0]*pow(pStarTilde/pTilde,1/gamma(left));
        if (left == 0){
            //result[1] = s[1] + (2*sqrt((gamma(left)*(s[2]+pInf(left))/s[0]))/(gamma(left)-1))*(pow(pStar/s[2],(gamma(left)-1)/(2*gamma(left)))-1);
            result[1] = s[1]+fRare(pStar,s,left);
        } else if (left == 1){
            //result[1] = s[1] - (2*sqrt((gamma(left)*(s[2]+pInf(left))/s[0]))/(gamma(left)-1))*(pow(pStar/s[2],(gamma(left)-1)/(2*gamma(left)))-1);
            result[1] = s[1]-fRare(pStar,s,left);
        }
    }

    return result;
}

std::array<double,3> riemann::rareFan(std::array<double,3> state, bool left, double x=0){
    std::array<double,3> result = {0,0,0};
    double eps = (x-x_disc)/endTime;
    //std::cout << "Eps = " << eps << std::endl;
    if (left == 1){
        result[1] = (state[1]*(gamma(left)-1)+2.0*(calcSoundSpeed(state,left)+eps))/(gamma(left)+1);
        //std::cout << result[1] << std::endl;
        result[0] = state[0]*pow((result[1]-eps)/calcSoundSpeed(state,left),2.0/(gamma(left)-1));
        result[2] = (state[2]+pInf(left))*pow((result[1]-eps)/calcSoundSpeed(state,left),2.0*gamma(left)/(gamma(left)-1))-pInf(left);
        return result;
    } else {
        result[1] = (state[1]*(gamma(left)-1)+2.0*(eps-calcSoundSpeed(state,left)))/(gamma(left)+1);
        //std::cout << result[1] << std::endl;
        result[0] = state[0]*pow((eps-result[1])/calcSoundSpeed(state,left),2.0/(gamma(left)-1));
        result[2] = (state[2]+pInf(left))*pow((eps-result[1])/calcSoundSpeed(state,left),2.0*gamma(left)/(gamma(left)-1))-pInf(left);
        return result;
    };
}

std::array<double,4> riemann::exctRiemann(){ // states will be in primitive variable form by default

    double pStar = 0.5*(state1[2]+state2[2]);
    double pOld = pStar;
    double epsilon = 1e-8;
    int cnt = 0;
    do {
        cnt++;
        pOld = std::max(1e-8,pStar);
        //std::cout << fRiemann(pOld) << " " << fPrm(pOld, state1,1)+fPrm(pOld,state2,0) << std::endl;
        pStar = pOld - (fRiemann(pOld)/(fPrm(pOld, state1,1)+fPrm(pOld,state2,0)));
        //std::cout << "iteration " << cnt << ", p* = " << pStar << ", pOld = " << pOld << std::endl;

    } while (std::fabs(pStar-pOld)/pOld > epsilon);

    //std::cout << "converged, p* = " << pStar << std::endl; 

    std::array<double,4> leftStar = {0,0,0,0};
    std::array<double,4> rightStar = {0,0,0,0};

    leftStar = solveRiem(pStar, state1, 1);
    rightStar = solveRiem(pStar, state2, 0);

    //std::cout << leftStar[0] << " " << leftStar[1] << " " << leftStar[2] << " " << leftStar[3] << std::endl;
    //std::cout << rightStar[0] << " " << rightStar[1] << " " << rightStar[2] << " " << rightStar[3] << std::endl;


    if (leftStar[3] > 0 || state1[1]-calcSoundSpeed(state1,1) > 0){ // right moving left shock or u_L - cs_L > 0
        return state1orig;
    } else if (rightStar[3] < 0 || state2[1]+calcSoundSpeed(state2,0) < 0){ // left moving right shock or (u_R+cs_R) < 0
        return state2orig;
    } else {
        double cL = calcSoundSpeed(state1,1);
        double cR = calcSoundSpeed(state2,0);
        double cLStar = sqrt((leftStar[2]+pInf(1))*gamma(1)/leftStar[0]);
        double cRStar = sqrt((rightStar[2]+pInf(0))*gamma(0)/rightStar[0]);

        if ((state1[1]-cL)*(leftStar[1]-cLStar) < 0){ // (u_L-c_L)*(u*L-c*L) < 0 - ie different signs; disc. is in left rarefaction
            std::array<double,3> sol = rareFan(state1,1);
            if (direction==0){ //y direction
                return {sol[0],state1orig[1],sol[1],sol[2]};
            } else {
                return {sol[0],sol[1],state1orig[2],sol[2]};
            }
        } else if ((state2[1]+cR)*(rightStar[1]+cRStar) < 0){ // (u_R+c_R)*(u*R+c*R) < 0 - ie different signs; disc. is in right rarefaction
            std::array<double,3> sol = rareFan(state2,0);
            if (direction==0){ //y direction
                return {sol[0],state2orig[1],sol[1],sol[2]};
            } else {
                return {sol[0],sol[1],state2orig[2],sol[2]};
            }
        } else {
            if (leftStar[1] > 0){ // u*L > 0 - right-moving contact disc.
                if (direction==0){
                    return {leftStar[0],state1orig[1],leftStar[1],leftStar[2]};
                } else {
                    return {leftStar[0],leftStar[1],state1orig[2],leftStar[2]};
                }
            } else { // u*L < 0 - left-moving contact disc.
                if (direction==0){
                    return {rightStar[0],state2orig[1],rightStar[1],rightStar[2]};
                } else {
                    std::array<double,4> result = {rightStar[0],rightStar[1],state2orig[2],rightStar[2]};

                    //std::cout << rightStar[0] << " " << rightStar[1] << " " << state2orig[2] << " " << rightStar[2] << std::endl;
                    return result;
                }
            }
        }
    }

}

std::array<double,4> riemann::interfaceRiemann(bool left){ // states will be in primitive variable form by default

    double pStar = 0.5*(state1[2]+state2[2]);
    double pOld = pStar;
    double epsilon = 1e-8;
    int cnt = 0;
    do {
        cnt++;
        pOld = std::max(1e-8,pStar);
        //std::cout << fRiemann(pOld) << " " << fPrm(pOld, state1,1)+fPrm(pOld,state2,0) << std::endl;
        pStar = pOld - (fRiemann(pOld)/(fPrm(pOld, state1,1)+fPrm(pOld,state2,0)));
        //std::cout << "iteration " << cnt << ", p* = " << pStar << ", pOld = " << pOld << std::endl;

    } while (std::fabs(pStar-pOld)/pOld > epsilon);

    //std::cout << "converged, p* = " << pStar << std::endl; 

    std::array<double,4> leftStar = {0,0,0,0};
    std::array<double,4> rightStar = {0,0,0,0};

    leftStar = solveRiem(pStar, state1, 1);
    rightStar = solveRiem(pStar, state2, 0);

    //std::cout << leftStar[0] << " " << leftStar[1] << " " << leftStar[2] << " " << leftStar[3] << std::endl;
    //std::cout << rightStar[0] << " " << rightStar[1] << " " << rightStar[2] << " " << rightStar[3] << std::endl;
    std::array<double,4> result = left ? std::array<double,4>{rightStar[0],rightStar[1],state2orig[2],rightStar[2]} : std::array<double,4>{leftStar[0],leftStar[1],state1orig[2],leftStar[2]};
    return result;
}

std::array<std::string,5> riemann::waveSignature(double pStar){
    std::array<std::string,5> result;
    if (pStar > state1[2] && pStar > state2[2]){
        result = {"NL","SL","CD","SR","NR"};
    } else if (pStar > state1[2]){
        result = {"NL","SL","CD","RR1","RR2"};
    } else if (pStar > state2[2]){
        result = {"RL1","RL2","CD","SR","NR"};
    } else {
        result = {"RL1","RL2","CD","RR1","RR2"};
    };

    return result;
}

std::array<double,5> riemann::wavePositions(std::array<std::string,5> waveSig, double x_disc, std::array<double,4> lSolution, std::array<double,4> rSolution){
    std::cout << "u*L = " << lSolution[1] << "endTime = " << endTime << "xDisc = " << x_disc << std::endl;
    std::unordered_map<std::string, double> calcWavePos = {
        {"NL", x0},
        {"NR", x1},

        {"SL", x_disc + lSolution[3]*endTime}, // uses calculated shock speed
        {"SR", x_disc + rSolution[3]*endTime},

        {"CD", x_disc + lSolution[1]*endTime}, // v*R = v*L

        {"RL1", x_disc + (state1[1]-calcSoundSpeed(state1,1))*endTime},                    // vL - cL
        {"RL2", x_disc + (lSolution[1]-sqrt(gamma(1)*(lSolution[2]+pInf(1))/lSolution[0]))*endTime}, // v*_L - c*_L
        {"RR1", x_disc + (rSolution[1]+sqrt(gamma(0)*(rSolution[2]+pInf(0))/rSolution[0]))*endTime}, // v*_R + c*_R
        {"RR2", x_disc + (state2[1]+calcSoundSpeed(state2,0))*endTime}                     // vR + cR
    };  

    std::array<double,5> result;

    for (size_t i=0; i<result.size(); ++i){
        std::string sig = waveSig[i];
        result[i] = calcWavePos[sig];
    };
    for (double i : result){
        std::cout << i << " ";
    }
    std::cout << std::endl;
    return result;
};

void riemann::exactRiemannSolution(double pStar, double x_disc){

    // calculating star states
    std::array<double,4> lSolution = solveRiem(pStar, state1, 1);
    std::cout << "lSolution is" << lSolution[0] << " " << lSolution[1] << " " << lSolution[2] << " " << lSolution[3] << std::endl;
    std::array<double,4> rSolution = solveRiem(pStar, state2, 0);
    std::cout << "rSolution is" << rSolution[0] << " " << rSolution[1] << " " << rSolution[2] << " " << rSolution[3] << std::endl;
    
    // finding wave signature
    std::array<std::string,5> waveSig = waveSignature(pStar);
    for (size_t i=0; i<waveSig.size();++i){
        std::cout << waveSig[i] << " ";
    };

    wavePos = wavePositions(waveSig, x_disc, lSolution, rSolution);
    std::cout << "calculated wavepos" << std::endl;
    for (size_t i=0; i<wavePos.size();++i){
        std::cout << wavePos[i] << " ";
    };
    
    std::vector< std::array<double,3> > uRiem;
    uRiem.resize(nCells);

    for (std::vector<double>::size_type i=0; i<uRiem.size();i++){
        double x = x0 + (i+0.5)*dx;
        if (x < wavePos[0]){
            uRiem[i] = state1;
        } else if (x < wavePos[1]){
            if (waveSig[1] == "SL"){
                uRiem[i] = state1;
                
            } else if (waveSig[1] == "RL2"){
                uRiem[i] = rareFan(state1, 1, x);
            }
        } else if (x < wavePos[2]){
            uRiem[i] = {lSolution[0],lSolution[1],lSolution[2]};
        } else if (x < wavePos[3]){
            uRiem[i] = {rSolution[0],rSolution[1],rSolution[2]};
        } else if (x < wavePos[4]){
            if (waveSig[4] == "NR"){
                uRiem[i] = state2;
            } else if (waveSig[4] == "RR2"){
                uRiem[i] = rareFan(state2, 0, x); //marker
            };
        } else if (x >= wavePos[4]){
            uRiem[i] = state2;
        } else {
            throw std::runtime_error("x or bounds not valid");
        };
    }
    std::cout << "printing dp #0" << std::endl;
    std::cout << uRiem[0][0] << " " << uRiem[0][1] << " " << uRiem[0][2] << " " << std::endl; 

    writeData("rhoExact",uRiem); writeData("vExact",uRiem); writeData("pExact",uRiem);
};

void riemann::writeData(std::string varName, std::vector<std::array<double,3>> data) const{
        // making data file
        std::string filename = dirName + "/" + varName + "/" + std::to_string(endTime).substr(0,4);
        std::cout << filename << std::endl;
        std::ofstream writeFile(filename);
        assert(writeFile.is_open());
        if (writeFile.is_open()){
            // writing data
            for (std::vector<double>::size_type i=0; i<data.size();i++){
                double x = x0 + (i+0.5)*dx;

                writeFile << x << " " << data[i][varMap.at(varName)] << std::endl; //consvToPrim(data[i],x<=wavePos[2])
                //std::cout << u[i] << " ";
            }
            writeFile.close();
            std::cout << filename << "written successfully" << std::endl;
        }else{
            std::cout << "failed to write." << std::endl;
        }

}

// Converters


double riemann::calcSoundSpeed(std::array<double,3> arr, bool left){ 
    return std::sqrt(gamma(left)*(arr[2]+pInf(left))/arr[0]);
    // takes the array in primitive variable form!!
}


// Constructor

riemann::riemann(double gL, double gR, std::array<double,4> stateL, std::array<double,4> stateR, bool dir, double x_0,double x_1, double xDisc, double t, int N, double pIL, double pIR)
    : gammaL(gL), gammaR(gR),
        state1orig(stateL), state2orig(stateR),
        direction(dir),
        x0(x_0), x1(x_1),
        x_disc(xDisc), endTime(t),
        nCells(N),
        pInfLeft(pIL), pInfRight(pIR)
{
    if (direction == 0){ // Y direction
        state1 = {stateL[0],stateL[2],stateL[3]};
        state2 = {stateR[0],stateR[2],stateR[3]};
    } else if (direction == 1){
        state1 = {stateL[0],stateL[1],stateL[3]};
        state2 = {stateR[0],stateR[1],stateR[3]};
    }

    dx = (x1-x0)/nCells;

    std::array<std::string,6> variables = {"rho","v","p","rhoExact","vExact","pExact"};

    //varMap["rho"] = 0; varMap["v"] = 1; varMap["p"] = 2;
    for (size_t i = 0; i<variables.size(); ++i){
        varMap[variables[i]] = i%3;
    }
}


