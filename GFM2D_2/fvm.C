#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <array>
#include "GFMSolver.H"


// g++ -o fvm fvm.C solver.C


int main(){
    
    // int nCellsX = 100;
    // int nCellsY = 100;
    // double x0{0}, x1{0.01}, y0{0}, y1{0.01};
    // double startTime = 0.0, endTime = 0.0002;
    
    // int nCellsX = 1300;
    // int nCellsY = 178;
    // double x0{0.0}, x1{0.325}, y0{-0.0445}, y1{0.0445};
    // double startTime = 0.0, endTime = 1200e-6;

    // int nCellsX = 1000;
    // int nCellsY = 1;
    // double x0{0}, x1{1}, y0{-0.5/nCellsX}, y1{0.5/nCellsX};
    // double startTime = 0.0, endTime = 150e-6;

    int nCellsX = 100;
    int nCellsY = 100;
    double x0{0}, x1{1}, y0{0}, y1{1};
    double startTime = 0.0, endTime = 0.25;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // construct EOS object and give to sovler
    //stiffenedGas stiffgas(7.15,3e8,0);
    idealGas idgas1(1.4); idealGas idgas2(1.4);

    std::array<EOS*,2> materials = {&idgas1,&idgas2};
    solver sim(x0,x1,y0,y1,startTime,endTime,nCellsX,nCellsY,2,cour);
    sim.setEOS(materials);

    enum fld{FLUID1=1,FLUID2};

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;

    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    sim.dirName = dirname;


    sim.setWriteInterval(0.01);
    sim.timeMulti = 1;

    sim.setBCs = [&sim](fluid& f) {
        sim.get_boundary().transmissiveTopBC(f);
        sim.get_boundary().transmissiveBottomBC(f);
        sim.get_boundary().transmissiveLeftBC(f);
        sim.get_boundary().transmissiveRightBC(f);
        // 1st bool: left/right = IS reflective, 2nd bool: top/bottom = IS reflective
        sim.get_boundary().updateBottomLeftCorner(f,false,false);
        sim.get_boundary().updateTopRightCorner(f,false,false);
        sim.get_boundary().updateBottomRightCorner(f,false,false);
        sim.get_boundary().updateTopLeftCorner(f,false,false);
    };

    std::cout << "set BCs successfully" << std::endl;

    // Set the flux function
    sim.flux = [&sim](std::array<double,4> input, EOS* eos, bool direction){
        return sim.fEuler(input,eos,direction);
    }; // because fEuler isnt static

    sim.slopeLim = [&sim](std::array<double,4> input){
        return sim.superbee(input);
    }; // because fEuler isnt static

    sim.PRIM = true; // primitive or conserved variable slope limiting
    sim.fluxMethod = [&sim](bool direction, fluid& f, EOS* eos) {
        sim.MUSCL(direction,f,eos); // Change this to other solver if desired
    };

    std::cout << "setting initial conditions" << sim.ghosts() << std::endl;

    // initialising arrays
    std::vector< std::vector<std::array<double,4> > > uInit1; resizeVector(uInit1,nCellsY+2*sim.ghosts(),nCellsX+2*sim.ghosts());
    std::vector< std::vector<std::array<double,4> > > uInit2; resizeVector(uInit2,nCellsY+2*sim.ghosts(),nCellsX+2*sim.ghosts());
    std::vector< std::vector<double>> phiInit; resizeVector(phiInit,nCellsY+2*sim.ghosts(),nCellsX+2*sim.ghosts());

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;
    for (size_t i=sim.ghosts(); i<uInit1.size()-sim.ghosts();++i){
        double y = y0 + (i-sim.ghosts()+0.5)*sim.get_dxdy()[1];
        for (size_t j=sim.ghosts(); j<uInit1[0].size()-sim.ghosts();j++){
            double discPos = 0.5;
            double x = x0 + (j-sim.ghosts()+0.5)*sim.get_dxdy()[0];
            // RIEMANN PROBLEMS //
            
            // toro3 -- 1.0,0.0,0.0,1000 -- 1.0, 0.0, 0.0, 0.01
            // toro5 -- 5.99924,0,19.5975,460.894 -- 5.99242,0,-6.19633,46.0950
            //phiInit[i][j] = (y + x)/sqrt(2);
            //phiInit[i][j] = (y-x)/sqrt(2);

            
            phiInit[i][j] = y-0.5;
                        
            if (y>0.5){
                uInit1[i][j] = sim.set_vals(FLUID1, 1,0.0,0.0,1); // 1,0.0,0.0,1 -- 0.125,0.0,0,0.1 -- toro1
                uInit2[i][j] = sim.set_vals(FLUID2, 1,0.0,0.0,1); // 1,0,-2,0.4 -- 1,0,2,0.4 -- toro2
            } else {
                uInit1[i][j] = sim.set_vals(FLUID1, 0.125,0.0,0,0.1);
                uInit2[i][j] = sim.set_vals(FLUID2, 0.125,0.0,0,0.1);
            }
                
            

            //  -- CYlindrical Sod -- //
            /*
            phiInit[i][j] = (sqrt((x)*(x)+(y)*(y))-0.2);
            
            // sod (x-1)*(x-1)+(y-1)*(y-1) < 0.4*0.4
            if (x*x+y*y < 0.2*0.2){
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(0.125,0.0,0,0.1);
            } else {
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(0.125,0.0,0,0.1);
            }
            */
            
                
            
            // ----------------- //

            // WANG TEST B
            
            //2D:
            /*
            if (x+y < 1){
                phiInit[i][j] = (y + x - (1-0.1*sqrt(2)))/sqrt(2);
                //phiInit[i][j] = (y-0.4);
            } else {
                //phiInit[i][j] = (0.6-y);
                phiInit[i][j] = -(y + x - (1+0.1*sqrt(2)))/sqrt(2);
            }

            
            
            if (x+y <= (1-0.25*sqrt(2))){
                uInit1[i][j] = sim.set_vals(FLUID1, 1.3765,0.3948/sqrt(2),0.3948/sqrt(2),1.57);
                uInit2[i][j] = sim.set_vals(FLUID2,1.3765,0.3948/sqrt(2),0.3948/sqrt(2),1.57);
            } else if (x+y <= (1-sqrt(2)*0.1) && x+y > (1-0.25*sqrt(2))){
                uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1);
                uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1);
            } else if (x+y <= (1+0.1*sqrt(2)) && x+y > (1-sqrt(2)*0.1)){
                uInit1[i][j] = sim.set_vals(FLUID1,0.1380,0,0,1);
                uInit2[i][j] = sim.set_vals(FLUID2,0.1380,0,0,1);
            } else {
                uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1);
                uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1);
            }
                */

            // // // Define the rotation angle in radians (35 degrees)
            // double theta = 90 * M_PI / 180.0;
            // // // Compute the projection along the interface normal (rotated by 35°)
            // double p = (x * cos(theta) + y * sin(theta));

            // // // Set the level-set variable (phi) with a smeared interface centered at p = 1
            // if (p < 0.5) {
            //     phiInit[i][j] = p - (0.5 - 0.1);
            // } else {
            //     phiInit[i][j] = -(p - (0.5 + 0.1));
            // }

            // // Physically realistic variables:

            // // Set the initial velocity fields based on the location relative to the interface
            // if (p <= (0.5 - 0.45)) {
            //     uInit1[i][j] = sim.set_vals(FLUID1, 1.3333, 0.3535*sqrt(1e5) * cos(theta), 0.3535*sqrt(1e5) * sin(theta), 1.5e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2, 1.3333, 0.3535*sqrt(1e5) * cos(theta), 0.3535*sqrt(1e5) * sin(theta), 1.5e5);
            // } else if (p <= (0.5 - 0.1) && p > (0.5 - 0.45)) {
            //     uInit1[i][j] = sim.set_vals(FLUID1, 1, 0, 0, 1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2, 1, 0, 0, 1e5);
            // } else if (p <= (0.5 + 0.1) && p > (0.5 - 0.1)) {
            //     uInit1[i][j] = sim.set_vals(FLUID1, 0.1379,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2, 0.1379,0,0,1e5);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1, 1, 0, 0, 1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2, 1, 0, 0, 1e5);
            // }
                
            
            // 1D:
            
            // if (x<0.5){
            //     phiInit[i][j] = (x-0.4);
            // } else {
            //     phiInit[i][j] = (0.6-x);
            // }
            
            // if (x<=0.05){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.3333,0.3535*sqrt(1e5),0,1.5e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.3333,0.3535*sqrt(1e5),0,1.5e5);
            // } else if (x <= 0.4){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1e5); 
            // } else if (x <= 0.6){
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.1379,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.1379,0,0,1e5);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1e5);
            // }
                
            
            // Mach 20 shock //

            // 1D:
            
            // if (x<0.5){
            //     phiInit[i][j] = (x-0.4);
            // } else {
            //     phiInit[i][j] = (0.6-x);
            // }
            
            // if (x<=0.05){
            //     uInit1[i][j] = sim.set_vals(FLUID1,5.92593,6220.51,0,4.665e7);
            //     uInit2[i][j] = sim.set_vals(FLUID2,5.92593,6220.51,0,4.665e7);
            // } else if (x <= 0.4){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1e5); 
            // } else if (x <= 0.6){
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.1379,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.1379,0,0,1e5);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1e5);
            // }
                
    
            
            
            

            // ----------------- //

            // ---------- Fedkiw B ---------- //
            
            // phiInit[i][j] = x-0.5;
            
            // if (x<=0.05){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.3333,0.3535*sqrt(1e5),0,1.5e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.3333,0.3535*sqrt(1e5),0,1.5e5);
            // } else if (x <= 0.5){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1e5); 
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.1379,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.1379,0,0,1e5);
            // }
                
            
            

            //std::cout << i << " " << j << std::endl;
            
            /*
            if (x > 0.5 && y > 0.5){ 
                uInit[i][j] = sim.set_vals(1,0.75,-0.5,1);
            } else if(x < 0.5 && y > 0.5){
                uInit[i][j] = sim.set_vals(2,0.75,0.5,1);
            } else if (x > 0.5 && y < 0.5){
                uInit[i][j] = sim.set_vals(3,-0.75,-0.5,1);
            } else {
                uInit[i][j] = sim.set_vals(1,-0.75,0.5,1);
            }
            */
            
            //std::cout << phiInit[i][j] << " ";
            
            //  SHOCK-BUBBLE TEST (fedkiw) //
            
            // phiInit[i][j] = -(sqrt((x-0.175)*(x-0.175)+(y)*(y))-0.025);
            
            // if ((x-0.175)*(x-0.175)+(y*y) < 0.025*0.025){
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.1380,0,0,1);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.1380,0,0,1);
            // } else if (x>=0.225){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.3764,-0.394,0,1.5698);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.3764,-0.394,0,1.5698);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,1,0,0,1);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1,0,0,1);
            // }

            // Marquina shock bubble

            // phiInit[i][j] = -(sqrt((x-0.175)*(x-0.175)+(y)*(y))-0.025);
            
            // if ((x-0.175)*(x-0.175)+(y*y) < 0.025*0.025){
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.2228,0,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.2228,0,0,101325);
            // } else if (x>=0.225){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.6861,-156.26,0,250638);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.6861,-156.26,0,250638);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.225,0,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.225,0,0,101325);
            // }

            // Haas-Sturtevant Shock Bubble
            
            // phiInit[i][j] = -(sqrt((x-0.175)*(x-0.175)+(y)*(y))-0.025);
            
            // if ((x-0.175)*(x-0.175)+(y*y) < 0.025*0.025){
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.223,0,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.223,0,0,101325);
            // } else if (x>=0.225){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.686,-113.5,0,159060);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.686,-113.5,0,159060);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.225,0,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.225,0,0,101325);
            // }

            // ------------- Richtmeyer - Meshkov instability -------------- //

            // phiInit[i][j] = x - 0.5 + 0.056*cos(2*M_PI*y);

            // if (x < 0.4){
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.67,161.68,0,1.5e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.67,161.68,0,1.5e5);
            // } else if (x < 0.5 - 0.056*cos(2*M_PI*y) ){ //
            //     uInit1[i][j] = sim.set_vals(FLUID1,0.5,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,0.5,0,0,1e5);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,2.5,0,0,1e5);
            //     uInit2[i][j] = sim.set_vals(FLUID2,2.5,0,0,1e5);
            // }


            // ------------- Water collision ------------------- // 

            // if (x>0.145){
            //     phiInit[i][j] = -(sqrt((x-0.175)*(x-0.175)+(y)*(y))-0.025);
            // } else {
            //     phiInit[i][j] = -(sqrt((x-0.115)*(x-0.115)+(y)*(y))-0.007);
            // }
            
            
            // if ((x-0.175)*(x-0.175)+(y*y) < 0.025*0.025){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1000,0,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1000,0,0,101325);
            // } else if ((x-0.115)*(x-0.115)+(y*y) < 0.007*0.007){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1000,3000,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1000,3000,0,101325);
            // } else if ((x-0.115)*(x-0.115)+(y*y) < 0.008*0.008){
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.225,3000,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.225,3000,0,101325);
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1,1.225,0,0,101325);
            //     uInit2[i][j] = sim.set_vals(FLUID2,1.225,0,0,101325);
            // }

            // ----------- Kelvin-Helmholtz instability ------------ //
            
            // phiInit[i][j] = (sqrt((x-0.005)*(x-0.005)+(y+017)*(y+0.017))-sqrt(0.0005));
            
            // if (sqrt((x-0.005)*(x-0.005)+(y+0.017)*(y+0.017)) > sqrt(0.0005)){
            //     uInit1[i][j] = sim.set_vals(FLUID1, 1.65,0.0,0,158000); // 1,0.0,0.0,1 -- 0.125,0.0,0,0.1 -- toro1
            //     uInit2[i][j] = sim.set_vals(FLUID2, 1.65,0.0,0,158000); // 1,0,-2,0.4 -- 1,0,2,0.4 -- toro2
            // } else {
            //     uInit1[i][j] = sim.set_vals(FLUID1, 0.288,-70,0.0,158000);
            //     uInit2[i][j] = sim.set_vals(FLUID2, 0.288,-70,0.0,158000);
            // }

                
                
                
            
            
            
        }

    }
    

    // Initialising u field: //

    sim.init(uInit1,sim.get_fluid(0)); sim.init(uInit2,sim.get_fluid(1));
    sim.phiInit(phiInit);

    
    // for (int j = 0; j < uInit.size(); j++){
    //     std::cout << uInit[j] << std::endl;
    // }

    sim.runGFM();
    //sim.run();
    std::cout << "run finished" << std::endl;

    // (optional) Calculating exact result with exact riemann solver

    /*
    std::array<double,3> sL = {1.333,0.3535*sqrt(1e5),1.5e5};//{1000,0,1e9};
    double gamma1 = 1.4;

    std::array<double,3> sR = {0.1379,0,1e5}; // {50,0,1e5};
    double gamma2 = 1.67;

    double time = 0.0002;
    double discPosition = 0.5;

    riemann solution(gamma1,gamma2,sL,sR,0,1,discPosition,time,100,0,0);
    solution.dirName = dirname;
    // initialise data locations
    std::vector<std::string> riemVariables = {"rhoExact","vExact","pExact"};
    for (size_t i = 0; i<riemVariables.size(); ++i){
        std::string name = dirname + "/" + riemVariables[i];
        std::filesystem::create_directory(name);
    }
        */

    // Exact Riemann solver
    //double pStar = solution.exctRiemann()[2];
    //solution.exactRiemannSolution(pStar,discPosition);

    // end exact calculation

    return 0;
}