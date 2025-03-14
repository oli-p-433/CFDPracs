#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <array>
#include "solver.H"

#include <thread>
#include <chrono>

// g++ -o fvm fvm.C solver.C


int main(){
    // int nCellsX = 325; // 325
    // int nCellsY = 90; // 90
    // int nGhost = 2;
    // double x0{0}, x1{0.325}; // 0.325
    // double y0{-0.045}, y1{0.045}; // 0.09
    // double startTime = 0.0, endTime = 0.3;

    int nCellsX = 100; // 325
    int nCellsY = 100; // 90
    int nGhost = 2;
    double x0{-0.5}, x1{0.5}; // 0.325
    double y0{-0.5}, y1{0.5}; // 0.09
    double startTime = 0.0, endTime = 0.25;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // constructing EOS objects
    idealGas idgas1(1.4); idealGas idgas2(1.4);
    std::array<EOS*,2> materials = {&idgas1,&idgas2};

    solver sim(x0,x1,y0,y1,startTime,endTime,nCellsX,nCellsY,nGhost,cour,1.4);
    sim.setEOS(materials);
    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;

    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    sim.dirName = dirname;

    // initialise data locations
    //sim.variables = {"rho","vx","vy","p"};
    std::cout << "Writing directories" << sim.variables.size() << std::endl;
    for (size_t i = 0; i<sim.variables.size(); ++i){
        std::string name = dirname + "/" + sim.variables[i];
        std::cout << name << std::endl;
        std::filesystem::create_directory(name);
    }

    sim.setWriteInterval(0.01);

    // Set the flux function
    sim.flux = [&sim](std::array<double,4> input, EOS* eos){
        return sim.fEuler(input,eos);
    }; // because fEuler isnt static

    sim.slopeLim = [&sim](std::array<double,4> input){
        return sim.vanLeer(input);
    }; // because fEuler isnt static

    sim.fluxMethod = [&sim](fluid& f, EOS* e) {
        sim.MUSCL(f, e); // Change this to sim.SLIC or sim.godunov as needed
    };

    sim.fluxRiemannSolver = [&sim](std::array<double,4> left,std::array<double,4> right, EOS* eos) {
        return sim.HLLC(left,right,eos);
    };

    sim.setBCs = [&sim](fluid& f) {
        sim.get_boundary().transmissiveTopBC(f);
        sim.get_boundary().transmissiveBottomBC(f);
        sim.get_boundary().transmissiveLeftBC(f);
        sim.get_boundary().transmissiveRightBC(f);
        // 1st bool: left/right reflective, 2nd bool: top/bottom reflective
        sim.get_boundary().updateBottomLeftCorner(f,false,false);
        sim.get_boundary().updateTopRightCorner(f,false,false);
        sim.get_boundary().updateBottomRightCorner(f,false,false);
        sim.get_boundary().updateTopLeftCorner(f,false,false);
    };

    std::vector< std::vector<std::array<double,4>>> uInit1, uInit2;
    std::vector< std::vector<double>> phiInit;
    // resizing:
    resize2D(nCellsY+2*nGhost,nCellsX+2*nGhost,uInit1); resize2D(nCellsY+2*nGhost,nCellsX+2*nGhost,uInit2); // resizes number of rows
    solver::resize2D(nCellsX+2*nGhost,nCellsY+2*nGhost,phiInit); 

    std::cout << "setting initial conditions=" << sim.ghosts() << std::endl;
    std::cout << uInit1.size() << uInit1[0].size() << std::endl;
    for (size_t i=0; i<uInit1.size();i++){
        double y = y0 + (i-static_cast<double>(sim.ghosts())+0.5)*sim.get_dy();
        for (size_t j=0; j<uInit1[0].size();j++){
            double x = x0 + (j-static_cast<double>(sim.ghosts())+0.5)*sim.get_dx();
            //std::cout << x << " " << y << std::endl;
            
            // RIEMANN PROBLEMS //
            
            // toro3 -- 1.0,0.0,0.0,1000 -- 1.0, 0.0, 0.0, 0.01
            // toro5 -- 5.99924,0,19.5975,460.894 -- 5.99242,0,-6.19633,46.0950
            //phiInit[i][j] = (y + x)/sqrt(2);
            //phiInit[i][j] = (y-x)/sqrt(2);
            
            phiInit[i][j] = y;
                        
            if (y<0){
                uInit1[i][j] = sim.set_vals(1,0.0,0.0,1); // 1,0.0,0.0,1 -- 0.125,0.0,0,0.1 -- toro1
                uInit2[i][j] = sim.set_vals(1,0.0,0.0,1); // 1,0,-2,0.4 -- 1,0,2,0.4 -- toro2
            } else {
                uInit1[i][j] = sim.set_vals(0.125,0.0,0,0.1);
                uInit2[i][j] = sim.set_vals(0.125,0.0,0,0.1);
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
            /*
            if (x+y < 1){
                phiInit[i][j] = (y + x - (1-0.1*sqrt(2)))/sqrt(2);
                //phiInit[i][j] = (y-0.4);
            } else {
                //phiInit[i][j] = (0.6-y);
                phiInit[i][j] = -(y + x - (1+0.1*sqrt(2)))/sqrt(2);
            }

            
            
            if (x+y <= (1-0.25*sqrt(2))){
                uInit1[i][j] = sim.set_vals(1.3765,0.3948/sqrt(2),0.3948/sqrt(2),1.57);
                uInit2[i][j] = sim.set_vals(1.3765,0.3948/sqrt(2),0.3948/sqrt(2),1.57);
            } else if (x+y <= (1-sqrt(2)*0.1) && x+y > (1-0.25*sqrt(2))){
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(1,0,0,1);
            } else if (x+y <= (1+0.1*sqrt(2)) && x+y > (1-sqrt(2)*0.1)){
                uInit1[i][j] = sim.set_vals(0.1380,0,0,1);
                uInit2[i][j] = sim.set_vals(0.1380,0,0,1);
            } else {
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(1,0,0,1);
            }
            */
            
            /*
            if (y <= 0.25){
                uInit1[i][j] = sim.set_vals(1.3765,0.3948/sqrt(2),0.3948/sqrt(2),1.57);
                uInit2[i][j] = sim.set_vals(1.3765,0.3948/sqrt(2),0.3948/sqrt(2),1.57);
            } else if (y <= 0.4 && y > 0.25){
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(1,0,0,1);
            } else if (y <= 0.6 && y > 0.4){
                uInit1[i][j] = sim.set_vals(0.1380,0,0,1);
                uInit2[i][j] = sim.set_vals(0.1380,0,0,1);
            } else {
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(1,0,0,1);
            }
            */
            
            

            // ----------------- //
            
            

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
            
            // SHOCK-BUBBLE TEST //
            
            // phiInit[i][j] = -(sqrt((x-0.15)*(x-0.15)+(y)*(y))-0.025);
            
            // if ((x-0.15)*(x-0.15)+(y*y) < 0.025*0.025){
            //     uInit1[i][j] = sim.set_vals(0.1380,0,0,1);
            //     uInit2[i][j] = sim.set_vals(0.1380,0,0,1);
            // } else if (x<=0.1){
            //     uInit1[i][j] = sim.set_vals(1.3765,0.3948,0,1.57);
            //     uInit2[i][j] = sim.set_vals(1.3765,0.3948,0,1.57);
            // } else {
            //     uInit1[i][j] = sim.set_vals(1,0,0,1);
            //     uInit2[i][j] = sim.set_vals(1,0,0,1);
            // }
            
                
            
            
            
        }

    }
    //std::cout << "Printing uInit" << std::endl;
    //sim.print_arr(uInit1,0);

    //throw std::runtime_error("have a look");
    std::cout << "set initial conditions" << std::endl;


    for (std::size_t var = 0; var < sim.variables.size(); ++var) {
        std::ofstream initFile(dirname + "/" + sim.variables[var] + "/0");
        for (std::vector<double>::size_type i = 0; i < uInit1.size(); i++) {
            double y = y0 + (i - static_cast<double>(sim.ghosts()) + 0.5) * sim.get_dy();
            for (std::vector<double>::size_type j=0; j<uInit1[0].size();j++){
                double x = x0 + (j-static_cast<double>(sim.ghosts()) + 0.5) * sim.get_dx();
                initFile << x << " " << y << " " << uInit1[i][j][var] << " " << uInit2[i][j][var] << " " << phiInit[i][j] << std::endl;
            }
            initFile << "\n";
        }
        initFile.close();
    }

    // Initialising u fields: //

    for (std::vector<double>::size_type i=0; i<uInit1.size();i++){
        for (std::vector<double>::size_type j=0; j<uInit1[0].size();j++){
            uInit1[i][j] = sim.eos[0]->primToConsv(uInit1[i][j]);
            uInit2[i][j] = sim.eos[1]->primToConsv(uInit2[i][j]);
        }
    }

    sim.init(uInit1,sim.get_fluid(0)); sim.init(uInit2,sim.get_fluid(1));
    sim.phiInit(phiInit);


    // Calculating exact result with exact riemann solver

    std::array<double,4> sL = {1,0,0,1};//{1000,0,1e9};
    double gamma1 = 1.4;

    std::array<double,4> sR = {0.125,0,0,0.1}; // {50,0,1e5};
    double gamma2 = 1.4;

    double time = 0.15;
    double discPosition = 0.5;

    riemann solution(gamma1,gamma2,sL,sR,1,0,1,discPosition,time,100,0,0);
    solution.dirName = dirname;
    // initialise data locations
    std::vector<std::string> riemVariables = {"rhoExact","vExact","pExact"};
    for (size_t i = 0; i<riemVariables.size(); ++i){
        std::string name = dirname + "/" + riemVariables[i];
        std::filesystem::create_directory(name);
    }

    // Exact Riemann solver
    
    std::array<double,4> sol = solution.exctRiemann();
    sim.print_state(sol);
    solution.exactRiemannSolution(sol[3],discPosition);
    
    // end exact calculation

    std::cout << "beginning simulation" << std::endl;
    sim.run();
    std::cout << "run finished" << std::endl;

    return 0;
}