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
    int nCells = 20;
    int nGhost = 2;
    double x0{0}, x1{1};
    double y0{0}, y1{1};
    double startTime = 0.0, endTime = 0.05;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // constructing EOS objects
    idealGas idgas1(1.4); idealGas idgas2(1.4);
    std::array<EOS*,2> materials = {&idgas1,&idgas2};

    solver sim(x0,x1,y0,y1,startTime,endTime,nCells,nGhost,cour,1.4);
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
        return sim.minbee(input);
    }; // because fEuler isnt static

    sim.fluxMethod = [&sim](fluid& f, EOS* e) {
        sim.godunov(f, e); // Change this to sim.SLIC or sim.godunov as needed
    };

    std::vector< std::vector<std::array<double,4>>> uInit1, uInit2;
    std::vector< std::vector<double>> phiInit;
    // resizing:
    resize2D(nCells+2*nGhost,nCells+2*nGhost,uInit1); resize2D(nCells+2*nGhost,nCells+2*nGhost,uInit2); // resizes number of rows
    solver::resize2D(nCells+2*nGhost,nCells+2*nGhost,phiInit); 

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;
    for (size_t i=0; i<uInit1.size();i++){
        double y = y0 + (i-static_cast<double>(sim.ghosts())+0.5)*sim.get_dy();
        for (size_t j=0; j<uInit1[0].size();j++){
            double x = x0 + (j-static_cast<double>(sim.ghosts())+0.5)*sim.get_dx();
            
            // RIEMANN PROBLEMS //

            //phiInit[i][j] = (y + x - 1)/sqrt(2);
            phiInit[i][j] = x - 0.5;
                        
            if (x<0.5){ 
                uInit1[i][j] = sim.set_vals(1,-2.0,0,0.4);
                uInit2[i][j] = sim.set_vals(1,-2.0,0,0.4);
            } else {
                uInit1[i][j] = sim.set_vals(1,2.0,0,0.4);
                uInit2[i][j] = sim.set_vals(1,2.0,0,0.4);
            }

            // ----------------- //

            // WANG TEST B
            /*
            if (x <= 0.5){
                phiInit[i][j] = (x-0.4);
            } else {
                phiInit[i][j] = (0.6-x);
            }
            */

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
            
            // SOD TEST //
            /*
            phiInit[i][j] = -(sqrt((x-1)*(x-1)+(y-1)*(y-1))-0.5);
            
            // sod (x-1)*(x-1)+(y-1)*(y-1) < 0.4*0.4
            if ((x-1)*(x-1)+(y-1)*(y-1) < 0.5*0.5){
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(1,0,0,1);
            } else {
                uInit1[i][j] = sim.set_vals(0.125,0,0,0.1);
                uInit2[i][j] = sim.set_vals(0.125,0,0,0.1);
            }
            */
            
            
            /*
            // WANG B //
            if (x <= 0.25){
                uInit1[i][j] = sim.set_vals(1.3765,0.3948,0,1.57);
                uInit2[i][j] = sim.set_vals(1.3765,0.3948,0,1.57);
            } else if (x <= 0.4 && x > 0.25){
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(1,0,0,1);
            } else if (x <= 0.6 && x > 0.4){
                uInit1[i][j] = sim.set_vals(0.1380,0,0,1);
                uInit2[i][j] = sim.set_vals(0.1380,0,0,1);
            } else {
                uInit1[i][j] = sim.set_vals(1,0,0,1);
                uInit2[i][j] = sim.set_vals(1,0,0,1);
            }
            */
            
        }

    }
    //std::cout << "Printing uInit" << std::endl;
    //sim.print_arr(uInit1,0);

    //throw std::runtime_error("have a look");


    for (std::size_t var = 0; var < sim.variables.size(); ++var) {
        std::ofstream initFile(dirname + "/" + sim.variables[var] + "/0");
        for (std::vector<double>::size_type i = sim.ghosts(); i < uInit1.size() - sim.ghosts(); i++) {
            double x = x0 + (i - sim.ghosts() + 0.5) * sim.get_dx();
            for (std::vector<double>::size_type j=sim.ghosts(); j<uInit1.size()-sim.ghosts();j++){
                double y = y0 + (j-sim.ghosts() + 0.5) * sim.get_dy();
                initFile << x << " " << y << " " << sim.eos[0]->consvToPrim(uInit1[i][j])[var] << " " << sim.eos[1]->consvToPrim(uInit2[i][j])[var] << " " << phiInit[i][j] << std::endl;
            }
            initFile << "\n";
        }
        initFile.close();
    }

    // Initialising u fields: //

    for (std::vector<double>::size_type i=sim.ghosts(); i<uInit1.size()-sim.ghosts();i++){
        for (std::vector<double>::size_type j=sim.ghosts(); j<uInit1.size()-sim.ghosts();j++){
            uInit1[i][j] = sim.eos[0]->primToConsv(uInit1[i][j]);
            uInit2[i][j] = sim.eos[1]->primToConsv(uInit2[i][j]);
        }
    }

    sim.init(uInit1,sim.get_fluid(0)); sim.init(uInit2,sim.get_fluid(1));
    sim.phiInit(phiInit);


    // Calculating exact result with exact riemann solver

    std::array<double,4> sL = {1,-2,0,0.4};//{1000,0,1e9};
    double gamma1 = 1.4;

    std::array<double,4> sR = {1,2,0,0.4}; // {50,0,1e5};
    double gamma2 = 1.4;

    double time = 0.025;
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

    sim.run();
    std::cout << "run finished" << std::endl;

    return 0;
}