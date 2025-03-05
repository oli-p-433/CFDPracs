#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <array>
#include "solver1D.H"


// g++ -o fvm fvm.C solver.C


int main(){
    int nCells = 200;
    double x0{0}, x1{1};
    double startTime = 0.0, endTime = 0.25;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // construct EOS object and give to sovler
    //stiffenedGas stiffgas(7.15,3e8,0);
    idealGas idgas1(1.4,1.4); idealGas idgas2(1.67,1.67);
    std::array<EOS*,2> materials = {&idgas1,&idgas2};
    solver sim(x0,x1,startTime,endTime,nCells,2,cour);
    sim.setEOS(materials);

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;

    std::string dirname;
    std::cout << std::endl << "Enter data output folder name:"; std::cin >> dirname;  
    std::filesystem::create_directory(dirname);
    sim.dirName = dirname;

    // initialise data locations
    std::vector<std::string> variables = {"rho","v","p"};
    for (size_t i = 0; i<variables.size(); ++i){
        std::string name = dirname + "/" + variables[i];
        std::filesystem::create_directory(name);
    }

    sim.setWriteInterval(0.01);

    // Set the flux function
    sim.flux = [&sim](std::array<double,5> input, EOS* eos){
        return sim.fEuler(input,eos);
    }; // because fEuler isnt static

    sim.slopeLim = [&sim](std::array<double,5> input){
        return sim.minbee(input);
    }; // because fEuler isnt static

    std::vector<std::array<double,5>> uInit; uInit.resize(sim.uPlus1.size());

    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;
    for (std::vector<double>::size_type i=sim.ghosts(); i<uInit.size()-sim.ghosts();i++){
        double discPos = 0.5;
        double x = x0 + (i-sim.ghosts()+0.5)*sim.get_dx();
        /*
        if (x <= 0.05){
            uInit[i] = sim.eos[0]->primToConsv({1,1.3333,0.3535*sqrt(1e5),1.5e5});
        } else if (x <= discPos && x > 0.05){
            uInit[i] = sim.eos[0]->primToConsv({1,1,0,1e5});
        } else {
            uInit[i] = sim.eos[0]->primToConsv({1,1,0,1});
        }
        */
        double vf = (1-(1e-3));
        double rho1{1}, rho2{1};
        if (x < 0.5){
            uInit[i] = sim.eos[0]->primToConsv({vf,rho1*vf,rho2*(1-vf),-2,0.4});
        } else {
            uInit[i] = sim.eos[0]->primToConsv({1-vf,rho1*(1-vf),rho2*vf,2,0.4});
        }
    }

    std::ofstream initFile(dirname + "/0");
    initFile << "x alpha rho1 rho2 v p" << std::endl;
    for (std::vector<double>::size_type i = sim.ghosts(); i < uInit.size() - sim.ghosts(); i++) {
        double x = x0 + (i - sim.ghosts() + 0.5) * sim.get_dx();
        std::array<double,5> prim = sim.eos[0]->consvToPrim(uInit[i]);
        initFile << x << " " << prim[0] << " " << prim[1] << " " << prim[2] << " " << prim[3] << " " << prim[4] << std::endl;
    }
    initFile.close();
    

    // Initialising u field: //

    sim.init(uInit);

    
    // for (int j = 0; j < uInit.size(); j++){
    //     std::cout << uInit[j] << std::endl;
    // }

    sim.run();
    std::cout << "run finished" << std::endl;

    // Calculating exact result with exact riemann solver

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

    // Exact Riemann solver
    //double pStar = solution.exctRiemann()[2];
    //solution.exactRiemannSolution(pStar,discPosition);

    // end exact calculation

    return 0;
}