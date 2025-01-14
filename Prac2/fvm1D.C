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
    int nCells = 500;
    double x0{0}, x1{1};
    double startTime = 0.0, endTime = 5e-5;

    double cour{0.8};
    std::cout << "Enter CFL number:"; std::cin >> cour;

    // construct EOS object and give to sovler
    stiffenedGas stiffgas(7.15,3e8,0);
    //idealGas idgas(1.4);
    solver sim(x0,x1,startTime,endTime,nCells,2,cour);
    sim.setEOS(&stiffgas);

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

    sim.setWriteInterval(5e-6);

    // Set the flux function
    sim.flux = [&sim](std::array<double,3> input){
        return sim.fEuler(input);
    }; // because fEuler isnt static

    sim.slopeLim = [&sim](std::array<double,3> input){
        return sim.minbee(input);
    }; // because fEuler isnt static

    std::vector<std::array<double,3>> uInit; uInit.resize(sim.uPlus1.size());
    std::cout << "n ghost cells=" << sim.ghosts() << std::endl;
    for (std::vector<double>::size_type i=sim.ghosts(); i<uInit.size()-sim.ghosts();i++){
        double x = x0 + (i-sim.ghosts()+0.5)*sim.get_dx();
        // for (int var = 0; var < 3; ++var){
        //     if (var == 0 || var == 2){
        //         uInit[i][var] = 1;
        //     } else {
        //         uInit[i][var] = 0;
        //     }
        // }
        if (x <= 0.5){
            for (int var = 0; var < 3; ++var){
                if (var == 0){
                    uInit[i][var] = 1500; //1500
                } else if (var == 1){
                    uInit[i][var] = 0;
                } else {
                    uInit[i][var] = 3000*101325; // 3000*101325
                }
            }
        } else {
            for (int var = 0; var < 3; ++var){
                if (var == 0){
                    uInit[i][var] = 1000;
                } else if (var == 1){
                    uInit[i][var] = 0;
                } else {
                    uInit[i][var] = 101325;
                }
            }
        }

    }


    std::ofstream initFile1(dirname + "/rho/0");
    for (std::vector<double>::size_type i=sim.ghosts(); i<uInit.size()-sim.ghosts();i++){
        double x = x0 + (i-sim.ghosts()+0.5)*sim.get_dx();
        initFile1 << x << " " << uInit[i][0] << std::endl;
    }
    initFile1.close();

    std::ofstream initFile2(dirname + "/v/0");
    for (std::vector<double>::size_type i=sim.ghosts(); i<uInit.size()-sim.ghosts();i++){
        double x = x0 + (i-sim.ghosts()+0.5)*sim.get_dx();
        initFile2 << x << " " << uInit[i][1] << std::endl;
    }
    initFile2.close();

    std::ofstream initFile3(dirname + "/p/0");
    for (std::vector<double>::size_type i=sim.ghosts(); i<uInit.size()-sim.ghosts();i++){
        double x = x0 + (i-sim.ghosts()+0.5)*sim.get_dx();
        initFile3 << x << " " << uInit[i][2] << std::endl;
    }
    initFile3.close();
    //std::cout << u[i] << " ";
    

    // Initialising u field: //

    for (std::vector<double>::size_type i=sim.ghosts(); i<uInit.size()-sim.ghosts();i++){
        uInit[i] = sim.eos->primToConsv(uInit[i]);
    }
    sim.init(uInit);

    
    // for (int j = 0; j < uInit.size(); j++){
    //     std::cout << uInit[j] << std::endl;
    // }

    sim.run();
    std::cout << "run finished" << std::endl;

    /*
    t = startTime;

    std::string errFileName = dirname + "/err";
    std::ofstream errFile(errFileName);
    if (errFile.is_open()){
        // writing data
        int counter = 0;
        do{
            errFile << t << " " << errs[counter] << std::endl;
            counter++;
            t = t + dt;
        } while (t < endTime);
    }else{
        errFile << "failed to write errors." << std::endl;
    }
    errFile.close();
    */

    return 0;
}