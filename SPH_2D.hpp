#ifndef SPH_2D_hpp
#define SPH_2D_hpp

#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

#define mu 0.001
#define G -9.81
#define PI 3.1415926
#define gama 7
#define C0 20
#define rho0 1000
#define boundary_influence_length_factor 0.5


#define q_factor 0.7
#define pow_index 3
#define total_factor 0.001

#define index_jones 5

using namespace std;

class SPH_main;

// 不用加东西
class SPH_particle
{
public:
    double x[2], v[2];                //position and velocity
    double rho, P;                    //density and pressure

    static SPH_main *main_data;        //link to SPH_main class so that it can be used in calc_index

    int list_num[2];                //index in neighbour finding array
    
    bool boundary_status = false;
    
    void calc_index(void);
    
    
//    SPH_particle();
//
//    // copy constructor
//    SPH_particle(SPH_particle & A);
//    // overload copy assignment
//    SPH_particle & operator= (SPH_particle & A);
//
//    SPH_particle(SPH_particle && A);
//
//    SPH_particle & operator= (SPH_particle && A);

};


class SPH_main
{
public:
    // 不用加东西
    SPH_main();

    // 不用加东西
    void set_values(void);
    void initialise_grid(void);

    // 不用加东西
    void place_points(double *min, double *max);
    
    // 加东西
    void allocate_to_grid(void);            //allocates all the points to the search grid (assumes that index has been appropriately updated)
    
    // 加东西
    void neighbour_iterate(SPH_particle *part);
    
    vector<SPH_particle> neighbour_particles;

    double h;                                //smoothing length
    double h_fac;
    double dx;                                //particle initial spacing

    double delta_t = 0.02;
    
    double min_x[2], max_x[2];                //dimensions of simulation region
    double inner_max_x[2], inner_min_x[2];

    int max_list[2];

    vector<SPH_particle> particle_list;                        //list of all the particles

    vector<vector<vector<SPH_particle*> > > search_grid;        //Outer 2 are the grid, inner vector is the list of pointers in each cell
    
    double extra_force = 0;
    double D;
//    -------------------------------------------------------------------------------
    double calculate_W(double r);
    
    double calculate_dW(double r);
    
    friend double distance(SPH_particle & i, SPH_particle & j);
    
    double calculate_rho_i(SPH_particle & part);
    
    double dv0_dt(SPH_particle & part);
    double dv1_dt(SPH_particle & part);
    
    double drho_dt(SPH_particle & part);
    
    void calculate_P();
    
    // updator
    SPH_particle whole_f(SPH_particle & part);
    
    void forward_euler();
    
    double Jones_potential(SPH_particle & part, string direction);
};



#endif /* SPH_2D_hpp */
