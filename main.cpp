#include "SPH_2D.hpp"
#include "file_writer.hpp"

SPH_main domain;

int main(void)
{
    domain.set_values();                                        //Set simulation parameters
    domain.initialise_grid();                                    //initialise simulation grid

    domain.place_points(domain.min_x,domain.max_x);                //places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
    
    write_file("original.vtp", &domain.particle_list);
    

//    domain.neighbour_iterate(&domain.particle_list[100]);        //finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle
    
    for ( int i = 0; i != domain.particle_list.size(); i++ )
    {
        domain.particle_list[i].rho = 1000;
    }
    
    domain.allocate_to_grid();                                    //needs to be called for each time step

    
    for (int i = 0; i != 5; i++)
    {
        if (i != 0 && i % 10 == 0)
        {
            for (int j = 0; j != domain.particle_list.size(); j++)
                domain.calculate_rho_i(domain.particle_list[j]);
        }
        domain.forward_euler();
//        if (i == 0)
//        {
//            for ( int j = 0; j != domain.particle_list.size(); j++ )
//                cout << domain.particle_list[j].rho << '\t';
//        }
//        cout << i << " iteration : \n";
    }
    
    write_file("example.vtp", &domain.particle_list);
    
    return 0;
}
