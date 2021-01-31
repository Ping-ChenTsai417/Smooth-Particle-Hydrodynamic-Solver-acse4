#include "SPH_2D.hpp"
#include "file_writer.hpp"

SPH_main domain;

int main(void)
{
    domain.set_values();                                        //Set simulation parameters
    domain.initialise_grid();                                    //initialise simulation grid

    domain.place_points(domain.min_x, domain.max_x);                //places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain

    write_file("original.vtp", &domain.particle_list);
    
    cout << domain.particle_list[152].x[0] << '\t' << domain.particle_list[152].x[1] << endl;

    //    domain.neighbour_iterate(&domain.particle_list[100]);        //finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle

    for (int i = 0; i != domain.particle_list.size(); i++)
    {
        domain.particle_list[i].rho = 1000;
    }

    domain.allocate_to_grid();                                    //needs to be called for each time step


    bool smooth;
    for (int i = 0; i != 5000; i++)
    {
        cout << "iteration " << i << endl;
        smooth = false;
        
        if (i % 10 == 0)
            smooth = true;
        domain.forward_euler(smooth);
//        cout << domain.particle_list[152].P << endl;
        
        if ( i % 50 == 0 )
        {
            string name = "example_" + to_string( (int) ( i / 50 ) ) + ".vtp";
            cout << name << endl;
            const char* cstr = name.c_str();
            write_file(cstr, &domain.particle_list);
        }
    }

    //write_file("example.vtp", &domain.particle_list);

    return 0;
}
