#include "SPH_2D.hpp"

SPH_main *SPH_particle::main_data;

// 计算每个黑线交点particle处于的蓝色grid里的grid的index
void SPH_particle::calc_index(void)
{
    for (int i = 0; i < 2; i++)
    {
        list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0*main_data->h));
    }
}

SPH_main::SPH_main()
{
    SPH_particle::main_data = this;
}

void SPH_main::set_values(void)
{
    inner_min_x[0] = min_x[0] = 0.0;
    inner_min_x[1] = min_x[1] = 0.0;

    inner_max_x[0] = max_x[0] = 20;
    inner_max_x[1] = max_x[1] = 10;

//    cout << "inner_min_x : " << inner_min_x[0] << " and " << inner_min_x[1] << endl;
//    cout << "inner_max_x : " << inner_max_x[0] << " and " << inner_max_x[1] << endl;
//
//    cout << "min_x : " << min_x[0] << " and " << min_x[1] << endl;
//    cout << "max_x : " << max_x[0] << " and " << max_x[1] << endl;
    
    dx = 0.2;
    
    h_fac = 1.3;
    h = dx*h_fac;
    t_max = 30;
    delta_t = 0.1 * h / C0;
}

void SPH_main::initialise_grid(void)
{
    for (int i = 0; i < 2; i++)
    {
        min_x[i] -= 3.0*h;
        max_x[i] += 3.0*h;                                                //add buffer for virtual wall particles

        // 蓝色网格每一维度可以有的最大限值
        max_list[i] = int((max_x[i] - min_x[i]) / (2.0*h) + 1.0);
    }

    // 建立grid，外层两个vector表示网格
    search_grid.resize(max_list[0]);
    for (int i=0;i<max_list[0];i++)
        search_grid[i].resize(max_list[1]);
    
//        cout << "inner_min_x : " << inner_min_x[0] << " and " << inner_min_x[1] << endl;
//        cout << "inner_max_x : " << inner_max_x[0] << " and " << inner_max_x[1] << endl;
//
//        cout << "min_x : " << min_x[0] << " and " << min_x[1] << endl;
//        cout << "max_x : " << max_x[0] << " and " << max_x[1] << endl;
}


// 添加初始化网格内的所有particles，以dx为间隔生成particle
void SPH_main::place_points(double *min, double *max)
{
    double x[2] = { min[0], min[1] };
    SPH_particle particle;

    
    while (x[0] <= max[0])
    {
        x[1] = min[1];
        while (x[1] <= max[1])
        {
            if (x[1] > 5 && x[1] < 10 && x[0] >= inner_min_x[0] && x[0] <= inner_max_x[0])
            {
                x[1] += dx;
                continue;
            }
            if (x[0] > 3 && x[0] < 20 && x[1] > 2 && x[1] < 5)
            {
                x[1] += dx;
                continue;
            }
            
            for (int i = 0; i < 2; i++)
                particle.x[i] = x[i];
            
            if ( ! ( x[0] > inner_min_x[0] && x[0] < inner_max_x[0] && x[1] > inner_min_x[1] &&  x[1] < inner_max_x[1] ) )
                particle.boundary_status = true;
            else
                particle.boundary_status = false;
            
            particle.calc_index();

            particle_list.push_back(particle);
            
//            if ( x[1] >= inner_min_x[1] && x[1] <= h)
//            {
//                cout << "-----------------------------\n";
//            }
//            if ( x[1] >= inner_max_x[1] - h && x[1] <= inner_max_x[1])
//            {
//                cout << "/////////////////////////////\n";
//            }

            x[1] += dx;
        }
        x[0] += dx;
    }
//    cout << "--------------------------size of particle is : " << particle_list.size() << endl;
}


void SPH_main::allocate_to_grid(void)                //needs to be called each time that all the particles have their positions updated
{
    // 删除全部蓝色grid里的所有黑色交点particle
    for (int i = 0; i < max_list[0]; i++)
        for (int j = 0; j < max_list[1]; j++)
            search_grid[i][j].clear();   // 删除这个网格点的particles

    for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
    {
        search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);    // 网格点加particles
    }
}




void SPH_main::neighbour_iterate(SPH_particle *part, bool smooth)                    //iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
    SPH_particle *other_part;
    double dist;            //distance between particles
    double dn[2];            //vector from 1st to 2nd particle

//    neighbour_particles.clear();
    
    
    double  dt_cfl = 10, dt_f = 10, dt_a = 10;

    // 只计算蓝色grid周围的一圈（8个）grid
    for (int i= part->list_num[0]-1;i<= part->list_num[0] + 1;i++)
        if (i>=0 && i<max_list[0])
            for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
                if (j >= 0 && j < max_list[1])
                {
                    // 周围的grid处理，计算particle和part grid里particle的距离
                    for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                    {
                        other_part = search_grid[i][j][cnt];

                        if (part != other_part)                    //stops particle interacting with itself
                        {
                            //Calculates the distance between potential neighbours
                            for (int n = 0; n < 2; n++)
                                dn[n] = part->x[n] - other_part->x[n];

                            dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
                            
                            if (dist < 2.*h)                    //only particle within 2h
                            {
//                                neighbour_particles.push_back( *other_part );
                                
                                // 计算
                                
                                double dW;
                                double m_j;
                                
                                double bracket_num1;
                                double bracket_num2;
                                
                                double e_ij;
                                
                                double first_item = 0;
                                double second_item = 0;
                                
                                m_j = dx * dx * rho0;
                                bracket_num1 = part->P / pow(part->rho, 2) + other_part->P / pow(other_part->rho, 2);
                                bracket_num2 = 1.0 / pow(part->rho, 2) + 1.0 / pow(other_part[j].rho, 2);
                                
                                double x_2 = pow((part->x[0] - other_part->x[0]), 2);
                                double y_2 = pow((part->x[1] - other_part->x[1]), 2);
                                
                                dist = sqrt(x_2 + y_2);
                                
                                dW = calculate_dW(dist);

                                e_ij = ( part->x[0] -  other_part->x[0] ) / dist;

                                first_item += m_j * bracket_num1 * dW * e_ij;
                                second_item += m_j * bracket_num2 * dW * (part->v[0] - other_part[j].v[0]) / dist;
                                
                                
                                part->a[0] += mu * second_item - first_item;
                                
                                bracket_num1 = part->P / pow(part->rho, 2) + other_part->P / pow(other_part->rho, 2);
                                bracket_num2 = 1.0 / pow(part->rho, 2) + 1.0 / pow(other_part[j].rho, 2);
                          
                                x_2 = pow((part->x[0] - other_part->x[0]), 2);
                                y_2 = pow((part->x[1] - other_part->x[1]), 2);
                              
                                dist = sqrt(x_2 + y_2);
                              
                                dW = calculate_dW(dist);

                                e_ij = ( part->x[1] -  other_part->x[1] ) / dist;


                                first_item = m_j * bracket_num1 * dW * e_ij;
                                 second_item = m_j * bracket_num2 * dW * (part->v[1] - other_part[j].v[1]) / dist;
                                
                                part->a[1] += mu * second_item - first_item;
                                
                                
                                dW = calculate_dW(dist);
                                    
                                double dot_product = ( (part->v[0] - other_part->v[0]) * (part->x[0] - other_part->x[0]) + (part->v[1] - other_part->v[1]) * (part->x[1] - other_part->x[1]) ) / dist;

                                part->D += m_j * dW * dot_product;
                                
                                // dt_cfl
                                double v_ij = sqrt( pow( part->v[0] - other_part->v[0], 2) + pow( part->v[1] - other_part->v[1], 2) );
                                double tmp_cfl = h / v_ij;
                                if (tmp_cfl < dt_cfl)
                                    dt_cfl = tmp_cfl;
                                
                                // dt_f
                                double a_i = sqrt( pow( part->a[0], 2) + pow( part->a[1], 2) );
                                double tmp_f = sqrt( h / a_i);
                                if (tmp_f < dt_f)
                                    dt_f = tmp_f;
                                
                                // dt_a
                                double tmp_a = h / C0 / sqrt( pow( part->rho / rho0, (gama - 1) / 2.0 ) );
                                if (tmp_a < dt_a)
                                    dt_a = tmp_a;
                                                                
                            }
                        }
                    }
                }
    
//    delta_t = min( min( dt_cfl, dt_f), dt_a );
}


//    -------------------------------------------------------------------------------
double SPH_main::calculate_W(double r)
{
    double q = r / h;
    double tmp;
    
    if ( q >= 0 && q <= 1)
        tmp =  1 - 1.5 * q * q + 0.75 * q * q * q;
    else if (q > 1 && q <= 2)
        tmp =  0.25 * pow(( 2 - q ), 3);
    else
        return 0;
    
    return (10 / 7.0 / PI / h / h) * tmp;
}

double SPH_main::calculate_dW(double r)
{
    double q = r / h;
    double tmp;
    
    if ( q >= 0 && q <= 1)
        tmp = - 3 * q + 2.25 * q * q;
    else if (q > 1 && q <= 2)
        tmp =  - 0.75 * pow(( 2 - q ), 2);
    else
        return 0;
    
    return (10 / 7.0 / PI / h / h / h) * tmp;
}

double distance(SPH_particle & i, SPH_particle & j)
{
    double x_2 = pow((i.x[0] - j.x[0]), 2);
    double y_2 = pow((i.x[1] - j.x[1]), 2);
    
    return sqrt(x_2 + y_2);
}

double SPH_main::calculate_rho_i(SPH_particle & part)
{
    double dist;
    double W;
    
    double numerator_sum = 0;
    double denominator_sum = 0;
    
    int num = (int) neighbour_particles.size();
    
    for (int j = 0; j != num; j++)
    {
        dist = distance( part, neighbour_particles[j] );
        W = calculate_W(dist);
        
        numerator_sum += W;
        denominator_sum += W / neighbour_particles[j].rho;
    }
    
    dist = 0;
    W = calculate_W(dist);
    
    numerator_sum += W;
    denominator_sum += W / part.rho;
    
    return ( numerator_sum / denominator_sum );
}



void SPH_main::forward_euler(bool smooth)
{
    allocate_to_grid();
    
    if (smooth)
    {
        SPH_particle * other_part;
        double dn[2];
        double dist;

        double W;
        for (int ii = 0; ii != particle_list.size(); ii++)
        {
            double numerator_sum = 0;
            double denominator_sum = 0;
            for (int i = particle_list[ii].list_num[0]-1;i<= particle_list[ii].list_num[0] + 1;i++)
                    if (i>=0 && i<max_list[0])
                        for (int j = particle_list[ii].list_num[1] - 1; j <= particle_list[ii].list_num[1] + 1; j++)
                            if (j >= 0 && j < max_list[1])
                            {
                                // 周围的grid处理，计算particle和part grid里particle的距离
                                for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                                {
                                    other_part = search_grid[i][j][cnt];
                                    //Calculates the distance between potential neighbours
                                    for (int n = 0; n < 2; n++)
                                        dn[n] = particle_list[ii].x[n] - other_part->x[n];

                                    dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
                                    
                                    if (dist < 2.*h)                    //only particle within 2h
                                    {
                                        W = calculate_W(dist);
                                        numerator_sum += W;
                                        denominator_sum += W / other_part->rho;
                                    }
                                    
                                }
                            }
            particle_list[ii].rho = numerator_sum / denominator_sum;
        }
    }
                                        

    for ( int i = 0; i != particle_list.size(); i++ )
    {
        particle_list[i].a[0] = 0;
        if (particle_list[i].boundary_status)
            particle_list[i].a[1] = 0;
        else
            particle_list[i].a[1] = G;
        particle_list[i].D = 0;
        neighbour_iterate( &( particle_list[i] ) , smooth);
    }
    
    for ( int i = 0; i != particle_list.size(); i++ )
    {
        if (particle_list[i].boundary_status == true)
        {
            particle_list[i].v[0] = 0;
            particle_list[i].v[1] = 0;
            particle_list[i].boundary_status = true;
        }
        else
        {
            bool vertical_wall = false, horizontal_wall = false;
            particle_list[i].x[0] = particle_list[i].x[0] + delta_t * particle_list[i].v[0];
            if (particle_list[i].x[0] < inner_min_x[0] || particle_list[i].x[0] > inner_max_x[0])
            {
                particle_list[i].x[0] = particle_list[i].x[0] - delta_t * particle_list[i].v[0];
                particle_list[i].v[0] = - velocity_lost_rate * particle_list[i].v[0];
                horizontal_wall = true;
            }
            particle_list[i].x[1] = particle_list[i].x[1] + delta_t * particle_list[i].v[1];
            if (particle_list[i].x[1] < inner_min_x[1] || particle_list[i].x[1] > inner_max_x[1])
            {
                particle_list[i].x[1] = particle_list[i].x[1] - delta_t * particle_list[i].v[1];
                particle_list[i].v[1] = - velocity_lost_rate * particle_list[i].v[1];
                vertical_wall = true;
            }
            
            if (!vertical_wall)
                particle_list[i].v[1] = particle_list[i].v[1] + delta_t * particle_list[i].a[1];
            if (!horizontal_wall)
                particle_list[i].v[0] = particle_list[i].v[0] + delta_t * particle_list[i].a[0];
            
        }
    //    cout << endl << particle_list[i].v[1] << endl;
        particle_list[i].rho = particle_list[i].rho + delta_t * particle_list[i].D;
        particle_list[i].calculate_P();
        particle_list[i].calc_index();
        
    }
    
    t_max += delta_t;
}


double SPH_main::average_neighbour_velocity(string direction)
{
    double sum = 0;
    for (int i = 0; i != neighbour_particles.size(); i++)
    {
        if (direction == "horizontal")
            sum += neighbour_particles[i].v[0];
        if (direction == "vertical")
            sum += neighbour_particles[i].v[1];
    }
    return (sum / neighbour_particles.size());
}


void SPH_main::predictor_corrector(bool smoo)
{
    cout << delta_t << endl;
    allocate_to_grid();

    for (int n = 0; n != 2; n++)
    {
//        for (int i = 0; i != max_list[0]; i++)
//            for (int j = 0; j != max_list[1]; j++)
//                search_grid[i][j].clear();
        
        if ( n == 0 )
        {
            
            if (smoo)
            {
                SPH_particle * other_part;
                double dn[2];
                double dist;

                double W;
                for (int ii = 0; ii != particle_list.size(); ii++)
                {
                    double numerator_sum = 0;
                    double denominator_sum = 0;
                    for (int i = particle_list[ii].list_num[0]-1;i<= particle_list[ii].list_num[0] + 1;i++)
                            if (i>=0 && i<max_list[0])
                                for (int j = particle_list[ii].list_num[1] - 1; j <= particle_list[ii].list_num[1] + 1; j++)
                                    if (j >= 0 && j < max_list[1])
                                    {
                                        // 周围的grid处理，计算particle和part grid里particle的距离
                                        for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                                        {
                                            other_part = search_grid[i][j][cnt];
                                            //Calculates the distance between potential neighbours
                                            for (int n = 0; n < 2; n++)
                                                dn[n] = particle_list[ii].x[n] - other_part->x[n];

                                            dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
                                            
                                            if (dist < 2.*h)                    //only particle within 2h
                                            {
                                                W = calculate_W(dist);
                                                numerator_sum += W;
                                                denominator_sum += W / other_part->rho;
                                            }
                                            
                                        }
                                    }
                    particle_list[ii].rho = numerator_sum / denominator_sum;
                }
            }
            
            for (int i = 0; i != particle_list.size(); i++)
            {
                particle_list[i].a[0] = 0;
                if (particle_list[i].boundary_status)
                    particle_list[i].a[1] = 0;
                else
                    particle_list[i].a[1] = G;
                particle_list[i].D = 0;
                neighbour_iterate( &particle_list[i], smoo);
            }
            
            for (int i = 0; i != particle_list.size(); i++)
            {

                
                
                if (particle_list[i].boundary_status == true)
                {
                    particle_list[i].v[0] = 0;
                    particle_list[i].v[1] = 0;
                }
                else
                {
                    particle_list[i].prev_x[0] = particle_list[i].x[0];
                    particle_list[i].prev_x[1] = particle_list[i].x[1];
                    particle_list[i].prev_v[0] = particle_list[i].v[0];
                    particle_list[i].prev_v[1] = particle_list[i].v[1];
                    particle_list[i].prev_rho = particle_list[i].rho;
                    
                    
                    bool vertical_wall = false, horizontal_wall = false;
                    particle_list[i].x[0] = particle_list[i].x[0] + 0.5 * delta_t * particle_list[i].v[0];
                    if ( particle_list[i].x[0] < inner_min_x[0] || particle_list[i].x[0] > inner_max_x[0] )
                    {
                         particle_list[i].x[0] =  particle_list[i].x[0] - 0.5 * delta_t * particle_list[i].v[0];
                        particle_list[i].v[0] = - velocity_lost_rate * particle_list[i].v[0];
                        horizontal_wall = true;
                    }
                    
                    particle_list[i].x[1] = particle_list[i].x[1] + 0.5 * delta_t * particle_list[i].v[1];
                    if ( particle_list[i].x[1] < inner_min_x[1] || particle_list[i].x[1] > inner_max_x[1] )
                    {
                        particle_list[i].x[1] =  particle_list[i].x[1] - 0.5 * delta_t * particle_list[i].v[1];
                        particle_list[i].v[1] = - velocity_lost_rate * particle_list[i].v[1];
                        vertical_wall = true;
                    }
                    
                    if (!horizontal_wall)
                        particle_list[i].v[0] = particle_list[i].v[0] + 0.5 * delta_t * particle_list[i].a[0];

                    if (!vertical_wall)
                        particle_list[i].v[1] = particle_list[i].v[1] + 0.5 * delta_t * particle_list[i].a[1];
                }
                
                particle_list[i].rho = particle_list[i].rho + 0.5 * delta_t * particle_list[i].D;
                particle_list[i].calc_index();
                particle_list[i].calculate_P();
            }
        }
        
        if ( n == 1 )
        {
            for (int i = 0; i != particle_list.size(); i++)
            {
                particle_list[i].a[0] = 0;
                if (particle_list[i].boundary_status)
                    particle_list[i].a[1] = 0;
                else
                    particle_list[i].a[1] = G;
                particle_list[i].D = 0;
                
            }
            
            
            for (int i = 0; i != particle_list.size(); i++)
            {
                if ( !particle_list[i].boundary_status )
                {
                    double x_[2], v_[2];
                    bool vertical_wall = false, horizontal_wall = false;
                    x_[0] = particle_list[i].prev_x[0] + 0.5 * delta_t * particle_list[i].v[0];
                    x_[1] = particle_list[i].prev_x[1] + 0.5 * delta_t * particle_list[i].v[1];
                    
                    particle_list[i].x[0] = 2 * x_[0] - particle_list[i].prev_x[0];
                    if ( particle_list[i].x[0] < inner_min_x[0] || particle_list[i].x[0] > inner_max_x[0] )
                    {
                         particle_list[i].x[0] =  particle_list[i].prev_x[0];
                        particle_list[i].v[0] = - velocity_lost_rate * particle_list[i].prev_v[0];
                        horizontal_wall = true;
                    }
                    
                    particle_list[i].x[1] = 2 * x_[1] - particle_list[i].prev_x[1];
                    if ( particle_list[i].x[1] < inner_min_x[1] || particle_list[i].x[1] > inner_max_x[1] )
                    {
                         particle_list[i].x[1] =  particle_list[i].prev_x[1];
                        particle_list[i].v[1] = - velocity_lost_rate * particle_list[i].prev_v[1];
                        horizontal_wall = true;
                    }
                    if (!horizontal_wall)
                    {
                        v_[0] = particle_list[i].prev_v[0] + 0.5 * delta_t * particle_list[i].a[0];
                        particle_list[i].v[0] = 2 * v_[0] - particle_list[i].prev_v[0];
                    }
                    if (!vertical_wall)
                    {
                        v_[1] = particle_list[i].prev_v[1] + 0.5 * delta_t * particle_list[i].a[1];
                        particle_list[i].v[1] = 2 * v_[1] - particle_list[i].prev_v[1];
                    }
                    
                }
                else
                {
                    particle_list[i].v[0] = 0;
                    particle_list[i].v[1] = 0;
                }
                double rho_;
                
                rho_ = particle_list[i].prev_rho + 0.5 * delta_t * particle_list[i].D;
                particle_list[i].rho = 2 * rho_ - particle_list[i].prev_rho;
                particle_list[i].calc_index();
                particle_list[i].calculate_P();

            }
        }
    }
}



















//
//double SPH_main::Jones_potential(SPH_particle & part, string direction)
//{
//    double extra_force = 0;
////    cout << part.x[0] << '\t' << part.x[1] << " : ---------------------------\n";
//    for (int i = 0; i != neighbour_particles.size(); i++)
//    {
////        cout << neighbour_particles[i].x[0] << '\t' << neighbour_particles[i].x[1] << '\t' << "staus is : " << neighbour_particles[i].boundary_status << endl;
//        if (neighbour_particles[i].boundary_status == true)
//        {
//            double r = distance(part, neighbour_particles[i]);
//            double H = 2 * h;
//            double r0 = 0.9 * h;
//
////            cout << "--------------------------------\n";
////            cout << r << endl;
//
////            cout << r0 / r << "--------------------------------";
//            if (direction == "vertical")
//            {
//                 double x = part.x[1] - neighbour_particles[i].x[1];
//                extra_force += index_jones * jone_G * H * ( pow( r0 / r, 4) - pow( r0 / r, 2) ) / pow(r, 2)  * x;
//            }
//            if (direction == "horizontal")
//            {
//                double x = part.x[0] - neighbour_particles[i].x[0];
//                extra_force += index_jones * G * H * ( pow( r0 / r, 4) - pow( r0 / r, 2) / pow(r, 2) ) * x;
//            }
//
//        }
//    }
//    return extra_force;
//}

//
//double SPH_main::nana_force(SPH_particle& part, string direction)
//{
//    double dist = 0;
//    //top and bottom
//    if (direction == "down")
//    {
//        dist = part.x[1] - inner_min_x[1];
//    }
//    if (direction == "top")
//    {
//        dist = part.x[1] - inner_max_x[1];
//    }
//    if (direction == "right")
//    {
//        dist = part.x[0] - inner_max_x[0];
//    }
//    if (direction == "left")
//    {
//        dist = part.x[0] - inner_min_x[0];
//    }
//    double f_dist = pow(0.5 * dx / dist, 6) - 1;
//    return Pref * dx * dx * f_dist;
//}
