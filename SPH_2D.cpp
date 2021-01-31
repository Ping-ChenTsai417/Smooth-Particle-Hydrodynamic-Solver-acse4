#include "SPH_2D.hpp"

SPH_main *SPH_particle::main_data;

// 计算每个黑线交点particle处于的蓝色grid里的grid的index
void SPH_particle::calc_index(void)
{
    for (int i = 0; i < 2; i++)
    {
        list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0*main_data->h));
//        cout << list_num[i] << '\t';
    }
//    cout << endl;
}

SPH_main::SPH_main()
{
    SPH_particle::main_data = this;
}

void SPH_main::set_values(void)
{
    inner_min_x[0] = min_x[0] = 0.0;
    inner_min_x[1] = min_x[1] = 0.0;

    inner_max_x[0] = max_x[0] = 2.0;
    inner_max_x[1] = max_x[1] = 1.0;

//    cout << "inner_min_x : " << inner_min_x[0] << " and " << inner_min_x[1] << endl;
//    cout << "inner_max_x : " << inner_max_x[0] << " and " << inner_max_x[1] << endl;
//
//    cout << "min_x : " << min_x[0] << " and " << min_x[1] << endl;
//    cout << "max_x : " << max_x[0] << " and " << max_x[1] << endl;
    
    dx = 0.02;
    
    h_fac = 1.3;
    h = dx*h_fac;
}

void SPH_main::initialise_grid(void)
{
    for (int i = 0; i < 2; i++)
    {
        min_x[i] -= 2.0*h;
        max_x[i] += 2.0*h;                                                //add buffer for virtual wall particles

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
            if (x[1] > 0.5 && x[1] < 1.0 && x[0] >= inner_min_x[0] && x[0] <= inner_max_x[0])
            {
                x[1] += dx;
                continue;
            }
            if (x[0] > 0.3 && x[0] < 2.0 && x[1] > 0.2 && x[1] <0.5)
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

//    cout << "--------------------------size of particle is : " << particle_list.size() << endl;
    for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
    {
       // 利用list_num表示grid的x，y坐标，然后对应每个grid添加新的particles
//        cout << "number of " << cnt << endl;
//        cout << particle_list[cnt].list_num[0] << '\t' << particle_list[cnt].list_num[1] << endl;
        search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);    // 网格点加particles
    }
}


void SPH_main::neighbour_iterate(SPH_particle *part)                    //iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
    SPH_particle *other_part;
    double dist;            //distance between particles
    double dn[2];            //vector from 1st to 2nd particle

    neighbour_particles.clear();
    
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
                                neighbour_particles.push_back( *other_part );
                            }
                            
//                            cout<< other_part->P << '\t';
                        }
                    }
                }
//    for ( int i = 0; i != neighbour_particles.size(); i++ )
//        cout << neighbour_particles[i].P << '\t';
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
    
    return (10 / 7.0 / PI / h / h) * tmp;
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


double SPH_main::dv0_dt(SPH_particle & part)
{
    double dist;
    double dW;
    double m_j;
    
    double bracket_num1;
    double bracket_num2;
    
    double e_ij;
    
    double first_item = 0;
    double second_item = 0;
    
    int num = (int) neighbour_particles.size();
    
    for (int j = 0; j != num; j++)
    {
//        cout << neighbour_particles[j].P << '\t';
        m_j = dx * dx * neighbour_particles[j].rho;
        bracket_num1 = part.P / pow(part.rho, 2) + neighbour_particles[j].P / pow(neighbour_particles[j].rho, 2);
        bracket_num2 = 1.0 / pow(part.rho, 2) + 1.0 / pow(neighbour_particles[j].rho, 2);
        
        dist = distance( part, neighbour_particles[j] );
        dW = calculate_dW(dist);

        e_ij = ( part.x[0] -  neighbour_particles[j].x[0] ) / dist;

        first_item += m_j * bracket_num1 * dW * e_ij;
        second_item += m_j * bracket_num2 * dW * (part.v[0] - neighbour_particles[j].v[0]) / dist;
    }
    
    double extra_a = 0;
//    cout << ( -first_item + mu * second_item ) << '\t';
    // 左边界
    if ( (part.x[0] <= inner_min_x[0] + boundary_influence_length_factor * dx) && ( part.x[0] >= inner_min_x[0]) )
    {
//        double abs_dist = part.x[0] - inner_min_x[0];
//        double q = abs_dist / ( q_factor * dx );
////        cout << q << '\t';
//
//        if (q < 0.04)
//            q = 0.04;
//
//        double Pref = ( rho0 * C0 * C0 / 7.0 ) * (pow(1.02, pow_index) - 1);
//        extra_a = total_factor * Pref / rho0 * ( pow((1/q),4) - pow(1/q,2) ) / abs_dist;
        extra_a = Jones_potential(part, "horizontal");
    }
    // 右边界
    else if ( (part.x[0] >= inner_max_x[0] - boundary_influence_length_factor * dx) && (part.x[0] <= inner_max_x[0]) )
    {
//        double abs_dist =  inner_max_x[0] - part.x[0];
//        double q = abs_dist / ( q_factor * dx );
////        cout << q << '\t';
//
//
//        if (q < 0.04)
//            q = 0.04;
//
//        double Pref = ( rho0 * C0 * C0 / 7.0 ) * (pow(1.02, pow_index) - 1);
//        extra_a = - total_factor * Pref / rho0 * ( pow((1/q),4) - pow(1/q,2) ) / abs_dist;
        extra_a = Jones_potential(part, "horizontal");
//        cout << extra_a << '\t';
    }
    
    return ( -first_item + mu * second_item + extra_a);
}


 double SPH_main::dv1_dt(SPH_particle & part)
 {
     double dist;
     double dW;
     double m_j;
     
     double e_ij;
     
     double bracket_num1;
     double bracket_num2;
     
     double first_item = 0;
     double second_item = 0;
     
     int num = (int) neighbour_particles.size();
     
     for (int j = 0; j != num; j++)
     {
         m_j = dx * dx * neighbour_particles[j].rho;
         
         bracket_num1 = part.P / pow(part.rho, 2) + neighbour_particles[j].P / pow(neighbour_particles[j].rho, 2);
         
         bracket_num2 = 1.0 / pow(part.rho, 2) + 1.0 / pow(neighbour_particles[j].rho, 2);
         
         dist = distance( part, neighbour_particles[j] );
         e_ij = (  part.x[1] - neighbour_particles[j].x[1]) / dist;
         
//         cout << e_ij << '\t' << endl;
         
         dW = calculate_dW(dist);
         
         first_item += m_j * bracket_num1 * dW * e_ij;
         second_item += m_j * bracket_num2 * dW * (part.v[1] - neighbour_particles[j].v[1]) / dist;
     }
//     cout << ( -first_item + mu * second_item + G) << '\t';
     // 底边界
     double extra_a = 0;
     
     if ( (part.x[1] <= boundary_influence_length_factor * dx) && (part.x[1] >= inner_min_x[0]) )
     {
//         double abs_dist = part.x[1] - inner_min_x[1];
//         double q = abs_dist / ( q_factor * dx );
////         cout << q << '\t';
//         if (q < 0.04)
//             q = 0.04;
//
//         double Pref = ( rho0 * C0 * C0 / 7.0 ) * (pow(1.02, pow_index) - 1);
//         extra_a = total_factor * Pref / rho0 * ( pow((1/q),4) - pow(1/q,2) ) / abs_dist;
         extra_a = Jones_potential(part, "vertical");
//              cout << extra_a << '\t';
     }
     // 上边界 一般不应该有值更新
     else if ( (part.x[1] >= inner_max_x[1] - boundary_influence_length_factor * dx) && (part.x[1] <= inner_max_x[1]) )
     {
//            double abs_dist =  inner_max_x[1] - part.x[1];
//            double q = abs_dist / ( q_factor * dx );
////         cout << q << '\t';
//
//          if (q < 0.04)
//              q = 0.04;
//
//          double Pref = ( rho0 * C0 * C0 / 7.0 ) * (pow(1.02, pow_index) - 1);
//          extra_a = - total_factor * Pref / rho0 * ( pow((1/q),4) - pow(1/q,2) ) / abs_dist;
          extra_a = Jones_potential(part, "vertical");
     }
//     cout << extra_a << '\t';
     return ( -first_item + mu * second_item + G + extra_a);
 }

double SPH_main::drho_dt(SPH_particle & part)
{
    double dist;
    double dW;
    double m_j;
    
    double sum = 0;
    
    double dot_product;
    
    int num = (int) neighbour_particles.size();
    
    for (int j = 0; j != num; j++)
    {
        m_j = dx * dx * neighbour_particles[j].rho;
//        cout << m_j << '\t';
        dist = distance( part, neighbour_particles[j] );
        dW = calculate_dW(dist);
        
        dot_product = ( (part.v[0] - neighbour_particles[j].v[0]) * (part.x[0] - neighbour_particles[j].x[0]) + (part.v[1] - neighbour_particles[j].v[1]) * (part.x[1] - neighbour_particles[j].x[1]) ) / dist;
//        cout << part.v[0] << '\t' << neighbour_particles[j].v[0] << '\t' ;//<< part.v[1] << '\t' << neighbour_particles[j].v[1] << '\t';
        sum += m_j * dW * dot_product;
    }
    
    
//    cout << sum << '\t';
    return sum;
}

void SPH_main::calculate_P()
{
    double B = rho0 * C0 * C0 / gama;
    for (int i = 0; i != particle_list.size(); i++)
    {
        particle_list[i].P = B * ( pow( (particle_list[i].rho / rho0) , gama) - 1 ) ;
    }
}


 SPH_particle SPH_main::whole_f(SPH_particle & part)
{
    SPH_particle new_part;
//    cout << part.boundary_status << endl;
    if (part.boundary_status == true)
    {
        new_part.x[0] = part.x[0];
        new_part.x[1] = part.x[1];
        new_part.v[0] = 0;
        new_part.v[1] = 0;
        new_part.boundary_status = true;
    }
    else
    {
        new_part.x[0] = part.x[0] + delta_t * part.v[0];
        new_part.x[1] = part.x[1] + delta_t * part.v[1];
        new_part.v[0] = part.v[0] + delta_t * dv0_dt(part);
        new_part.v[1] = part.v[1] + delta_t * dv1_dt(part);
    }
    new_part.P = part.P;
//    cout << endl << part.v[1] << endl;
    new_part.rho = part.rho + delta_t * drho_dt(part);

    return new_part;
}


void SPH_main::forward_euler()
{
    vector<SPH_particle> new_particle_list;
    
    // 必须在每次update前更新P
    calculate_P();
//    for ( int i = 0; i != particle_list.size(); i++ )
//    {
//        particle_list[i].P = 10;
//    }
    
    for ( int i = 0; i != particle_list.size(); i++ )
    {
//        cout << particle_list[i].boundary_status << '\t';
//        cout << i << "iteration --------------------------";
//        cout << particle_list[i].x[0] << '\t' << particle_list[i].x[1] << endl;
        neighbour_iterate( &( particle_list[i] ) );
        
//        if ( (particle_list[i].x[1] <= boundary_influence_length_factor * h) && (particle_list[i].x[1] >= inner_min_x[0]) )
        {
//            for (int j = 0; j != neighbour_particles.size(); j++)
//            {
//                cout << neighbour_particles[j].v[1] << '\t';
//            }
//            cout << endl;
        }
//        cout << particle_list[i].boundary_status << '\t';
        SPH_particle new_particle = whole_f( particle_list[i] );
        new_particle_list.push_back( new_particle );
    }
    
    
    
    for ( int i = 0; i != particle_list.size(); i++ )
    {
        particle_list[i] = new_particle_list[i];
        particle_list[i].calc_index();
    }
    
    // 目前不能加，因为我们的particle会越界，所以越界的计算index，不能正常运行allocate_to_grid
    // 也就是说，allocate_to_grid可以当作检验我们的particle是否越界的准则
//    allocate_to_grid();
    
}


double SPH_main::Jones_potential(SPH_particle & part, string direction)
{
    double extra_force = 0;
//    cout << part.x[0] << '\t' << part.x[1] << " : ---------------------------\n";
    for (int i = 0; i != neighbour_particles.size(); i++)
    {
//        cout << neighbour_particles[i].x[0] << '\t' << neighbour_particles[i].x[1] << '\t' << "staus is : " << neighbour_particles[i].boundary_status << endl;
        if (neighbour_particles[i].boundary_status == true)
        {
            double r = distance(part, neighbour_particles[i]);
            double H = 2 * h;
            double r0 = boundary_influence_length_factor * dx;
            
//            cout << "--------------------------------\n";
//            cout << r << endl;
            
            if (r < r0)
            if (direction == "vertical")
            {
                 double x = part.x[1] - neighbour_particles[i].x[1];
                extra_force += index_jones * G * H * ( pow( r0 / r, 4) - pow( r0 / r, 2) / pow(r, 2) ) * x;
            }
            if (direction == "horizontal")
            {
                double x = part.x[0] - neighbour_particles[i].x[0];
                extra_force += index_jones * G * H * ( pow( r0 / r, 4) - pow( r0 / r, 2) / pow(r, 2) ) * x;
            }

        }
    }
    return extra_force;
}

