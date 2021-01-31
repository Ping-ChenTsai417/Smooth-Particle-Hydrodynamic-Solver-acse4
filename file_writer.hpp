//
//  file_writer.hpp
//  ACSE4_Group_Project3
//
//  Created by 谢逊 on 2020/2/24.
//  Copyright © 2020 谢逊. All rights reserved.
//

#ifndef file_writer_hpp
#define file_writer_hpp

#pragma once

#include <vector>
#include "SPH_2D.hpp"

int write_file(const char* filename,
           std::vector<SPH_particle> *particle_list);


#endif /* file_writer_hpp */
