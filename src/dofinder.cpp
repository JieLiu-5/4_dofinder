#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <math.h>
#include <iomanip>

#include "auxiliary_alg.h"
#include "function_alg.h"

using namespace std;

int main(int argc, char *argv[])
{
    int id_dim;
    int id_method;
    int id_case;
    int deg_start, deg_end;
    int var_start, var_end;
    
    double tolerance_target = 1e-12;
    int max_refine = 25;
    

    if (argc != 10)
    {
        std::cout<<"usage: "<< argv[0] <<" <id_dim> <id_problem> <id_case> <deg_start> <deg_end> <var_start> <var_end> <tolerance> <max_refine>\n";
        exit(EXIT_FAILURE);
    } else
    {
        id_dim = atoi(argv[1]);
        id_method = atoi(argv[2]);
        id_case = atoi(argv[3]);
        deg_start = atoi(argv[4]);
        deg_end = atoi(argv[5]);
        var_start = atoi(argv[6]);
        var_end = atoi(argv[7]);
        tolerance_target = atof(argv[8]);
        max_refine = atoi(argv[9]);
    }
            
    Dofinder obj_dofinder(id_dim, id_method, id_case, deg_start, deg_end, var_start, var_end, tolerance_target, max_refine);
    
    int id_success_FEM = obj_dofinder.run();
    
//     if(id_success_FEM==0)                       // we try the mixed FEM once if the standard FEM does not work
//     {
//         obj_dofinder.run();
//     }
    
    return 0;

}
