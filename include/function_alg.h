
#include <iostream>
#include <vector>

#include <class_alg.h>
#include <color_setting.h>

using namespace std;

Dofinder::Dofinder(const int id_dim, 
                   int id_method, 
                   int id_case,
                   const int id_deg_initial, 
                   const int id_deg_last, 
                   const int var_start, 
                   const int var_end, 
                   double tolerance_target, 
                   const int max_refine):
id_dim(id_dim),
id_method(id_method),
id_case(id_case),
id_deg_initial(id_deg_initial),
id_deg_last(id_deg_last),
var_start(var_start),
var_end(var_end),
tolerance_target(tolerance_target),
max_refine(max_refine)
{
    cout << "arguments received: \n";
    
    cout << "  id_dim: " << id_dim << "\n";
    cout << "  id_method: " << id_method << "\n";
    cout << "  id_case: " << id_case << "\n";
    cout << "  id_deg_initial: " << id_deg_initial << "\n";
    cout << "  id_deg_last: " << id_deg_last << "\n";
    cout << "  var_start: " << var_start << "\n";
    cout << "  var_end: " << var_end << "\n";
    cout << "  tolerance_target: " << tolerance_target << "\n";
    cout << "  max_refine: " << max_refine << "\n";  
}


void Dofinder::step_1_setting_up()
{
    std::cout << "step 1: setting up\n";    
    
    obj_string_file_name = "data_error_"+vec_dim[id_dim]+"_"+vec_problem[id_method]+".txt";
    
    cout << "  obj_string_file_name: " << obj_string_file_name << "\n";
    
    const char * filename = obj_string_file_name.c_str();
    if( remove( filename ) != 0 )
        perror( "  error deleting file" );
    else
        puts( "  file successfully deleted" );  
    
//     vector_alpha_R_reservoir.clear();
//     vector_deg_reservoir.clear();
//     vector_var_reservoir.clear();
//     data_vector.clear();
    
    for(unsigned int i=0; i<vector_tolerance_reservior.size(); ++i)
    {
      vector_tolerance_reservior[i]=tolerance_target;
    }
    
    data_vector.resize(max_refine, std::vector<double>());
    
    var_total = var_end-var_start+1;
    deg_total = id_deg_last-id_deg_initial+1;
    
    vector_CPU_of_R_c.resize(var_total);
    vector_R_c.resize(var_total);
    
    vector_N_opt_of_tolerance.resize(var_total);
    vector_accuracy_of_tolerance.resize(var_total);
    
    vector_N_opt_of_E_min.resize(var_total);
    vector_E_min.resize(var_total);
    
    
    id_ndofs = 2;
    
    if(id_method == 0 or id_method==2)            
    {
        vector_deg_reservoir={1,2,3,4,5};            // {1,2,4,8,16};
        vector_var_reservoir=vector_var_sm;
        vector_alpha_R_reservoir=vector_alpha_R_sm;
        
        beta_R = 2.0;
        
        id_error_solu = 3;
        
        id_l2_solu = 6;
        
    }else if(id_method == 1 or id_method==3)
    {
        vector_deg_reservoir={1,2,3,4,5};            // {2,3,5,9,17};
        vector_var_reservoir=vector_var_mm;
        vector_alpha_R_reservoir=vector_alpha_R_mm;
        
        beta_R = 1.0;
        
        id_error_solu = 4;
        
        id_l2_solu = 8;
    }
    
    id_loc_CPU=id_error_solu+3;
    
    matrix_R_c.resize(deg_total, std::vector<int>());
    matrix_CPU_of_R_c.resize(deg_total, std::vector<double>());
    
    matrix_N_opt_of_E_min.resize(deg_total, std::vector<double>());
    matrix_E_min.resize(deg_total, std::vector<double>());
    
    matrix_N_opt_of_tolerance.resize(deg_total, std::vector<double>());
    matrix_accuracy_of_tolerance.resize(deg_total, std::vector<double>());

    if (id_dim==0)
    {
        R_min_for_l2_norm = 2;
    }else if (id_dim==1)
    {
        R_min_for_l2_norm = 2;
    }
    
    criterion_of_l2_norm = 0.01;
    
    for(unsigned int i=0; i<deg_total; ++i)
    {
        vector_deg_of_interest.push_back(vector_deg_reservoir[id_deg_initial+i]);
    }
    
    for(unsigned int i=0; i<var_total; ++i)
    {
        vector_var_of_interest.push_back(vector_var_reservoir[var_start+i]);
        vector_tolerance_of_interest.push_back(vector_tolerance_reservior[var_start+i]);
        vector_alpha_R_of_interest.push_back(vector_alpha_R_reservoir[var_start+i]);
    }
    
    print_input_info_general();
    
}

void Dofinder::print_input_info_general()
{   
    
    cout << "  Input information\n";   
    
    TablePrinter tp_custom(&std::cout);
    tp_custom.AddColumn("Custom item", 25);
    tp_custom.AddColumn("Value", 20);
    tp_custom.PrintHeader();
    tp_custom << "Dimension" << vec_dim[id_dim];
    tp_custom << "Problem" << vec_problem[id_method];
    tp_custom << "Variables of interest" << concatenate_string(vector_var_of_interest);
    tp_custom << "associated tolerances" << concatenate_string(vector_tolerance_of_interest);
    tp_custom << "element degrees of interest" << concatenate_string(vector_deg_of_interest);
    tp_custom.PrintFooter();
    
    TablePrinter tp_dynamic(&std::cout);
    tp_dynamic.AddColumn("Dynamic item", 25);
    tp_dynamic.AddColumn("Value", 20);
    tp_dynamic.PrintHeader(); 
    tp_dynamic << "beta_R" << transform_to_string(beta_R);
    tp_dynamic << "alpha_R's" << concatenate_string(vector_alpha_R_of_interest);    
    tp_dynamic.PrintFooter();
    
    TablePrinter tp_default(&std::cout);
    tp_default.AddColumn("Default item", 25);
    tp_default.AddColumn("Value", 20);  
    tp_default.PrintHeader();
    tp_default << "R_min for l2 norm" << R_min_for_l2_norm;
    tp_default << "criterion for l2 norm" << transform_to_string(criterion_of_l2_norm);
    tp_default << "space_before_error_solu" << id_error_solu;
    tp_default.PrintFooter();   
}


void Dofinder::updating_data_vector(std::vector<std::vector<double> > &data_vector)
{
    
    ifstream input_file;
    std::string data_line;
    std::string::size_type sz;
    
    const char * filename = obj_string_file_name.c_str();
    input_file.open(filename);
    
    if (!input_file)
    {
        cout << "Unable to open file\n";
        exit(1);
    }

    vector<string> data_per_line;
    string element_of_data_per_line;
    
    unsigned int id_row=0;
    
    while (getline(input_file, data_line))                                  // note: this only applies to the situation that there is only one line in a .txt file
    {
        
        string input(data_line);
        for (string::iterator it = input.begin(); it!=input.end(); ++it)
        {   
            if (isspace(*it))                                           // *it represents one character
            {
//                 cout << "element_of_data_per_line: " << element_of_data_per_line << "\n";
                data_per_line.push_back(element_of_data_per_line);                                       // attaching an integral part of data per line
                element_of_data_per_line="";
            }else
            {
                element_of_data_per_line = element_of_data_per_line+*it;
            }
        }
        
        data_per_line.push_back(element_of_data_per_line);
        
        data_vector[real_refine-1].clear();
        for (unsigned int i = 0; i < data_per_line.size()-1; ++i)
        {
            data_vector[real_refine-1].push_back(stod(data_per_line[i], &sz));
        }
    }
    input_file.close();
}


void Dofinder::performing_real_computation(unsigned int degree)
{
    stringstream ss;

    ss << "../external/"+vec_dim[id_dim]+"/step-super "
        << id_method << " "
        << id_case << " "
        << 1 << " "
        << 2 << " "
        << 0 << " "
        << 1 << " "
        << 1 << " "
        << 0 << " "
        << degree << " "
        << real_refine;

    system(ss.str().c_str());
}


void Dofinder::initializing_data_vector(unsigned int degree_being_used, unsigned int R_min)
{
    cout << "  Initializing (degree: " << degree_being_used << ", R_min: " << R_min << ")\n";
    
    data_vector.clear();
    data_vector.resize(max_refine, std::vector<double>());
    
    real_refine = 0;
    
//     cout << "  #REFINE: 1 --> " << R_min << " ";
    for (unsigned int i = 0; i < R_min; ++i)
    {
        ++real_refine;
        performing_real_computation(degree_being_used);
        updating_data_vector(data_vector);
    }
    
    cout << "  data_vector after initialization:\n";
    print_vector_of_vector(data_vector);
}

void Dofinder::retrieving_l2_solu()
{
    vector_l2_solu.clear();
    
    for(unsigned int j=0; j<real_refine; ++j)
    {
        vector_l2_solu.push_back(data_vector[j][id_l2_solu]);
    }
    
//     cout << "l2 norm of the solution\n";
//     print_vector(vector_l2_solu);
}



void Dofinder::h_refinement_for_l2_norm()
{
    performing_real_computation(vector_deg_reservoir[id_deg_initial]);

    updating_data_vector(data_vector);
    
    vector_l2_solu.push_back(data_vector[real_refine-1][id_l2_solu]);
    
    l2_solu=vector_l2_solu.back();          

    quality_l2_solu = fabs(vector_l2_solu[real_refine-1] -vector_l2_solu[real_refine-2])/vector_l2_solu[real_refine-1];
}



void Dofinder::step_2_seeking_l2_norm_and_determining_alpha_R()
{
    cout << "step 2: seeking l2 norm\n";    
    
    initializing_data_vector(vector_deg_reservoir[id_deg_initial], R_min_for_l2_norm);
    retrieving_l2_solu();    
    
#if 1

    quality_l2_solu = fabs(vector_l2_solu[real_refine-1] -vector_l2_solu[real_refine-2])/vector_l2_solu[real_refine-1];
            
    for (unsigned int i = 0; i<1; ++i)
    {
        if (quality_l2_solu > criterion_of_l2_norm)			//
        {
//             cout << BOLDRED << " not reached" << "\n" << RESET;
            ++real_refine;
            h_refinement_for_l2_norm();
        }else
        {
//             cout << BOLDGREEN << "  reached" << "\n\n" << RESET;
            l2_solu = vector_l2_solu.back();
            
//             if(id_method == 1)
//             {
//             l2_grad = data_vector[real_refine-1][id_l2_solu+1];
//             }
            break;
        }  
    }
    
    cout << "  l2_solu: " << l2_solu << " (last refinement: "<< real_refine << ")\n";
    
    if(id_method==0)
    {
        for (int i = 0; i<=var_total; ++i)
        {
            vector_alpha_R_of_interest[i]=vector_alpha_R_of_interest[i]*l2_solu;
        }
    }else if(id_method==1)
    {
        cout << "  l2_grad: " << l2_grad << "\n"; 
        
        for (int i = var_start; i<=var_end; ++i)
        {
            if(i==0 or i==1)
            {
                vector_alpha_R_of_interest[i-var_start]=vector_alpha_R_of_interest[i-var_start]*l2_solu;
            }else if(i==2)
            {
                vector_alpha_R_of_interest[i-var_start]=vector_alpha_R_of_interest[i-var_start]*l2_grad;
            }
        }        
    }
    
    cout << "  alpha_R's: ";
    print_vector(vector_alpha_R_of_interest);
    
#endif
    
}

void Dofinder::searching_over_mixed_FEM()
{
    id_method=1;
    
    cout << "  Transiting to ";
    cout << vec_problem[id_method] << endl;
    cout << "\n\n\n\n";
}


void Dofinder::determining_R_min_for_pred_and_c_relax()
{
    
    if (degree_being_used<6)
    {
        if (id_dim==0)
        {
            R_min_for_prediction = 2;                    // 9-degree_being_used
        }else if (id_dim==1)
        {
            R_min_for_prediction = 4;
        }
    }else
    {
        R_min_for_prediction = 3;
    }
    
    c_relax_for_convergence_max = 1.1;         // 1.05, this parameter makes a difference for P1/P0 and P2/P1 elements, Sep. 1, 2019
    
    if(id_method == 0)
    {
        if (degree_being_used<4)
        {
            c_relax_for_convergence_min = 0.9;
        }else
        {
            c_relax_for_convergence_min = 0.7;      // 0.7    0.5,  to obtain an acceptable order of convergence, the relaxation coefficient is necessary, Aug. 15, 2019
        }
        
    }else if(id_method == 1)
    {
//         if(degree_being_used==1)
//         {
//             R_min = 8;
//         }
        
        if(degree_being_used<4)
        {
            c_relax_for_convergence_min = 0.9;         // 0.9
        }else if ((degree_being_used >= 4) && (degree_being_used < 10))
        {
            c_relax_for_convergence_min = 0.7;
        }else
        {
            c_relax_for_convergence_min = 0.5;
        }
    }
    
        cout << "  R_min for prediction: " << R_min_for_prediction << "\n";
//         cout << "  c_relax: (" << c_relax_for_convergence_min << ", " << c_relax_for_convergence_max << ")\n";    
}
        

void Dofinder::determining_beta_T_alpha_R_etc_for_one_var()
{
    if(id_method == 0 or id_method==2)
    {
        beta_T=degree_being_used+1-id_var_being_investigated; 
    }else if(id_method == 1 or id_method==3)
    {
        switch(id_var_being_investigated)                                
        {
            case 0:
                beta_T = degree_being_used;
                break;
            case 1:
                beta_T = degree_being_used+1;
                id_ndofs = 3;
                break;
            case 2:
                beta_T = degree_being_used;
                id_ndofs = 3;
                break;
        }        
    }
    
//     if(id_dim==1)
//     {
//         beta_T = beta_T/2;
//     }

    cout << "beta_T: " << beta_T << "\n";
    
    beta_T_bound_max = beta_T * c_relax_for_convergence_max;
    beta_T_bound_min = beta_T * c_relax_for_convergence_min; 
    
    alpha_R=vector_alpha_R_of_interest[id_var_being_investigated-var_start];
}


void Dofinder::step_3_initializing_data_of_one_degree(unsigned int input_degree)
{
    
    degree_being_used = input_degree;
    cout << BOLDCYAN << "  degree: " << degree_being_used << "\n" << RESET;
    n_degrees_used += 1;                                // only for determing the number of element degrees used for prediction
    
    real_refine=0;
    id_success_one_deg=1;   
    
    vector_accuracy_of_tolerance.clear();    
    vector_N_opt_of_tolerance.clear();
    vector_E_min.clear();  
    vector_N_opt_of_E_min.clear();
    vector_R_c.clear();
    vector_CPU_of_R_c.clear();

    vector_accuracy_of_tolerance.resize(var_total);
    vector_N_opt_of_tolerance.resize(var_total);
    vector_E_min.resize(var_total);   
    vector_N_opt_of_E_min.resize(var_total);
    vector_R_c.resize(var_total);
    vector_CPU_of_R_c.resize(var_total);   
    
    determining_R_min_for_pred_and_c_relax();
    initializing_data_vector(degree_being_used, R_min_for_prediction);
}

void Dofinder::step_4_looping_over_all_var()
{
    for (int i = var_start; i<=var_end; ++i)
    {
        id_var_being_investigated = i;
        cout << BOLDYELLOW << "\n  \"" << vector_var_reservoir[id_var_being_investigated] << "\"\n" << RESET;
            
        step_5_determing_n_opt_for_one_var();
        
        /*if(id_success_one_deg==0)
        {
            break;              // terminate the calculation for the next var
        } */
    }     
}


void Dofinder::print_results_of_one_degree(unsigned int id_degree)
{
    cout << "\n";
    
    TablePrinter tp_one_degree(&std::cout);
    tp_one_degree.AddColumn("Item", 25);
    tp_one_degree.AddColumn("Value", 13*var_total);

    tp_one_degree.PrintHeader();
    tp_one_degree << "failed variable" << vector_failed_var[id_degree];
    tp_one_degree << "accuracy of tolerance" << concatenate_string(matrix_accuracy_of_tolerance[id_degree]);
    tp_one_degree << "N_opt of tolerance" << concatenate_string(matrix_N_opt_of_tolerance[id_degree]);
    tp_one_degree << "E_min" << concatenate_string(matrix_E_min[id_degree]);
    tp_one_degree << "N_opt of E_min" << concatenate_string(matrix_N_opt_of_E_min[id_degree]);
    tp_one_degree << "R_c" << concatenate_string(matrix_R_c[id_degree]);
    tp_one_degree << "CPU of R_c" << concatenate_string(matrix_CPU_of_R_c[id_degree]);
    tp_one_degree.PrintFooter();     
    
    cout << "\n";
    
}

void Dofinder::store_data_of_one_degree(unsigned int id_degree)
{
    matrix_R_c[id_degree] = vector_R_c;
    matrix_CPU_of_R_c[id_degree] = vector_CPU_of_R_c;
    
    matrix_N_opt_of_tolerance[id_degree] = vector_N_opt_of_tolerance;
    matrix_accuracy_of_tolerance[id_degree] = vector_accuracy_of_tolerance;
    
    matrix_N_opt_of_E_min[id_degree] = vector_N_opt_of_E_min;
    matrix_E_min[id_degree] = vector_E_min;       
    
//     print_results_of_one_degree(id_degree);
}

void Dofinder::step_5_determing_n_opt_for_one_var()
{
    if(id_method==0 && degree_being_used==1 && id_var_being_investigated==2)
    {
        cout << "  nonexistent variable\n";
        N_opt_of_tolerance = 1e100;
        E_R_of_N_opt_of_tolerance = 1e100;
        N_opt_of_E_min = 1e100;
        E_min = 1e100;
//                 break;
    }else
    {
        cout << "\n";
        
        tolerance_being_investigated=vector_tolerance_of_interest[id_var_being_investigated-var_start];
            
        id_loc_err=id_error_solu+id_var_being_investigated;       
        determining_beta_T_alpha_R_etc_for_one_var();                   // parameters completed
        
        cout << left << setw(15) << "  tolerance_being_investigated" << tolerance_being_investigated << "\n";
        cout << setw(15) << "  alpha_R" << alpha_R << "\n";
        cout << setw(15) << "  beta_T" << beta_T_bound_min << "~" << beta_T_bound_max << "\n";   
        
    //     if(id_method == 1 && id_var_being_investigated>1)
    //     {
    //         cout << "  l2_solu: " << l2_solu << ", l2_grad: " << l2_grad << "\n";               
    //         initializing_data_vector(degree_being_used, R_min_for_prediction);
    //     }        
        
        for (int id_refine=real_refine; id_refine<max_refine; ++id_refine)
        {
            E_second_last = data_vector[id_refine-2][id_loc_err];
            E_last=data_vector[id_refine-1][id_loc_err];
            
            E_R = alpha_R*pow(data_vector[id_refine-1][2], beta_R);
            cout << "  E_R: " << E_R << "\n";            
        
            if (E_last > E_R)            
            {
                cout << "  E_last: " << E_last << " > " << "E_R" << ", ";
                
                if (E_last > tolerance_being_investigated)
                {
                    cout << "> " << "tolerance_being_investigated" << "\n";
                    
                    tilde_P = log2(E_second_last/E_last);
                    std::cout << "  numerical beta_T: " << tilde_P << "";
                    
                    if (tilde_P < beta_T_bound_min)
                    {
                        cout << " < " << "beta_T_bound_min" << "\n\n";
                //                     if (E_last < E_second_last)
                //                     {
                //                         cout << "  E_last still decreasing" << endl;
                        h_refinement_for_convergence();
                //                     }else
                //                     {
                //                         vector_R_c[id_var_being_investigated]=id_refine-1;
                //                         break;                            
                //                     }

                    }else
                    {
                        cout << " >= " << "beta_T_bound_min" << ", ";
                        if (tilde_P < beta_T_bound_max)
                        {
                            cout << BOLDGREEN << "< " << "beta_T_bound_max\n\n" << RESET;
                            
                            alpha_T= E_last/(pow(data_vector[real_refine-1][2],-beta_T));           // 
                            std::cout << "  alpha_T: " << alpha_T << std::endl;                            
                            
                            collecting_data_for_R_c(real_refine);                                     // 
                            
                            predicting_n_opt_using_formular();
                            
                            break;                        
                        }else
                        {
                            cout << ">= " << "beta_T_bound_max\n\n";
                            cout << "  #REFINE: ";
                            cout << id_refine+1 << " --> ";
                            
                            h_refinement_for_convergence();
                            
                            cout << "\n";
                        }
                    }                       
                    
                    
                }else           // this means the tolerance can be reached by only seeking the desired l2 norm of the solution
                {
                    cout << BOLDGREEN << "< tolerance_being_investigated" << ", \n" << RESET;
    //                     cout << "  optimal refinement level: " <<   [real_refine-1][0] << std::endl;

                    N_opt_of_tolerance = data_vector[real_refine-1][id_ndofs];
                    E_R_of_N_opt_of_tolerance=alpha_R*pow(N_opt_of_tolerance, beta_R);
                    
                    E_min = tolerance_being_investigated;                          // 
                    
                    collecting_data_for_R_c(real_refine);     
                
                    break;
                }
            }else
            {
                cout << "  E_last: " << E_last << ", E_R: " << E_R;
                cout << " --> h_refinement stops because E_R > E_last" << std::endl;

                N_opt_of_E_min = data_vector[real_refine-2][id_ndofs];
                E_min = data_vector[real_refine-2][id_loc_err];   
                
                collecting_data_for_R_c(real_refine);                    
    
                break;
            }
            
        }        
    }

    
    accuracy_of_tolerance = tolerance_being_investigated+E_R_of_N_opt_of_tolerance;
    cout << "  accuracy of tolerance: " << accuracy_of_tolerance << " ";
    if (E_R_of_N_opt_of_tolerance<0.5*tolerance_being_investigated)
    {
        cout << BOLDGREEN << "successful" << "\n" << RESET;
    }else
    {
        cout << BOLDRED << "failed" << "\n" << RESET;        
        id_success_one_deg=0;
    }
    
    vector_N_opt_of_tolerance[id_var_being_investigated-var_start]=N_opt_of_tolerance;
    vector_accuracy_of_tolerance[id_var_being_investigated-var_start]=accuracy_of_tolerance;
    
    vector_N_opt_of_E_min[id_var_being_investigated-var_start]=N_opt_of_E_min;
    vector_E_min[id_var_being_investigated-var_start]=E_min;   
        
    
//     cout << "\n";
//     cout << "  N_opt of tolerance: " << N_opt_of_tolerance << ", ";   
//     cout  << "E_R of N_opt of tolerance: " << E_R_of_N_opt_of_tolerance << " ";
    
//     cout << "  N_opt of prediction: " << N_opt_of_E_min << ", ";
//     cout << "E_min: " << E_min << "\n\n"; 
    
    N_opt_of_E_min=0;
    E_min=0;       
}


void Dofinder::collecting_data_for_R_c(int id_refine)
{
    vector_R_c[id_var_being_investigated-var_start]=id_refine;
    
    cout << "\n";
    
    cout << "  CPU time to reach R_c: ";
    for(unsigned int i = 0; i<id_refine; ++i)
    {
        cout << data_vector[i][id_loc_CPU] << ", ";
        vector_CPU_of_R_c[id_var_being_investigated-var_start] += data_vector[i][id_loc_CPU];
    }
    cout << "\n";        
    cout << "  total CPU time: " << vector_CPU_of_R_c[id_var_being_investigated] << "\n";
}


void Dofinder::h_refinement_for_convergence()
{
    real_refine++;         // real_refine does not increases when the analytical convergence rate is reached or the tolerance is satisfied
    performing_real_computation(degree_being_used);
    updating_data_vector(data_vector);
}

void Dofinder::predicting_n_opt_using_formular()
{
/*    if(id_dim==1)                                // when using formulas to predict N_opt, the slope halves
    {
        beta_T = beta_T/2;
    }  */  
        
    N_opt_of_tolerance = pow(alpha_T/(1.0*tolerance_being_investigated), 1/beta_T);
    E_R_of_N_opt_of_tolerance=alpha_R*pow(N_opt_of_tolerance, beta_R);
    
    
    N_opt_of_E_min = pow(alpha_T*beta_T/(alpha_R*beta_R), 1.0/(beta_R+beta_T));
    E_min = alpha_T*pow(N_opt_of_E_min,-beta_T) + alpha_R*pow(N_opt_of_E_min,beta_R);

}

void Dofinder::success_indicator_for_all_degrees_used()
{
    int is_success_this_degree;
    
    for (unsigned int i = 0; i < min(n_degrees_used,deg_total); ++i)
    {
        is_success_this_degree = 1;
        
        for (unsigned int j = 0; j < var_total; ++j)
        {
//             cout << vector_var_reservoir[var_start+j] << " ";
            
            if (matrix_E_min[i][j]<vector_tolerance_reservior[var_start+j])
            {
//                 n_success_var++;
                
//                 cout << BOLDGREEN << "satisfied \n" << RESET;
                vector_n_success_var[i]++;
            }
            else
            {
//                 cout << BOLDRED << "not satisfied " << RESET;
                vector_failed_var[i]= vector_failed_var[i] + " " + vector_var_reservoir[var_start+j];
                is_success_this_degree = 0;
//                 break;
            }
        }
        
//         cout << "\n";
//         cout << "degree_being_used " << vector_deg_reservoir[id_deg_initial+i] << ": ";
//         
//         if(is_success_this_degree==1)
//         {
//             cout << BOLDGREEN << "successful\n" << RESET;
//         }
//         else
//         {
//             cout << BOLDRED << "failed\n\n" << RESET;
//         }
    }
}

void Dofinder::printing_final_results()
{
    vector_n_success_var.resize(n_degrees_used);
    vector_n_all_var.resize(n_degrees_used);
    vector_failed_var.resize(n_degrees_used);
    
    for (unsigned int i=0;i<deg_total;++i)
    {
        vector_n_all_var[i]=var_total;
    }    
    
    success_indicator_for_all_degrees_used();
    
    if(id_success_one_deg==1)
    {
        cout << BOLDGREEN << "  current setting available\n" << RESET;
        
        
        TablePrinter tp_report(&std::cout);
        tp_report.AddColumn("Item", 25);
        tp_report.AddColumn("Value", 20);  
        tp_report.PrintHeader();
        tp_report << "available method" << vec_problem[id_method];
        tp_report << "available minimal p" << degree_being_used;
        tp_report.PrintFooter(); 
        
//         cout << "  number of successful var of each degree_being_used\n";
//         for (unsigned int i=0; i<n_degrees_used; ++i)
//         {
//             cout << "  degree_being_used " << vector_deg_reservoir[id_deg_initial+i] << ": " << vector_n_success_var[i] << "\n";
//         }
    //     print_vector(vector_n_success_var);
    //     print_vector(vector_n_all_var);     
        
    }else if(id_success_one_deg==0)
    {
        id_success_FEM=0;
        
        cout << BOLDRED << "  current setting unavailable\n" << RESET;
        
//         if(id_method==0)
//         {
//             searching_over_mixed_FEM();
//         }
    }
    
    cout << "\n";
    cout << "  using " << n_degrees_used << " (of " << deg_total << ")" << " degree(s)\n";
    for(unsigned int id_degree=0;id_degree<n_degrees_used;++id_degree)
    {
        cout << "  degree_being_used: " << vector_deg_reservoir[id_deg_initial+id_degree];
        print_results_of_one_degree(id_degree);
    }    
    cout << "\n\n";
}


int Dofinder::run()
{  
    step_1_setting_up();
    
    step_2_seeking_l2_norm_and_determining_alpha_R();
    
#if 1
    cout << "\nstep 3: prediction\n";
    for (unsigned int id_degree = 0; id_degree < deg_total; ++id_degree)
    {
        step_3_initializing_data_of_one_degree(vector_deg_reservoir[id_deg_initial+id_degree]);     
        
        step_4_looping_over_all_var();

        store_data_of_one_degree(id_degree);
        if(id_success_one_deg==1)
        {
            break;
        }
    }
    
    
    cout << "\nstep 4: output" << endl;
    printing_final_results();
    
    obj_string_assistant = "data_E_min";
    save_vector_of_vector_to_txt(obj_string_assistant, matrix_E_min);
    obj_string_assistant = "data_N_opt_of_E_min";
    save_vector_of_vector_to_txt(obj_string_assistant, matrix_N_opt_of_E_min);
    obj_string_assistant = "data_CPU_of_R_c";
    save_vector_of_vector_to_txt(obj_string_assistant, matrix_CPU_of_R_c);
    obj_string_assistant = "data_R_c";
    save_vector_of_vector_to_txt(obj_string_assistant, matrix_R_c);
    
#endif
    return id_success_FEM;
    
}

