
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Dofinder
{
public:
    Dofinder(const int id_dim, 
             int id_method, 
             int id_case,
             const int id_deg_initial, 
             const int id_deg_last, 
             const int var_start, 
             const int var_end, 
             double tolerance_target, 
             const int max_refine);
    
    int run();
    
//     ~Dofinder();
        
private:
    
    
    std::vector<std::vector<double> > data_vector;
    void updating_data_vector(std::vector<std::vector<double> > &data_vector);
    void initializing_data_vector(unsigned int degree_being_used, unsigned int R_min);
    
    void performing_real_computation(unsigned int degree_being_used);
    
    
    
    void step_1_setting_up();
    void print_input_info_general();
    
    
    void step_2_seeking_l2_norm_and_determining_alpha_R();
    int R_min_for_l2_norm;
    void retrieving_l2_solu();
    void h_refinement_for_l2_norm();
    vector<double> vector_l2_solu;
    double l2_solu, l2_grad;
    double quality_l2_solu, criterion_of_l2_norm;
    std::vector<double> vector_alpha_R_of_interest;
        
    
    void step_3_initializing_data_of_one_degree(unsigned int input_degree);
    int R_min_for_prediction;
    void determining_R_min_for_pred_and_c_relax();                               // updating values based on degree
    void collecting_data_for_R_c(int refine);
    void store_data_of_one_degree(unsigned int id_degree);
    void print_results_of_one_degree(unsigned int id_degree);

    void step_4_looping_over_all_var();
    vector<std::string> vector_var_of_interest;
    void determining_beta_T_alpha_R_etc_for_one_var();                                       // updating values based on degree and var              
    
    void step_5_determing_n_opt_for_one_var();
    int id_var_being_investigated;
    double beta_T, beta_T_bound_min, beta_T_bound_max;    
    double c_relax_for_convergence_min, c_relax_for_convergence_max;
    void h_refinement_for_convergence();
    double alpha_T;
    void predicting_n_opt_using_formular();
    
    
    void printing_final_results();
    void success_indicator_for_all_degrees_used();    
    
    void searching_over_mixed_FEM();
    
    
    
    const int id_dim=0;                   // '0' for 1D, '1' for '2D'
                                          // this parameter is only for this program not for the deal.ii program
    vector<string> vec_dim = {"1d", "2d"};

    int id_method = 0; 
    vector<std::string> vec_problem={"sm_real", "mm_real", "sm_complex", "mm_complex"};
    
    const unsigned int id_case;
    
    vector<int> vector_deg_reservoir;
    vector<int> vector_deg_of_interest;
    int degree_being_used;
    unsigned int id_deg_initial = 0, id_deg_last;
    unsigned int deg_total = 1;
    
    

    vector<std::string> vector_var_reservoir;
    vector<std::string> vector_var_sm={"u","u_x","u_xx"};
    vector<std::string> vector_var_mm={"u","v","v_x"};
    
    
    unsigned int var_start = 0, var_end;    
    unsigned int var_total = 1;
    
    int id_error_solu;
    int id_ndofs, id_loc_err, id_loc_CPU, id_l2_solu;
    
    const int max_refine;                           // denoting the maximum rows of data_vector
    int real_refine = 0;
    
    string obj_string_file_name, obj_string_assistant;
    
    double alpha_R, beta_R;
    double E_second_last, E_last, tilde_P, E_R; 
    
    
    double N_opt_of_tolerance, E_R_of_N_opt_of_tolerance, accuracy_of_tolerance;
    double N_opt_of_E_min, E_min;

    std::vector<double> vector_alpha_R_reservoir;
    std::vector<double> vector_alpha_R_sm={2e-17,5e-17,1e-15};
    std::vector<double> vector_alpha_R_mm={2e-17,2e-16*0.19,1e-15};                      // we calculate ||d||_2 manually
    
    double tolerance_target;
    double tolerance_being_investigated;
    vector<double> vector_tolerance_reservior={1e-2,1e-2,1e-2};
    vector<double> vector_tolerance_of_interest;    
    
    vector<int> vector_R_c;
    vector<double> vector_CPU_of_R_c;
    
    vector<double> vector_N_opt_of_E_min, vector_E_min;
    vector<double> vector_N_opt_of_tolerance, vector_accuracy_of_tolerance;
    
    vector<vector<int>> matrix_R_c;
    vector<vector<double>> matrix_CPU_of_R_c;
    
    vector<vector<double>> matrix_N_opt_of_E_min;
    vector<vector<double>> matrix_E_min;
    
    vector<vector<double>> matrix_N_opt_of_tolerance;
    vector<vector<double>> matrix_accuracy_of_tolerance;
    
    int id_success_one_deg=0, id_success_FEM=0;
    
    unsigned int n_degrees_used=0;
    
    vector<int> vector_n_success_var;
    vector<int> vector_n_all_var;
    vector<string> vector_failed_var;
    
};




