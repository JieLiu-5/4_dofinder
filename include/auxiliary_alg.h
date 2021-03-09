
#include <iostream>
#include <vector>

using namespace std;

#include "bprinter/table_printer.h"
#include "PolynomialRegression.h"

#define _USE_MATH_DEFINES
#include <math.h>

#if defined(USE_BOOST_KARMA)
#include <boost/spirit/include/karma.hpp>
namespace karma = boost::spirit::karma;
#endif
using bprinter::TablePrinter;

// std::cout.precision(2);

template <typename T>
void print_vector(std::vector<T> &vector)
{
    for (T i:vector)
    {
        cout << i << " ";
    }
    cout << "\n";
}


template <typename T>
void print_vector_of_vector(std::vector<std::vector<T> > &vector_of_vector)
{
    unsigned int column_no= vector_of_vector[0].size(); 
    unsigned int row_no = vector_of_vector.size();
    
//     cout << "row_no: " << row_no << "\n";
    
    for (unsigned int k=0; k<row_no; ++k)
    {
        if(vector_of_vector[k].empty())
        {
            break;
        }else
        {
            cout << "  [" << k << "]: ";
            for (unsigned int j=0; j<column_no; ++j)
            {
                std::cout << vector_of_vector[k][j] << " ";          // << std::setprecision(2) 
            }            
        }
        cout << "\n";
    }
    cout << "\n";
}


template <typename T>
void save_vector_of_vector_to_txt(string& obj_string, std::vector<std::vector<T> > &vector_of_vector)
{
    int column_no= vector_of_vector[0].size(); 
    int row_no = vector_of_vector.size();
    
    ofstream fid;
    fid.open(obj_string+".txt");
    
    for (unsigned int k=0; k<row_no; ++k)
    {
        for (int j=0; j<column_no; ++j)
        {
            fid << vector_of_vector[k][j] << " ";
        }
        fid << "\n";
    }
    fid.close();
}


template <typename T>
string concatenate_string(std::vector<T> &vector, int id_start, int id_end)
{
    std::ostringstream streamObj;
    string string_Obj;
    
    for (int i=id_start; i<=id_end; ++i)
    {
        streamObj << vector[i];
        string_Obj=string_Obj+" "+streamObj.str();
        streamObj.str("");
    }
    return string_Obj;
}

template <typename T>
string concatenate_string(std::vector<T> &vector)
{
    std::ostringstream streamObj;
    string string_Obj;
    
    for (int i=0; i<vector.size(); ++i)
    {
        streamObj << vector[i];
        string_Obj=string_Obj+" "+streamObj.str();
        streamObj.str("");
    }
    
    return string_Obj;
}


template <typename T>
string transform_to_string(T Obj)
{
    std::ostringstream streamObj;
    streamObj << Obj;
    return streamObj.str();
}


template <typename T>
void print_table(string& string_Obj, vector<int>& vector_deg, unsigned int & deg_start, vector<vector<T>> & matrix, vector<int>& vector_n_success_var)
{
    ostringstream streamObj;
    string string_extract="";
    
    TablePrinter tp(&std::cout);
    tp.AddColumn("Degree", 10);
    if(matrix[0].size()==1)
    {
        tp.AddColumn(string_Obj, 25);
    }else
    {
        tp.AddColumn(string_Obj, 15*matrix[0].size());
    }
    tp.PrintHeader(); 
    for (unsigned int i=0; i<vector_n_success_var.size(); ++i)
    {
        for (unsigned int j=0; j<vector_n_success_var[i]; ++j)
        {   
            streamObj << matrix[i][j];
            string_extract=string_extract+" "+ streamObj.str();
            streamObj.str("");
        }
//         cout << string_extract << "\n";
        tp << to_string(vector_deg[deg_start-1+i]) << string_extract;             //   string_extract   concatenate_string(matrix[i])
        string_extract="";
    }
    tp.PrintFooter();  
}

//     cout << "\n";
//     string string_Obj;
//     string_Obj="Accuracy of tolerance";
//     print_table(string_Obj,vector_deg,deg_start,matrix_accuracy_of_tolerance,vector_n_success_var);
//     string_Obj="N_opt of tolerance";
//     print_table(string_Obj,vector_deg,deg_start,matrix_N_opt_of_tolerance,vector_n_success_var);
//     string_Obj="E_min";
//     print_table(string_Obj,vector_deg,deg_start,matrix_E_min,vector_n_all_var);
//     string_Obj="N_opt of E_min";
//     print_table(string_Obj,vector_deg,deg_start,matrix_N_opt_of_E_min,vector_n_all_var);  
//     string_Obj="critical refinement level";
//     print_table(string_Obj,vector_deg,deg_start,matrix_R_c,vector_n_all_var);
//     string_Obj="CPU of R_c";
//     print_table(string_Obj,vector_deg,deg_start,matrix_CPU_of_R_c,vector_n_all_var);
    
    





