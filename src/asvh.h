#ifndef  ASVH_H
#define ASVH_H

#include <string>
#include <fstream>
#include <ios>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <armadillo>

using namespace std;

using namespace arma;
using namespace std;

class asvh
{
    public:

        enum { best1bin, rand1bin, randbest1bin, best2bin, rand2bin, currrand2bin, curbest2bin,\
                      vg_best1bin, vg_rand1bin, vg_randbest1bin, vg_best2bin, vg_rand2bin, vg_currrand2bin, vg_curbest2bin
        } SCDE_type;
        enum { INV_emp_dist} quantile_type;

        asvh( string settings_file);
        ~asvh();
        asvh(const asvh& other);
        asvh& operator=(const asvh& other);

        unsigned int N_Complexes;//!< Number of complexes
        unsigned int Dim;//!< Problem dimension -- number of model parameters

        double lo_limit;//!< lower limit for initialization
        double up_limit;//!< upper limit for initialization;
        bool reject_outside_bounds_constrait;//!< Rejection the outside bounds individuals

        double CR;//!< CR parameter
        double F;//!< F parameter
        double K;//!< K parameter

        unsigned int all_n_Population;//!< The total sum of all Populations = (nPop) * N_Complexes
        unsigned int n_one_population;//!< Total number of one population members
        unsigned int all_n_Param;//!< The total number of all parameters in Population = Dim * all_n_Population

        uvec nPop;//!< Vector of population numbers, its size controlled by the number complexes (case unequal complexes)

        mat Model_param_Parents;//!< matrix(1:Dim +1 , 1:all_n_Population) -- stored parameter values of Parents models, LAST ROW FITTNES

        colvec best_model;//!< the best model (1:Dim +1)
        colvec best_model_ens;//!< the best model (1:Dim +1)
        colvec fitness;//!< stored fitness of Parents models (1all_n_Population), LAST ROW FITTNES
        unsigned int number_func_evaluation;//!< Number of function evaluations
        unsigned int Number_of_generations_in_one_complex;//!< Maximum number of allowed function evaluations in each partial complex population
        unsigned int max_number_of_shuffles;//!< Maximimu number of shuffling
        unsigned int max_function_eval;

        cube Model_Comp;//!< data in complexes 3d (1:Dim +1, 1:all_n_Population, N_Complexes), LAST ROWs FITTNES

  //      unsigned int PRoblem_Number;//!< Problem selection for calculating the fitness function
  //      unsigned int Nfunc;//!< Number of basic functions from which the fitness fucntion ois computed

        unsigned int ensemble;//!< Number of simulation runs
        unsigned int SCDE_selection;//!< SCDE_selection_0-best1bin_1-rand1bin_2-randbest1bin_3-best2bin_4-rand2bin_5-currrand2bin_6-curbest2bin

//        test_functions OBj_function;//!< Instances of the optimization problem
        bool converged;//!< Problem converged default false
        double converge_limit;//!< Precision for checking the convergence

        string directory;//!< Base directory for results output
        string path_out_file;//!< Directory for ensemble results outcome
        string SCDE_NAME;//!< Indetification of SCDE type
        stringstream EN_INFO;//!< Used parameters

        void check_convergence();//!< Method for checking the converge

        mat DE_best_one_bin(mat POpulation);//!< DE/best/1/bin
        mat DE_rand_one_bin(mat POpulation);//!< DE/rand/1/bin
        mat DE_rand_to_best_one_bin(mat POpulation);//!< DE/rand_to_best/1/bin
        mat DE_best_two_bin(mat POpulation);//!< DE/best/2/bin
        mat DE_rand_two_bin(mat POpulation);//!< DE/rand/2/bin
        mat DE_current_to_rand_2_bin(mat POpulation);//!< DE/current_to_rand/2/bin
        mat DE_current_to_best_2_bin(mat POpulation);//!< DE/current_to_best/2/bin

        void initialize_all_POpulation();//!< inicilaization of all members of Population
        void make_mat_from_Compl();//!< Combining the complexes to the parents matrix
        void sort_model_param_parents_mat();//!< Sorting the parent matrix according to fittness function
        void make_Compl_from_mat();//!< Creating the complexes from matix (P1,P2,P3,...,PN,P1,P2,P3,...,Pn,...)

        unsigned int rand_int(unsigned int n);//!< generates integer random number
        double randu_interval();//!< Random number generator (0,1)
        template<class my_vec> my_vec random_perm(unsigned int n);//! random permutations
        //uvec random_perm(unsigned int n)
        void SC_DE_best_one_bin( unsigned int ensemble_number);//!< SCDE DE/best/1/bin
        void SC_DE_rand_one_bin( unsigned int ensemble_number);//!< SCDE DE/rand/1/bin
        void SC_DE_rand_to_best_one_bin( unsigned int ensemble_number);//!< SCDE DE/rand_to_best/1/bin
        void SC_DE_best_two_bin(unsigned int ensemble_number);//!< SCDE DE/best/2/bin
        void SC_DE_rand_two_bin(unsigned int ensemble_number);//!< SCDE DE/rand/2/bin
        void SC_DE_current_to_rand_two_bin(unsigned int ensemble_number);//!< SCDE DE/current_to_rand/2/bin
        void SC_DE_current_to_best_two_bin(unsigned int ensemble_number);//!< SCDE DE/current_to_best/2/bin


        unsigned int vg_A;//! Parameter controling the geometric series of generations
        unsigned int vg_R;//! Parameter controling the geometric series of generations
        void vg_SC_DE_best_one_bin( unsigned int ensemble_number);//!< SCDE DE/best/1/bin
        void vg_SC_DE_rand_one_bin( unsigned int ensemble_number);//!< SCDE DE/rand/1/bin
        void vg_SC_DE_rand_to_best_one_bin( unsigned int ensemble_number);//!< SCDE DE/rand_to_best/1/bin
        void vg_SC_DE_best_two_bin(unsigned int ensemble_number);//!< SCDE DE/best/2/bin
        void vg_SC_DE_rand_two_bin(unsigned int ensemble_number);//!< SCDE DE/rand/2/bin
        void vg_SC_DE_current_to_rand_two_bin(unsigned int ensemble_number);//!< SCDE DE/current_to_rand/2/bin
        void vg_SC_DE_current_to_best_two_bin(unsigned int ensemble_number);//!< SCDE DE/current_to_best/2/bin

        void run_ensemble();//!< Compute ensemble run

        double fittnes(colvec X);//!< Compute fittnes
        double searched_minimum;//!< Forgiven problem searched value of global minimum

        void out_fittness(double value);//!< Writing the fittness to the file
        void out_fittness(double value, colvec X);//!<  Writing the fittnes and params to the file

        unsigned int quantile_selection;//!< Type of quantile estimates
        template<class my_vec>  my_vec quantiles_(my_vec data, my_vec probs);//!< Template for quantile estimation
        template<class my_vec> my_vec unique_(my_vec data);//!< Method for quantile estimation of 7 type Hydrman et al see cal_quantile description

        void calc_quantiles();//!< Calculation the quantiles and printing them to the file
        void calc_quantiles_fast_small();
        colvec contraint_check_stop_on_boundary(colvec MODEL);//!< checking the constraint


        /** Recesions stuff */


        enum { LinRes, TresLinRes, NonLinRes, ExpRes} res_type;

        unsigned int ResType;//!< The type of reservoirs
        unsigned int Number_of_Recessions;//!< The number of recessions
        vec ndat_inRec;
        unsigned int NRECdata;

        rowvec low_limits;
        rowvec upp_limits;

        mat recesions_data;
        mat sim_obs_recession_data;
        mat sim_obs_recession_data_ens;
        mat PARAM_ens;
        string input_data_file;

        void get_recessions_data();//!< Uploads the recesion data
        double res_outflow(colvec Par, double S);//!< the outflow from a resevoir
        double compute_ss_outflow_1Reccesion(mat Qdata, colvec Par, double So);//!< Copute the sum of sqaures of Resids the Outflow from the one recession
        double compute_fittnes(colvec Par);//!< The sum of squares
        mat get_one_recession(unsigned int Ind_rec);//!< filters form recession matrrix one recession
        colvec params_RESERVOIR(colvec Par);//!< Give the reservoirs params
        mat compute_outflow_1Reccesion(mat Qdata, colvec Par, double So);//!< Copute the Resids the Outflow from the one recession
        void out_sim_data(colvec Par, string FILET);//!< Compute and write out sim data
        mat out_sim_data_MAT(colvec Par);//!< Compute matrix of outputs



    protected:
    private:
};


#endif // ASVH_H
