#include "asvh.h"

asvh::asvh(string settings_file)
{
    string trash;

    ifstream asvh_settings(settings_file.c_str());
    if(!asvh_settings) { // osetreni jestli je soubor mozne otevrit
      cout << "\nUnable to open file :"<< settings_file.c_str() <<" with initial settings for computations.";
      exit(EXIT_FAILURE);
    }
    asvh_settings >> scientific;
    asvh_settings >> trash >> trash >> N_Complexes >> trash >> n_one_population;
//   asvh_settings >> trash >> lo_limit >> trash >> up_limit;
//    asvh_settings >> trash >> PRoblem_Number;
//    if((PRoblem_Number >=0)&&(PRoblem_Number <=10)) Nfunc = 1;
    asvh_settings >> trash >> Number_of_generations_in_one_complex;
    asvh_settings >> trash >> max_number_of_shuffles;
    asvh_settings >> trash >> max_function_eval;
    asvh_settings >> trash >> F >> CR >> K;
    asvh_settings >> trash >> converge_limit;
    asvh_settings >> trash >> ensemble;
    asvh_settings >> trash >> SCDE_selection;
    asvh_settings >> trash >> trash;
    asvh_settings >> trash >> vg_A >> vg_R;
    asvh_settings >> trash >> quantile_selection;
    asvh_settings >> trash >> reject_outside_bounds_constrait;
    asvh_settings >> trash >> directory;
    asvh_settings >> trash >> ResType;
    asvh_settings >> trash >> Number_of_Recessions;
//    asvh_settings >> trash >> path_out_file;

    switch (ResType){
        case LinRes:
            Dim = 1 + Number_of_Recessions;
          break;
        case TresLinRes:
            Dim = 1 + Number_of_Recessions;
          break;
        case NonLinRes:
            Dim = 2 + Number_of_Recessions;
         break;
       case ExpRes:
           Dim = 2 + Number_of_Recessions;
         break;
 //       default:
 //        break;
    }

    if(Dim<=1) { // osetreni jestli je soubor mozne otevrit
      cout << "\nSmall number of parameters  in "<< settings_file.c_str() <<".";
      exit(EXIT_FAILURE);
    }
    low_limits.set_size(Dim);
    upp_limits.set_size(Dim);
    low_limits.fill(-99999);
    upp_limits.fill(99999);

    asvh_settings >> trash;
    double help_var=0.0;
      switch (ResType){
        case LinRes:
            asvh_settings >> help_var;
            low_limits.fill(help_var);
            asvh_settings >> low_limits(Dim-1);
          break;
        case TresLinRes:
            asvh_settings >> help_var;
            low_limits.fill(help_var);
            asvh_settings >> low_limits(Dim-1);
          break;
        case NonLinRes:
            asvh_settings >> help_var;
            low_limits.fill(help_var);
            asvh_settings >> low_limits(Dim-2);
            asvh_settings >> low_limits(Dim-1);
         break;
       case ExpRes:
            asvh_settings >> help_var;
            low_limits.fill(help_var);
            asvh_settings >> low_limits(Dim-2);
            asvh_settings >> low_limits(Dim-1);
        break;
 //       default:
 //        break;
    }


    asvh_settings >> trash;
   // cout << endl << low_limits << endl;
     switch (ResType){
        case LinRes:
            asvh_settings >> help_var;
            upp_limits.fill(help_var);
            asvh_settings >> upp_limits(Dim-1);
          break;
        case TresLinRes:
            asvh_settings >> help_var;
            upp_limits.fill(help_var);
            asvh_settings >> upp_limits(Dim-1);
          break;
        case NonLinRes:
            asvh_settings >> help_var;
            upp_limits.fill(help_var);
            asvh_settings >> upp_limits(Dim-2);
            asvh_settings >> upp_limits(Dim-1);
         break;
       case ExpRes:
            asvh_settings >> help_var;
            upp_limits.fill(help_var);
            asvh_settings >> upp_limits(Dim-2);
            asvh_settings >> upp_limits(Dim-1);
        break;
 //       default:
 //        break;
    }

//
//
//    for (unsigned int i=0; i<1 ;i++ ){
//            asvh_settings  >> upp_limits(i);
//      }
   // cout << endl << upp_limits << endl;

    asvh_settings >> trash >> input_data_file;
    asvh_settings.close();



    all_n_Population = N_Complexes * n_one_population;
    all_n_Param = Dim *all_n_Population;

    cout << endl << "**** REVIEW of ASVH settings ****"<< endl;
    cout << scientific;
    cout << "\nDE initialization with following settings from file: " << settings_file.c_str();
    cout << "\nNumber of complexes: " << N_Complexes << "\nNumber of members in one Population: " << n_one_population << endl;
    cout << "The total population: " << N_Complexes * n_one_population << "\nDimension of the problem: " << Dim  <<endl;
    cout << "Total number of real parameters in population: " << all_n_Param;// << "\nLower bound for parameter values: " << lo_limit << endl;
    //cout << "Upper bound for parameter limit: " << up_limit<< endl;
 //   cout << "Problem optimization 0 -> sphere etc: " << PRoblem_Number << endl;
    cout << "\nFuntion evaluations in one complex is " <<    Number_of_generations_in_one_complex << endl;
    cout << "Maximum allowed shuffling  is " <<    max_number_of_shuffles << endl;
    cout << "Converge limit is " << converge_limit << endl;
    cout << "Max_func_eval: " << max_function_eval << endl;
    cout << "Read parameter values F - CR - K: " << F << " "<< CR << " " << K << endl;
    cout << "The ensemble " << ensemble << endl;
    cout << "Algorithm_selection: " << SCDE_selection << endl;
    cout << "SCDE_selection_0-best1bin_1-rand1bin_2-randbest1bin_3-best2bin_4-rand2bin_5-currrand2bin_6-curbest2bin " << endl;
    cout <<  "vg_SCDE_selection_7-best1bin_8-rand1bin_9-randbest1bin_10-best2bin_11-rand2bin_12-currrand2bin_13-curbest2bin" << endl;
    cout << "vg_Parameters_vgA_vgR: " << vg_A << " " << vg_R << endl;
    cout << "Quantile estimates according to Hyndman 0-Inverse of Empirical function: " << quantile_selection << endl;
    cout << "Constrains_O-noconstr_1-inside-bounds: " << reject_outside_bounds_constrait << endl;
    cout << "Directory_for_output: " << directory.c_str() << endl;
    cout <<  "The_Resevoir_Type_0-kLR_1-TresLR_2-NR_3-EXP: " << ResType << endl;
    cout <<  "The_number_of_Recessions:  " << Number_of_Recessions << endl;
//    cout <<  "The low limits are: "  <<  low_limits;
//    cout <<  "The Upper limits are: "  << upp_limits;
    cout << "The input data file is: "  << input_data_file;

    ndat_inRec.set_size(Number_of_Recessions);
    ndat_inRec.fill(0);
    int info_constr= 0;
    if(reject_outside_bounds_constrait) info_constr =1;

    Model_Comp.set_size(Dim+1,n_one_population,N_Complexes);
    Model_Comp.fill(99999999);

    Model_param_Parents.set_size(Dim +1, all_n_Population);
    Model_param_Parents.fill(99999999);


    best_model.set_size(Dim+1);
    best_model_ens.set_size(Dim+1);

//    test_functions help_fittness(Dim,Nfunc,PRoblem_Number);
 //   help_fittness.initialize();

//    OBj_function = help_fittness;

    number_func_evaluation = 0;
    converged = false;
    path_out_file += "INITIAL";

    switch (SCDE_selection)
    {
        case best1bin:
            SCDE_NAME = "_BEST1BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_ONE-GEN-" << Number_of_generations_in_one_complex << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  << "_Constraint-" << info_constr;
            break;
        case rand1bin:
            SCDE_NAME =  "_RAND1BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_ONE-GEN-" << Number_of_generations_in_one_complex << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR << "_Constraint-" << info_constr;
            break;
        case randbest1bin:
            SCDE_NAME =  "_RANDBEST1BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_ONE-GEN-" << Number_of_generations_in_one_complex << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr;
            break;
        case best2bin:
            SCDE_NAME =  "_BEST2BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_ONE-GEN-" << Number_of_generations_in_one_complex << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr;
            break;
        case rand2bin:
            SCDE_NAME = "_RAND2BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_ONE-GEN-" << Number_of_generations_in_one_complex << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr;
            break;
        case currrand2bin:
            SCDE_NAME = "_CURRAND2BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_ONE-GEN-" << Number_of_generations_in_one_complex << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr;
            break;
        case curbest2bin:
            SCDE_NAME = "_CURBEST2BIN_";
            EN_INFO  << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_ONE-GEN-" << Number_of_generations_in_one_complex << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr;
            break;
         case vg_best1bin:
            SCDE_NAME = "_vg_BEST1BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO <<  "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  << "_Constraint-" << info_constr <<"_vgA-" << vg_A << "_vgR-" << vg_R;
            break;
        case vg_rand1bin:
            SCDE_NAME =  "_vg_RAND1BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR << "_Constraint-" << info_constr<<"_vgA-" << vg_A << "_vgR-" << vg_R;;
            break;
        case vg_randbest1bin:
            SCDE_NAME =  "_vg_RANDBEST1BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr<<"_vgA-" << vg_A << "_vgR-" << vg_R;;
            break;
        case vg_best2bin:
            SCDE_NAME =  "_vg_BEST2BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr<<"_vgA-" << vg_A << "_vgR-" << vg_R;;
            break;
        case vg_rand2bin:
            SCDE_NAME = "_vg_RAND2BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr<<"_vgA-" << vg_A << "_vgR-" << vg_R;;
            break;
        case vg_currrand2bin:
            SCDE_NAME = "_vg_CURRAND2BIN_";
            EN_INFO << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr<<"_vgA-" << vg_A << "_vgR-" << vg_R;;
            break;
        case vg_curbest2bin:
            SCDE_NAME = "_vg_CURBEST2BIN_";
            EN_INFO  << "_DIM-" << Dim << "_ResType-" << ResType <<"_N-COMP-" << N_Complexes << "_SHUF-GEN-" << max_number_of_shuffles;
            EN_INFO << "_POP1COM-" << n_one_population << "_F-" << F << "_CR-" << CR  \
            << "_K-" << K << "_Constraint-" << info_constr<<"_vgA-" << vg_A << "_vgR-" << vg_R;;
            break;
    }


    recesions_data.set_size(1,1);
    sim_obs_recession_data.set_size(1,1);
    sim_obs_recession_data_ens.set_size(1,1);
    PARAM_ens.set_size((Dim+1),ensemble);
    NRECdata = 0;
//    //Sugnathan values of global minimas
//    switch (PRoblem_Number)
//    {
//        case 0:
//            searched_minimum = -450.0;
//            break;
//         case 1:
//            searched_minimum = -450.0;
//            break;
//        case 2:
//            searched_minimum = -450.0;
//            break;
//        case 3:
//            searched_minimum = -450.0;
//            break;
//        case 4:
//            searched_minimum = -310.0;
//            break;
//         case 5:
//            searched_minimum = 390.0;
//            break;
//        case 6:
//            searched_minimum = -180.0;
//            break;
//        case 7:
//            searched_minimum = -140.0;
//            break;
//        case 8:
//            searched_minimum = -330.0;
//            break;
//        case 9:
//            searched_minimum = -330.0;
//            break;
//        case 10:
//            searched_minimum = 90.0;
//            break;
//        default:
//            break;
//    }
//

}

asvh::~asvh()
{
    //dtor
}

asvh::asvh(const asvh& other)
{
    //copy ctor
        N_Complexes = other.N_Complexes;
        Dim = other.Dim;

        lo_limit = other.lo_limit;
        up_limit = other.up_limit;

        CR = other.CR;
        F = other.F;
        K = other.K;

        all_n_Population = other.all_n_Population;
        n_one_population = other.n_one_population;
        all_n_Param = other.all_n_Param;

        nPop = other.nPop;

        Model_param_Parents = other.Model_param_Parents;

        best_model = other.best_model;
        best_model_ens = other.best_model_ens;
        fitness = other.fitness;
        number_func_evaluation = other.number_func_evaluation;
        Number_of_generations_in_one_complex = other.Number_of_generations_in_one_complex;
        max_number_of_shuffles = other.max_number_of_shuffles;

        Model_Comp = other.Model_Comp;

//        PRoblem_Number = other.PRoblem_Number;
//        Nfunc = other.Nfunc;

        ensemble = other.ensemble;
        SCDE_selection = other.SCDE_selection;

//        OBj_function = other.OBj_function;
        converged = other.converged;
        converge_limit = other.converge_limit;

        directory = other.directory;
        path_out_file = other.path_out_file;

        quantile_selection = other.quantile_selection;

        searched_minimum = other.searched_minimum;

        vg_A = other.vg_A;
        vg_R = other.vg_R;

        ResType =other.ResType;
        Number_of_Recessions =other.Number_of_Recessions;
        low_limits = other.low_limits;
        upp_limits = other.upp_limits;
        recesions_data = other.recesions_data;
        input_data_file = other.input_data_file;
        ndat_inRec = other.ndat_inRec;
        NRECdata = other.NRECdata;
        sim_obs_recession_data =other.sim_obs_recession_data;
        PARAM_ens = other.PARAM_ens;
        sim_obs_recession_data_ens = other.sim_obs_recession_data_ens;
}

asvh& asvh::operator=(const asvh& other)
{    //assignment operator
    if (this == &other) {
            return *this; // handle self assignment
    }  else {
        N_Complexes = other.N_Complexes;
        Dim = other.Dim;

        lo_limit = other.lo_limit;
        up_limit = other.up_limit;

        CR = other.CR;
        F = other.F;
        K = other.K;

        all_n_Population = other.all_n_Population;
        n_one_population = other.n_one_population;
        all_n_Param = other.all_n_Param;

        nPop = other.nPop;

        Model_param_Parents = other.Model_param_Parents;

        best_model = other.best_model;
        best_model_ens = other.best_model_ens;
        fitness = other.fitness;
        number_func_evaluation = other.number_func_evaluation;
        Number_of_generations_in_one_complex = other.Number_of_generations_in_one_complex;
        max_number_of_shuffles = other.max_number_of_shuffles;

        Model_Comp = other.Model_Comp;

//        PRoblem_Number = other.PRoblem_Number;
//        Nfunc = other.Nfunc;

        ensemble = other.ensemble;
        SCDE_selection = other.SCDE_selection;

//        OBj_function = other.OBj_function;
        converged = other.converged;
        converge_limit = other.converge_limit;

        directory = other.directory;
        path_out_file = other.path_out_file;

        quantile_selection = other.quantile_selection;
        searched_minimum = other.searched_minimum;

        vg_A = other.vg_A;
        vg_R = other.vg_R;

        ResType =other.ResType;
        Number_of_Recessions =other.Number_of_Recessions;
        low_limits = other.low_limits;
        upp_limits = other.upp_limits;
        recesions_data = other.recesions_data;
        input_data_file = other.input_data_file;
        ndat_inRec = other.ndat_inRec;
        NRECdata = other.NRECdata;
        sim_obs_recession_data =other.sim_obs_recession_data;
        PARAM_ens = other.PARAM_ens;
        sim_obs_recession_data_ens = other.sim_obs_recession_data_ens;
    }
    return *this;
}

/** The initialization of all random Population */
void asvh::initialize_all_POpulation()
{
    //Latin Hypercube Sampling

//Initializing the Model param Parents matrix
       rowvec help_vec;
       rowvec my_init;
       help_vec.set_size(all_n_Population);
       my_init.set_size(all_n_Population);
       double NNd;
       NNd = static_cast<double>(all_n_Population);
       for (unsigned int par =0; par<Dim ;par++ ){
         my_init = random_perm<rowvec>(all_n_Population);
         help_vec.set_size(all_n_Population);
         help_vec.randu();
         my_init = my_init - help_vec;
         my_init = (as_scalar(upp_limits(par)) - as_scalar(low_limits(par))) * my_init / NNd + as_scalar(low_limits(par));
         Model_param_Parents.row(par) = my_init;
         help_vec.reset();
         my_init.reset();
         }

//    cout << endl <<endl<< Model_param_Parents; // tested on LR LR and EXP and NR
    colvec help_col;
    help_col.set_size(Dim);
    for(unsigned int i=0; i< all_n_Population;i++){
        help_col = Model_param_Parents.submat(0,i,Dim-1,i);
        //cout <<endl <<help_col << endl;
        Model_param_Parents(Dim,i) = fittnes(help_col);
    }

    make_Compl_from_mat();
}

/**
  * Random permutation according to Richard Durstenfeld in 1964 in Communications of the ACM volume 7, issue 7, as "Algorithm 235: Random permutation"
  */
  template <class my_vec> my_vec asvh::random_perm(unsigned int n)
{
  unsigned int j=0;
  my_vec shuffled_index, tmp;

  shuffled_index.set_size(n);
  for (unsigned int i=0;i<n ;i++ ){
    shuffled_index(i) = i+1;
    }

   tmp.set_size(1);
   tmp.fill(1);

  for (unsigned int i = n - 1; i > 0; i--) {
    j = rand_int(i);
    tmp(0) = shuffled_index(j);
    shuffled_index(j) = shuffled_index(i);
    shuffled_index(i) = tmp(0);
  }

  return(shuffled_index);
}

/** Random integer number generator */
unsigned int asvh::rand_int(unsigned int n)
{
  unsigned int limit = RAND_MAX - RAND_MAX % n;
  unsigned int rnd;

//  srand((unsigned)time(0));
  //cout << "\n time "<< time(0)<< endl;

  do {
    rnd = rand();
  } while (rnd >= limit);

  return rnd % n;
}

/** After sorting the population according to its fittness -- shuffling the population from Parent matrix to complexes */
void asvh::make_Compl_from_mat()
{
    unsigned int help_var=0;
    for (unsigned int j=0; j< n_one_population ; j++){
       for (unsigned int i =0; i< N_Complexes ;i++ ){
         Model_Comp.subcube(0,j,i,Dim,j,i) = Model_param_Parents.col(help_var);
         help_var++;
         }
      }
}

/** Creating the Parent matrix by combining the Complexes of small Populations */
void asvh::make_mat_from_Compl()
{
    unsigned int left_col=0, right_col;
    right_col = n_one_population -1;
    for (unsigned int j=0; j<N_Complexes ;j++ ){
      Model_param_Parents.submat(0,left_col,Dim,right_col) = Model_Comp.slice(j);
      left_col = right_col+1;
      right_col += n_one_population;
      }

//      cout << endl <<"make_mat_from_compl\n" << Model_param_Parents<< endl ;
}

/** Wrapper for test_fuction class */
double asvh::fittnes(colvec X)
{
    double value =0.0;

    value = compute_fittnes(X);

    number_func_evaluation++;

//    out_fittness(value,X);
    out_fittness(value);

    return(value);
}

/** Sorting columns of Parent matrix according to fittnes */
void asvh::sort_model_param_parents_mat()
{
    rowvec fit_ness;
    fit_ness = Model_param_Parents.row(Dim);

    //cout << fit_ness;
    umat indexes = sort_index(fit_ness);
//    cout <<  "\n sorted indexes\n"  << indexes;

    mat Help_Model_param_parents;
    Help_Model_param_parents.set_size(Dim + 1, all_n_Population);

    for (unsigned int i = 0; i< all_n_Population ;i++ ){
       Help_Model_param_parents.col(i) = Model_param_Parents.col(indexes(i));
      }
    Model_param_Parents = Help_Model_param_parents;
//     cout << Model_param_Parents;
    best_model = Model_param_Parents.col(0);

//    check_convergence();

}

/** Method for checking the converge of computations */
void asvh::check_convergence()
{
    if (as_scalar(best_model(Dim)) <= converge_limit) {
            converged = true;
          //  cout << "\nThe problem has converged.\n";
    }
}

/** Shuffled DE best one */
void asvh::SC_DE_best_one_bin( unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();
    initialize_all_POpulation();

//    cout << endl << path_out_file;
    unsigned int kk=0;
    while( (kk < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_best_one_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            kk++;
    }//end while loop

}

/** Shuffled vgSDE best one vg */
void asvh::vg_SC_DE_best_one_bin( unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/vgSC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();
    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    unsigned int help_vg_A, help_vg_R,coef_multi=0;// vg_A_next;
    help_vg_A = vg_A;
    help_vg_R = vg_R;
//    vg_A_next = vg_A;
//    vg_A_next +=  vg_A_next * pow(vg_R, coef_multi);
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            vg_A = vg_A + vg_R*coef_multi;
//            vg_A_next = vg_A + vg_R*(coef_multi+1);
            Number_of_generations_in_one_complex = vg_A;
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_best_one_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
//            if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
//cout << endl << kk << " " << Number_of_generations_in_one_complex << "  vga " << vg_A <<  " vganext  " << vg_A_next << " " << max_function_eval << " n_ feval " << Number_of_generations_in_one_complex * N_Complexes  * n_one_population <<  " " <<  vg_A_next * N_Complexes  * n_one_population - vg_A * N_Complexes  * n_one_population << endl;
            coef_multi++;
            k++;
    }//end while loop
    vg_A = help_vg_A;
    vg_R = help_vg_R;
}

/** DE/BEST/1/bin */
mat asvh::DE_best_one_bin(mat Population)
{
    mat offsprings;
    mat POpulation;

    POpulation = Population;

    unsigned int number_rows, number_cols;
    unsigned int r1, r2, index_rand;
    double rr = 0.0;

    number_cols = POpulation.n_cols;
    number_rows = POpulation.n_rows;

    offsprings.set_size(number_rows,number_cols);
    offsprings.fill(44);
    //    if(Dim == 2) {
//        cout <<"\nNumber of partial offsprings in one complex is too low: " << n_one_population;
//        exit(EXIT_FAILURE);
//    }
    for (unsigned int j=0;j<Number_of_generations_in_one_complex;j++ ){
            for (unsigned int i=0; i<n_one_population ; i++){
                //generation of neede random indexes
                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r1 = static_cast<int>(floor( (rr * n_one_population)));
                        if (r1 == n_one_population) r1--;
                }while(r1 == i);

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r2 = static_cast<int>(floor( (rr * n_one_population)));
                        if (r2 == n_one_population) r2--;
                }while((r2 == i)||(r2==r1));
            //    cout << r1 <<" i " << i << " population member " << r2 <<endl;
                index_rand = rand_int(Dim);
                for (unsigned int s =0;s <Dim ;s++ ){
                   if((randu_interval()<CR)|| (s==index_rand)){
                     offsprings(s,i) = best_model(s) + F * (POpulation(s,r1) - POpulation(s,r2));
                     if(reject_outside_bounds_constrait){

                        if((as_scalar(offsprings(s,i))< as_scalar(low_limits(s))) || (as_scalar(offsprings(s,i))>as_scalar(upp_limits(s))))
                           offsprings(s,i) = POpulation(s,i);
                     }
                    } else {
                      offsprings(s,i) = POpulation(s,i);
                    }//end if of statement
                  }//end Dim loop
                    colvec help_col;
                    help_col.set_size(Dim);
                    help_col = offsprings.submat(0,i,Dim-1,i);
                    offsprings(Dim,i) = fittnes(help_col);
                  if(as_scalar(offsprings(Dim,i)) < as_scalar(POpulation(Dim,i))){
                    POpulation.col(i) = offsprings.col(i);
                     if(as_scalar(offsprings(Dim,i))<=best_model(Dim)) best_model = offsprings.col(i);
                  }//end of evaluation of fittnes condition
              }//end of offspring population loop
      }//end of max number of function evaluation loop
    return POpulation;
}

/** Random number generator provides double value in interval (0.1) */
double asvh::randu_interval()
{
    double value = 0.0;
    unsigned int LIM_MAX;

    LIM_MAX =   1+ (static_cast<unsigned int>( RAND_MAX ) );
    value = static_cast<double>(rand() ) / (static_cast<double>(LIM_MAX));

    return (value);
}

/** Writing the fitness to the out file */
void asvh::out_fittness(double value)
{
      ofstream out_stream(path_out_file.c_str(), ios::app);
      if (!out_stream) {
       cout << "\nIt is impossible to write to file  " << path_out_file;
       exit(EXIT_FAILURE);
     }
    out_stream << number_func_evaluation << "\t" << (value) <<"\n";
    out_stream.close();
}

/** Writing the fitness to the out file */
void asvh::out_fittness(double value, colvec X)
{
      ofstream out_stream(path_out_file.c_str(), ios::app);
      if (!out_stream) {
       cout << "\nIt is impossible to write to file  " << path_out_file;
       exit(EXIT_FAILURE);
     }
    out_stream << number_func_evaluation << "\t" << (value) << trans(X);
    out_stream.close();
}



/** Shuffled Drandst one bin */
void asvh::SC_DE_rand_one_bin( unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_rand_one_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            k++;
    }//end while loop

}

/** Shuffled vgSDE rand to one bin */
void asvh::vg_SC_DE_rand_one_bin( unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    unsigned int help_vg_A, help_vg_R,coef_multi=0;// vg_A_next;
    help_vg_A = vg_A;
    help_vg_R = vg_R;
//    vg_A_next = vg_A;
//    vg_A_next +=  vg_A_next * pow(vg_R, coef_multi);
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            vg_A = vg_A + vg_R*coef_multi;
//            vg_A_next = vg_A + vg_R*(coef_multi+1);
            Number_of_generations_in_one_complex = vg_A;
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_rand_one_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            coef_multi++;
            k++;
    }//end while loop
    vg_A = help_vg_A;
    vg_R = help_vg_R;
}

/** DE/RAND/1/bin */
mat asvh::DE_rand_one_bin(mat POpulation)
{
    mat offsprings;
    unsigned int number_rows, number_cols;
    unsigned int r1, r2, r3, index_rand;
    double rr = 0.0;

    number_cols = POpulation.n_cols;
    number_rows = POpulation.n_rows;

    offsprings.set_size(number_rows,number_cols);
    offsprings.fill(44);
//
//    if(Dim == 2) {
//        cout <<"\nNumber of partial offsprings in one complex is too low: " << n_one_population;
//        exit(EXIT_FAILURE);
//    }
    for (unsigned int j=0;j<Number_of_generations_in_one_complex;j++ ){
            for (unsigned int i=0; i<n_one_population ; i++){
                //generation of neede random indexes
                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r1 = static_cast<int>(floor( (rr * n_one_population)));
                }while(r1 == i);

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r2 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r2 == i)||(r2==r1));

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r3 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r3 == i)||(r3==r1)||(r3==r2));

            //    cout << r1 <<" i " << i << " population member " << r2 <<endl;
                index_rand = rand_int(Dim);
                for (unsigned int s =0;s <Dim ;s++ ){
                   if((randu_interval()<CR)||  (s==index_rand)){
                     offsprings(s,i) = POpulation(s,r1) + F * (POpulation(s,r2) - POpulation(s,r3));
                     if(reject_outside_bounds_constrait){
                        if((as_scalar(offsprings(s,i))< as_scalar(low_limits(s))) || (as_scalar(offsprings(s,i))>as_scalar(upp_limits(s))))
                           offsprings(s,i) = POpulation(s,i);
                     }
                    } else {
                      offsprings(s,i) = POpulation(s,i);
                    }//end if of statement
                  }//end Dim loop
                    colvec help_col;
                    help_col.set_size(Dim);
                    help_col = offsprings.submat(0,i,Dim-1,i);
                    offsprings(Dim,i) = fittnes(help_col);
                  if(as_scalar(offsprings(Dim,i)) < as_scalar(POpulation(Dim,i))){
                    POpulation.col(i) = offsprings.col(i);
                     if(as_scalar(offsprings(Dim,i))<=best_model(Dim)) best_model = offsprings.col(i);
                  }//end of evaluation of fittnes condition
              }//end of offspring population loop
      }//end of max number of function evaluation loop
    return POpulation;
}

/** DE/rand_to_best/1 */
mat asvh::DE_rand_to_best_one_bin(mat POpulation)
{
    mat offsprings;
    unsigned int number_rows, number_cols;
    unsigned int r1, r2, r3, r4, index_rand;
    double rr = 0.0;

    number_cols = POpulation.n_cols;
    number_rows = POpulation.n_rows;

    offsprings.set_size(number_rows,number_cols);
    offsprings.fill(44);
//
//    if(Dim == 2) {
//        cout <<"\nNumber of partial offsprings in one complex is too low: " << n_one_population;
//        exit(EXIT_FAILURE);
//    }
    for (unsigned int j=0;j<Number_of_generations_in_one_complex;j++ ){
            for (unsigned int i=0; i<n_one_population ; i++){
                //generation of neede random indexes
                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r1 = static_cast<int>(floor( (rr * n_one_population)));
                }while(r1 == i);

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r2 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r2 == i)||(r2==r1));

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r3 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r3 == i)||(r3==r1)||(r3==r2));

                 do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r4 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r4 == i)||(r4==r1)||(r4==r2)||(r4==r3));

            //    cout << r1 <<" i " << i << " population member " << r2 <<endl;
                index_rand = rand_int(Dim);
                for (unsigned int s =0;s <Dim ;s++ ){
                   if((randu_interval()<CR)||  (s==index_rand)){
                     offsprings(s,i) = POpulation(s,r1) + K * (best_model(s) - POpulation(s,r4)) + F * (POpulation(s,r2) - POpulation(s,r3));
                     if(reject_outside_bounds_constrait){
                        if((as_scalar(offsprings(s,i))< as_scalar(low_limits(s))) || (as_scalar(offsprings(s,i))>as_scalar(upp_limits(s))))
                           offsprings(s,i) = POpulation(s,i);
                     }
                    } else {
                      offsprings(s,i) = POpulation(s,i);
                    }//end if of statement
                  }//end Dim loop
                    colvec help_col;
                    help_col.set_size(Dim);
                    help_col = offsprings.submat(0,i,Dim-1,i);
                    offsprings(Dim,i) = fittnes(help_col);
                  if(as_scalar(offsprings(Dim,i)) < as_scalar(POpulation(Dim,i))){
                    POpulation.col(i) = offsprings.col(i);
                     if(as_scalar(offsprings(Dim,i))<=best_model(Dim)) best_model = offsprings.col(i);
                  }//end of evaluation of fittnes condition
              }//end of offspring population loop
      }//end of max number of function evaluation loop
    return POpulation;
}

/** Pure SDE DE/RANDTOBEST/1/BIN */
void asvh::SC_DE_rand_to_best_one_bin( unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();
    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_rand_to_best_one_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();

//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            k++;
    }//end while loop
}

/** Pure vgSDE DE/RANDTOBEST/1/BIN */
void asvh::vg_SC_DE_rand_to_best_one_bin( unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();
    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    unsigned int help_vg_A, help_vg_R,coef_multi=0;// vg_A_next;
    help_vg_A = vg_A;
    help_vg_R = vg_R;
//    vg_A_next = vg_A;
//    vg_A_next +=  vg_A_next * pow(vg_R, coef_multi);
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            vg_A = vg_A + vg_R*coef_multi;
//            vg_A_next = vg_A + vg_R*(coef_multi+1);
            Number_of_generations_in_one_complex = vg_A;
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_rand_to_best_one_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();

//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            coef_multi++;
            k++;
    }//end while loop
    vg_A = help_vg_A;
    vg_R = help_vg_R;
}

/** DE/BEST2RAND/2/BIN */
mat asvh::DE_best_two_bin(mat POpulation)
{
    mat offsprings;
    unsigned int number_rows, number_cols;
    unsigned int r1, r2, r3, r4, index_rand;
    double rr = 0.0;

    number_cols = POpulation.n_cols;
    number_rows = POpulation.n_rows;

    offsprings.set_size(number_rows,number_cols);
    offsprings.fill(44);
//
//    if(Dim == 2) {
//        cout <<"\nNumber of partial offsprings in one complex is too low: " << n_one_population;
//        exit(EXIT_FAILURE);
//    }
    for (unsigned int j=0;j<Number_of_generations_in_one_complex;j++ ){
            for (unsigned int i=0; i<n_one_population ; i++){
                //generation of neede random indexes
                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r1 = static_cast<int>(floor( (rr * n_one_population)));
                }while(r1 == i);

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r2 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r2 == i)||(r2==r1));

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r3 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r3 == i)||(r3==r1)||(r3==r2));

                 do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r4 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r4 == i)||(r4==r1)||(r4==r2)||(r4==r3));

            //    cout << r1 <<" i " << i << " population member " << r2 <<endl;
                index_rand = rand_int(Dim);
                for (unsigned int s =0;s <Dim ;s++ ){
                   if((randu_interval()<CR)||  (s==index_rand)){
                     offsprings(s,i) = best_model(s) + K * ( POpulation(s,r1) - POpulation(s,r4)) + F * (POpulation(s,r2) - POpulation(s,r3));
                     if(reject_outside_bounds_constrait){
                        if((as_scalar(offsprings(s,i))< as_scalar(low_limits(s))) || (as_scalar(offsprings(s,i))>as_scalar(upp_limits(s))))
                           offsprings(s,i) = POpulation(s,i);
                     }
                    } else {
                      offsprings(s,i) = POpulation(s,i);
                    }//end if of statement
                  }//end Dim loop
                    colvec help_col;
                    help_col.set_size(Dim);
                    help_col = offsprings.submat(0,i,Dim-1,i);
                    offsprings(Dim,i) = fittnes(help_col);
                  if(as_scalar(offsprings(Dim,i)) < as_scalar(POpulation(Dim,i))){
                    POpulation.col(i) = offsprings.col(i);
                     if(as_scalar(offsprings(Dim,i))<=best_model(Dim)) best_model = offsprings.col(i);
                  }//end of evaluation of fittnes condition
              }//end of offspring population loop
      }//end of max number of function evaluation loop
    return POpulation;
}

/** Pure SDE DE/BEST2RAND/2/BIN */
void asvh::SC_DE_best_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_best_two_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            k++;
    }//end while loop
}

/** Pure vgSDE DE/BEST2RAND/2/BIN */
void asvh::vg_SC_DE_best_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    unsigned int help_vg_A, help_vg_R,coef_multi=0;// vg_A_next;
    help_vg_A = vg_A;
    help_vg_R = vg_R;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            vg_A = vg_A + vg_R*coef_multi;
//            vg_A_next = vg_A + vg_R*(coef_multi+1);
            Number_of_generations_in_one_complex = vg_A;
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_best_two_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            coef_multi++;
            k++;
    }//end while loop
    vg_A = help_vg_A;
    vg_R = help_vg_R;
}

/** DE/RAND2/2/BIN */
mat asvh::DE_rand_two_bin(mat POpulation)
{
    mat offsprings;
    unsigned int number_rows, number_cols;
    unsigned int r1, r2, r3, r4, r5, index_rand;
    double rr = 0.0;

    number_cols = POpulation.n_cols;
    number_rows = POpulation.n_rows;

    offsprings.set_size(number_rows,number_cols);
    offsprings.fill(44);
//    if(Dim == 2) {
//        cout <<"\nNumber of partial offsprings in one complex is too low: " << n_one_population;
//        exit(EXIT_FAILURE);
//    }
    for (unsigned int j=0;j<Number_of_generations_in_one_complex;j++ ){
            for (unsigned int i=0; i<n_one_population ; i++){
                //generation of neede random indexes
                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r1 = static_cast<int>(floor( (rr * n_one_population)));
                }while(r1 == i);

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r2 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r2 == i)||(r2==r1));

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r3 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r3 == i)||(r3==r1)||(r3==r2));

                 do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r4 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r4 == i)||(r4==r1)||(r4==r2)||(r4==r3));

                 do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r5 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r5 == i)||(r5==r1)||(r5==r2)||(r5==r3)||(r5==r4));
            //    cout << r1 <<" i " << i << " population member " << r2 <<endl;
                index_rand = rand_int(Dim);
                for (unsigned int s =0;s <Dim ;s++ ){
                   if((randu_interval()<CR)||  (s==index_rand)){
                     offsprings(s,i) = POpulation(s,r1) + K * ( POpulation(s,r5) - POpulation(s,r4)) + F * (POpulation(s,r2) - POpulation(s,r3));
                     if(reject_outside_bounds_constrait){
                        if((as_scalar(offsprings(s,i))< as_scalar(low_limits(s))) || (as_scalar(offsprings(s,i))>as_scalar(upp_limits(s))))
                           offsprings(s,i) = POpulation(s,i);
                     }
                    } else {
                      offsprings(s,i) = POpulation(s,i);
                    }//end if of statement
                  }//end Dim loop
                    colvec help_col;
                    help_col.set_size(Dim);
                    help_col = offsprings.submat(0,i,Dim-1,i);
                    offsprings(Dim,i) = fittnes(help_col);
                  if(as_scalar(offsprings(Dim,i)) < as_scalar(POpulation(Dim,i))){
                    POpulation.col(i) = offsprings.col(i);
                     if(as_scalar(offsprings(Dim,i))<=best_model(Dim)) best_model = offsprings.col(i);
                  }//end of evaluation of fittnes condition
              }//end of offspring population loop
      }//end of max number of function evaluation loop
    return POpulation;
}

/** Pure SDE DE/RAND2/2/BIN */
void asvh::SC_DE_rand_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_rand_two_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            k++;
    }//end while loop
}

/** Pure vgSDE DE/RAND2/2/BIN */
void asvh::vg_SC_DE_rand_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    unsigned int help_vg_A, help_vg_R,coef_multi=0;// vg_A_next;
    help_vg_A = vg_A;
    help_vg_R = vg_R;
//    vg_A_next = vg_A;
//    vg_A_next +=  vg_A_next * pow(vg_R, coef_multi);
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            vg_A = vg_A + vg_R*coef_multi;
//            vg_A_next = vg_A + vg_R*(coef_multi+1);
            Number_of_generations_in_one_complex = vg_A;
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_rand_two_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
    //        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            coef_multi++;
            k++;
    }//end while loop
    vg_A = help_vg_A;
    vg_R = help_vg_R;
}

/** DE/CURRENTTORAND/2/BIN */
mat asvh::DE_current_to_rand_2_bin(mat POpulation)
{
    mat offsprings;
    unsigned int number_rows, number_cols;
    unsigned int r1, r2, r3, index_rand;
    double rr = 0.0;

    number_cols = POpulation.n_cols;
    number_rows = POpulation.n_rows;

    offsprings.set_size(number_rows,number_cols);
    offsprings.fill(44);
//
//    if(Dim == 2) {
//        cout <<"\nNumber of partial offsprings in one complex is too low: " << n_one_population;
//        exit(EXIT_FAILURE);
//    }
    for (unsigned int j=0;j<Number_of_generations_in_one_complex;j++ ){
            for (unsigned int i=0; i<n_one_population ; i++){
                //generation of neede random indexes
                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r1 = static_cast<int>(floor( (rr * n_one_population)));
                }while(r1 == i);

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r2 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r2 == i)||(r2==r1));

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r3 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r3 == i)||(r3==r1)||(r3==r2));

            //    cout << r1 <<" i " << i << " population member " << r2 <<endl;
                index_rand = rand_int(Dim);
                for (unsigned int s =0;s <Dim ;s++ ){
                   if((randu_interval()<CR)|| (s==index_rand)){
                     offsprings(s,i) = POpulation(s,i) + K*(POpulation(s,i) - POpulation(s,r1)) + F * (POpulation(s,r2) - POpulation(s,r3));
                     if(reject_outside_bounds_constrait){
                        if((as_scalar(offsprings(s,i))< as_scalar(low_limits(s))) || (as_scalar(offsprings(s,i))>as_scalar(upp_limits(s))))
                           offsprings(s,i) = POpulation(s,i);
                     }
                    } else {
                      offsprings(s,i) = POpulation(s,i);
                    }//end if of statement
                  }//end Dim loop
                    colvec help_col;
                    help_col.set_size(Dim);
                    help_col = offsprings.submat(0,i,Dim-1,i);
                    offsprings(Dim,i) = fittnes(help_col);
                  if(as_scalar(offsprings(Dim,i)) < as_scalar(POpulation(Dim,i))){
                    POpulation.col(i) = offsprings.col(i);
                     if(as_scalar(offsprings(Dim,i))<=best_model(Dim)) best_model = offsprings.col(i);
                  }//end of evaluation of fittnes condition
              }//end of offspring population loop
      }//end of max number of function evaluation loop
    return POpulation;
}

/** Pure SDE DE/CURRENTTORAND/2/BIN */
void asvh::SC_DE_current_to_rand_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_current_to_rand_2_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
//        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            k++;
    }//end while loop
}

/** Pure vgSDE DE/CURRENTTORAND/2/BIN */
void asvh::vg_SC_DE_current_to_rand_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    unsigned int help_vg_A, help_vg_R,coef_multi=0;// vg_A_next;
    help_vg_A = vg_A;
    help_vg_R = vg_R;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            vg_A = vg_A + vg_R*coef_multi;
//            vg_A_next = vg_A + vg_R*(coef_multi+1);
            Number_of_generations_in_one_complex = vg_A;
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_current_to_rand_2_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
//        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            coef_multi++;
            k++;
    }//end while loop
    vg_A = help_vg_A;
    vg_R = help_vg_R;
}


/** DE/CURRENTTOBEST/2/BIN */
mat asvh::DE_current_to_best_2_bin(mat POpulation)
{
     mat offsprings;
    unsigned int number_rows, number_cols;
    unsigned int r1, r2, index_rand;
    double rr = 0.0;

    number_cols = POpulation.n_cols;
    number_rows = POpulation.n_rows;

    offsprings.set_size(number_rows,number_cols);
    offsprings.fill(44);
//    if(Dim == 2) {
//        cout <<"\nNumber of partial offsprings in one complex is too low: " << n_one_population;
//        exit(EXIT_FAILURE);
//    }
    for (unsigned int j=0;j<Number_of_generations_in_one_complex;j++ ){
            for (unsigned int i=0; i<n_one_population ; i++){
                //generation of neede random indexes
                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r1 = static_cast<int>(floor( (rr * n_one_population)));
                }while(r1 == i);

                do{
                        rr = static_cast<double>(rand() ) / (static_cast<double>(RAND_MAX));
                        r2 = static_cast<int>(floor( (rr * n_one_population)));
                }while((r2 == i)||(r2==r1));
            //    cout << r1 <<" i " << i << " population member " << r2 <<endl;
                index_rand = rand_int(Dim);
                for (unsigned int s =0;s <Dim ;s++ ){
                   if((randu_interval()<CR)|| (s==index_rand)){
                     offsprings(s,i) = POpulation(s,i) + K * ( best_model(s) - POpulation(s,i)) + F * (POpulation(s,r1) - POpulation(s,r2));
                     if(reject_outside_bounds_constrait){
                        if((as_scalar(offsprings(s,i))< as_scalar(low_limits(s))) || (as_scalar(offsprings(s,i))>as_scalar(upp_limits(s))))
                           offsprings(s,i) = POpulation(s,i);
                     }
                    } else {
                      offsprings(s,i) = POpulation(s,i);
                    }//end if of statement
                  }//end Dim loop
                    colvec help_col;
                    help_col.set_size(Dim);
                    help_col = offsprings.submat(0,i,Dim-1,i);
                    offsprings(Dim,i) = fittnes(help_col);
                  if(as_scalar(offsprings(Dim,i)) < as_scalar(POpulation(Dim,i))){
                    POpulation.col(i) = offsprings.col(i);
                     if(as_scalar(offsprings(Dim,i))<=best_model(Dim)) best_model = offsprings.col(i);
                  }//end of evaluation of fittnes condition
              }//end of offspring population loop
      }//end of max number of function evaluation loop
    return POpulation;
}

/** Pure SDE DE/CURRENTTOBEST/2/BIN */
void asvh::SC_DE_current_to_best_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_current_to_best_2_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
//        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            k++;
    }//end while loop
}

/** Pure vgSDE DE/CURRENTTOBEST/2/BIN */
void asvh::vg_SC_DE_current_to_best_two_bin(unsigned int ensemble_number)
{
    stringstream EN_INFO1;
    EN_INFO1 << "_ens_" << ensemble_number;

    path_out_file +=  "/SC";
    path_out_file += SCDE_NAME.c_str() +  EN_INFO.str() + EN_INFO1.str();

    initialize_all_POpulation();
//    cout << endl << path_out_file;
    unsigned int k=0;
    unsigned int help_vg_A, help_vg_R,coef_multi=0;// vg_A_next;
    help_vg_A = vg_A;
    help_vg_R = vg_R;
    while( (k < max_number_of_shuffles)){
            sort_model_param_parents_mat();
            make_Compl_from_mat();
            vg_A = vg_A + vg_R*coef_multi;
//            vg_A_next = vg_A + vg_R*(coef_multi+1);
            Number_of_generations_in_one_complex = vg_A;
            for (unsigned int i =0;i<N_Complexes ; i++){
              Model_Comp.slice(i) = DE_current_to_best_2_bin(Model_Comp.slice(i));
            }
            make_mat_from_Compl();
//            cout << boolalpha << "\nSC DE/best/1 -- " << k <<" iteration of the main loop. \nConvergence Status --> "<< converged;
//        if(k==4) converged = true;
//            if(converged) {
//                cout <<"\n\nThe computation reached the convergence with number of function eval: " << number_func_evaluation << "."<< endl;
//                cout << Model_Comp << endl << Model_param_Parents;
//                break;
//            }//end for loop
            coef_multi++;
            k++;
    }//end while loop
    vg_A = help_vg_A;
    vg_R = help_vg_R;
}

/** The general interface to the problem testings
 SCDE_selection_0-best1bin_1-rand1bin_2-randbest1bin_3-best2bin_4-rand2bin_5-currrand2bin_6-curbest2bin*/
void asvh::run_ensemble()
{
    cout << "\n\n**** ENSEMBLE RUN ****\n";
    string filet_names , filet_out_names;//file with names of ensemble outcomes

    filet_names +=directory + "/convergence/File_names" + SCDE_NAME.c_str() + EN_INFO.str();
    filet_out_names +=directory + "/outflows/AAL_ENS_BEST_OUT_DAT" + SCDE_NAME.c_str() + EN_INFO.str();

    ofstream filet_name_stream(filet_names.c_str(), ios::app);
      if (!filet_name_stream) {
       cout << "\nIt is impossible to write to file with nams " << filet_names.c_str()<< endl;
       exit(EXIT_FAILURE);
     }

     best_model_ens.fill(9999999);
     string filet_out_names_BE, filet_out_names_PAR;
     filet_out_names_BE +=directory + "/outflows/BEST_models_ENSEMBLE";
     filet_out_names_PAR +=directory + "/params/BEST_models_ENSEMBLE_PARAMS";

     mat MAT_SIM_REC;

    for(unsigned int i=0; i<ensemble;i++){
            srand(i+1);
            if (!path_out_file.empty()) path_out_file.clear();
            path_out_file += directory + "/convergence";
            cout << "\nEnsemble simulation: " << i+1;
            number_func_evaluation = 0;
            switch (SCDE_selection)
            {
                case best1bin:
                    SC_DE_best_one_bin(i);
                    break;
                case rand1bin:
                    SC_DE_rand_one_bin(i);
                    break;
                case randbest1bin:
                    SC_DE_rand_to_best_one_bin(i);
                    break;
                case best2bin:
                    SC_DE_best_two_bin(i);
                    break;
                case rand2bin:
                    SC_DE_rand_two_bin(i);
                    break;
                case currrand2bin:
                    SC_DE_current_to_rand_two_bin(i);
                    break;
                case curbest2bin:
                    SC_DE_current_to_best_two_bin(i);
                    break;
                case vg_best1bin:
                    vg_SC_DE_best_one_bin(i);
                    break;
                case vg_rand1bin:
                    vg_SC_DE_rand_one_bin(i);
                    break;
                case vg_randbest1bin:
                    vg_SC_DE_rand_to_best_one_bin(i);
                    break;
                case vg_best2bin:
                    vg_SC_DE_best_two_bin(i);
                    break;
                case vg_rand2bin:
                    vg_SC_DE_rand_two_bin(i);
                    break;
                case vg_currrand2bin:
                    vg_SC_DE_current_to_rand_two_bin(i);
                    break;
                case vg_curbest2bin:
                    vg_SC_DE_current_to_best_two_bin(i);
                    break;
                 default:
                    break;
            }
            filet_name_stream << path_out_file.c_str() << endl;
            cout << endl <<"Best model: " << (as_scalar(best_model(Dim))) << endl;
//            out_sim_data(best_model, filet_out_names_BE);
            MAT_SIM_REC = out_sim_data_MAT(best_model);
            sim_obs_recession_data_ens.col(i) = MAT_SIM_REC.col(3);
            PARAM_ens.col(i) = best_model;
            if(as_scalar(best_model(Dim))  < as_scalar(best_model_ens(Dim)) ) best_model_ens = best_model;
            best_model.fill(9999999);
        }
        filet_name_stream.close();

      // cout << endl << "\n***The ensemble run has been finished.***" << endl;

     //  calc_quantiles_fast_small();

      // cout << endl << "***The quantiles has been computed.***\n" << endl;

       cout << "\n*************The BEST MODEL ENSEMBLE, reservoir parameters********\n";
    //   cout << best_model_ens;

       out_sim_data(best_model_ens, filet_out_names);


     ofstream out_DATA(filet_out_names_BE.c_str(), ios::app);
      if (!out_DATA) {
       cout << "\nIt is impossible to write to file  " << filet_out_names_BE;
       exit(EXIT_FAILURE);
     }
     sim_obs_recession_data_ens.reshape(NRECdata,(ensemble+3));
     sim_obs_recession_data_ens.col(ensemble) = sim_obs_recession_data.col(2);
     sim_obs_recession_data_ens.col(ensemble+1) = sim_obs_recession_data.col(1);
     sim_obs_recession_data_ens.col(ensemble+2) = sim_obs_recession_data.col(0);
     out_DATA << sim_obs_recession_data_ens;
     out_DATA.close();

    ofstream out_PAR(filet_out_names_PAR.c_str(), ios::app);
      if (!out_PAR) {
       cout << "\nIt is impossible to write to file  " << filet_out_names_PAR;
       exit(EXIT_FAILURE);
     }
     out_PAR << PARAM_ens;
     out_PAR.close();

}

/** Removes the duplicated elements in my_vec
 the implementation for type 7 of estimation of quantiles
TODO implementa rest quantile estimators
*/
template<class my_vec>  my_vec asvh::unique_(my_vec data)
{
   data = sort(data);

   my_vec return_vec;
   return_vec.set_size(data.n_elem);
   return_vec.fill(999999);

   unsigned int number_d =0, help_ind =0, help_var =0;
   number_d = data.n_elem;

    for (unsigned int i=0;i< number_d;i++ ){
           while (as_scalar(data(help_ind)) == as_scalar(data(help_var)) ){
               help_var++;
           }
           help_ind += (help_var - help_ind) + 1;
           help_var = help_ind;
           return_vec(help_ind) = data(help_ind);
      }

   return_vec.reshape(help_ind);

   return(return_vec);
}

/** Returns vector of quantiles for given data
IMPLEMENTED
  Discontinuous:
       INV_emp_dist - inverse of the empirical distribution function

 NOT IMPLEMENTED
  //   INV_emp_dist_av - like type 1, but with averaging at discontinuities (g=0)
  //   SAS_near_int - SAS definition: nearest even order statistic
   //  Piecwise linear continuous:
   //    In this case, sample quantiles can be obtained by linear interpolation
   //    between the k-th order statistic and p(k).
   //    type=4 - linear interpolation of empirical cdf, p(k)=k/n;
   //    type=5 - hydrolgical a very popular definition, p(k) = (k-0.5)/n;
   //    type=6 - used by Minitab and SPSS, p(k) = k/(n+1);
   //    type=7 - used by S-Plus and R, p(k) = (k-1)/(n-1);
   //    type=8 - resulting sample quantiles are approximately median unbiased
   //             regardless of the distribution of x. p(k) = (k-1/3)/(n+1/3);
   //    type=9 - resulting sample quantiles are approximately unbiased, when
   //             the sample comes from Normal distribution. p(k)=(k-3/8)/(n+1/4);
   //
   //    default type = 7
   //
    References:
    1) Hyndman, R.J and Fan, Y, (1996) "Sample quantiles in statistical packages"
                                        American Statistician, 50, 361-365

*/
template<class my_vec> my_vec  asvh::quantiles_(my_vec data, my_vec probs)
{
    my_vec sorted_data;
    my_vec quantile;

    sorted_data = sort(data);
    quantile.set_size(probs.n_elem);

    uvec index;
    index.set_size(probs.n_elem);

    switch (quantile_selection)
    {
        case INV_emp_dist:
            for (unsigned int i=0; i<probs.n_elem  ;i++ ){
                    if(as_scalar(probs(i)) == 0.0) index(i) = 0;
                        else index(i) = static_cast<unsigned int> (ceil((as_scalar(probs(i) * sorted_data.n_elem)))) -1;
                   }
            for (unsigned int i=0; i<probs.n_elem  ;i++ ){
                      quantile(i) = sorted_data(index(i));
                  }
            break;

//        case INV_emp_dist_av:
//
//            break;
//        default:
//            break;
    }

 //   cout << index;
    return(quantile);
}

/** Slow Quantile estimation for large problems*/
void asvh::calc_quantiles()
{
/*     working example */
//     colvec data;
//     data << 1<< 2<< 3<< 4<< 5 << 6 << 7 << 8<< 9 << 10 <<11 << 12 << 13 << 14 << 15 << 16 << 17 << 18 << 19 << 20 << 21 << 22 << 23 << 24 << 25;
     rowvec prst;
     prst << 0.<< 0.25<< 0.5 << 0.75 << 1.0;

//     cout << data;
//     cout << endl << quantiles_(data, prst);
   string filet_names;//file with names of ensemble outcomes

    filet_names +=directory + "/File_names" + SCDE_NAME.c_str() + EN_INFO.str();
    //filet_names += directory + "/File_names" + "_SC" + SCDE_NAME.c_str() + EN_INFO.str();

 //   cout << filet_names;

    string *file_en, trash;
    unsigned int number_of_files = 0;

    ifstream de_file_ens_names(filet_names.c_str());
    while(de_file_ens_names.good()){
    de_file_ens_names >> trash;
    number_of_files++;

    }
    number_of_files--;
    cout << endl << "Opening the "<< number_of_files << endl;
    file_en = new string[number_of_files];

    de_file_ens_names.clear();
    de_file_ens_names.seekg (0, ios::beg);

    for (unsigned int i =0; i< number_of_files ;i++ ){
      de_file_ens_names >> file_en[i] ;
    //  cout << file_en[i] << endl;
      }
    de_file_ens_names.close();

    ifstream first_filet(file_en[0].c_str());
    unsigned int num_rows =0;
    while(first_filet.good()){
        getline(first_filet,trash);
        num_rows++;
        }
    first_filet.close();
    num_rows --;

    rowvec Quantiles_ens;
    Quantiles_ens.set_size(prst.n_elem+2);

    rowvec help_vec1;
    mat  data_one_file;

    string filet_quantiles;
    filet_quantiles +=directory + "/Ensemble_quantiles_SC_" +SCDE_NAME.c_str() + EN_INFO.str();

    ofstream out_quantile_stream(filet_quantiles.c_str(), ios::app);
      if (!out_quantile_stream) {
       cout << "\nIt is impossible to write to file  " << filet_quantiles;
       exit(EXIT_FAILURE);
      }

    cout << "\nComputing the Quantiles\n";
    for (unsigned int j=0; j<num_rows ; j++){
      help_vec1.set_size(number_of_files);
      for (unsigned int i =0; i<number_of_files; i++){
             data_one_file.load(file_en[i].c_str());
             help_vec1(i) = data_one_file(j,1);
             if(j==0) cout <<  i;
             data_one_file.reset();
             }
             cout << endl<< help_vec1;
      Quantiles_ens.subvec(0,prst.n_elem-1) = quantiles_<rowvec>(help_vec1, prst);
      Quantiles_ens(prst.n_elem) = mean(help_vec1);
      Quantiles_ens(prst.n_elem+1) = stddev(help_vec1);
      for(unsigned int qu =0; qu < (prst.n_elem+2); qu++){
        out_quantile_stream << Quantiles_ens(qu) << "\t";
      }
      out_quantile_stream << "\n";
      help_vec1.reset();
     }

    delete [] file_en;
    out_quantile_stream.close();
}

/** Fast Quantile estimation for middle problems*/
void asvh::calc_quantiles_fast_small()
{
/*     working example */
//     colvec data;
//     data << 1<< 2<< 3<< 4<< 5 << 6 << 7 << 8<< 9 << 10 <<11 << 12 << 13 << 14 << 15 << 16 << 17 << 18 << 19 << 20 << 21 << 22 << 23 << 24 << 25;
     rowvec prst;
     prst << 0.0 << 0.25<< 0.5 << 0.75 << 1.0;

//     cout << data;
//     cout << endl << quantiles_(data, prst);
   string filet_names;//file with names of ensemble outcomes

    filet_names +=directory + "/File_names" + SCDE_NAME.c_str() + EN_INFO.str();
    //filet_names += directory + "/File_names" + "_SC" + SCDE_NAME.c_str() + EN_INFO.str();

 //   cout << filet_names;

    string *file_en, trash;
    unsigned int number_of_files = 0;

    ifstream de_file_ens_names(filet_names.c_str());
    while(de_file_ens_names.good()){
       de_file_ens_names >> trash;
       number_of_files++;
      }
    number_of_files--;
//    cout << endl << "Opening the "<< number_of_files << endl;
    file_en = new string[number_of_files];

    de_file_ens_names.clear();
    de_file_ens_names.seekg (0, ios::beg);

    for (unsigned int i =0; i< number_of_files ;i++ ){
      de_file_ens_names >> file_en[i] ;
    //  cout << file_en[i] << endl;
      }
    de_file_ens_names.close();

    ifstream first_filet(file_en[0].c_str());
    unsigned int num_rows =0;
    while(first_filet.good()){
        getline(first_filet,trash);
        num_rows++;
        }
    first_filet.close();
    num_rows --;
   // cout << num_rows << endl;

    mat Quantiles_ens;
    Quantiles_ens.set_size(num_rows, prst.n_elem+2);



    mat  data_all_file;
    data_all_file.set_size(num_rows,number_of_files);

      for (unsigned int i =0; i<number_of_files; i++){
             mat data_one_file;
             data_one_file.load(file_en[i].c_str());
             data_all_file.submat(0,i,num_rows-1,i) = data_one_file.col(1);
             }

      //       cout <<"\nAll mat finallize.d\n";
    string filet_quantiles;
    filet_quantiles +=directory + "/Ensemble_quantiles_SC_" +SCDE_NAME.c_str() + EN_INFO.str();

//    ofstream out_quantile_stream(filet_quantiles.c_str(), ios::app);
//      if (!out_quantile_stream) {
//       cout << "\nIt is impossible to write to file  " << filet_quantiles;
//       exit(EXIT_FAILURE);
//      }

    cout << "\nComputing the Quantiles\n";
    for (unsigned int j=0; j<num_rows ; j++){
      rowvec help_vec1;
      help_vec1.set_size(number_of_files);
      help_vec1 = data_all_file.row(j);
      Quantiles_ens.submat(j,0,j,prst.n_elem-1) = quantiles_<rowvec>(help_vec1, prst);
      Quantiles_ens(j,prst.n_elem) = mean(help_vec1);
      Quantiles_ens(j,prst.n_elem+1) = stddev(help_vec1);
     }

    Quantiles_ens.save(filet_quantiles.c_str(),raw_ascii);
    delete [] file_en;
//    out_quantile_stream.close();
}

/**
 *  The stop on border constraint handling
 *@param colvec of model params
 */
colvec asvh::contraint_check_stop_on_boundary(colvec MODEL)
{
    for (unsigned int par=0; par<Dim ; par++ ){
       if(as_scalar(MODEL(par)) < as_scalar(low_limits(par)))  MODEL(par) = low_limits(par);
       if(as_scalar(MODEL(par)) > as_scalar(upp_limits(par))) MODEL(par) = low_limits(par);
      }

      return MODEL;
}


/**
 * Estimates the outflow from single resevoir
 * @param reservoir params
  *@param storage of reservoir
 */
double asvh::res_outflow(colvec Par, double S)
{
    double out=0.0;
    switch (ResType){
        case LinRes:
            out = (as_scalar(Par(0)) )* S;
          break;
        case TresLinRes:
            out = (1/(as_scalar(Par(0)) )) * S;
          break;
        case NonLinRes:
            out = (1/(as_scalar(Par(0)) )) * pow (S, (as_scalar(Par(1)) ));
         break;
       case ExpRes:
           out = (1/ (as_scalar(Par(0)) )) * exp (S / (as_scalar(Par(1)) ));
         break;
 //       default:
 //        break;
    }

   return(out);
}

/**
 * Give the reservoirs params
 *@param clovec of model params
 */
colvec asvh::params_RESERVOIR(colvec Par)
{
    colvec params;

        switch (ResType){
        case LinRes:
            params.set_size(1);
            params(0) = Par(Dim-1) ;
          break;
        case TresLinRes:
            params.set_size(1);
            params(0) = Par(Dim-1) ;
          break;
        case NonLinRes:
            params.set_size(2);
            params(0) = Par(Dim-2) ;
            params(1) = Par(Dim-1);
         break;
       case ExpRes:
            params.set_size(2);
            params(0) = Par(Dim-2) ;
            params(1) = Par(Dim-1);
         break;
 //       default:
 //        break;
    }
    return( params);
}

/**
 * Computes the Sum of Squares for 1 recession
  *@param matrxi with data from one recession
  *@param vector of reservoir parameters selected form param colvec
  *@param intial storage for given recession
  */
double asvh::compute_ss_outflow_1Reccesion(mat Qdata, colvec Par, double So)
{
    double S;
    colvec pars;
    unsigned int nr;

    S = So;
    nr = Qdata.n_rows;

    colvec QQ_rec;//recesion outflow
    QQ_rec.set_size(nr);
    QQ_rec.fill(-99999999.9);
    double helpQ=0.0;

    for (unsigned int qdat =0 ;qdat <nr ; qdat++){
      if( S>0) {
            helpQ = res_outflow(Par, S);
            if(helpQ < S) {
                QQ_rec(qdat) = helpQ;
            } else {
                QQ_rec(qdat) = S;
            }
      } else QQ_rec(qdat) = 0.0;
      if((S - as_scalar(QQ_rec(qdat)))>0){
       S = S - as_scalar(QQ_rec(qdat));
      } else {
       S = 0.0;
      }
     }
     vec Resids;

     Resids = log(Qdata.col(0)) - log(QQ_rec);

     double SS=0.0;
     for(unsigned int dd=0;dd<nr;dd++){
        SS= SS + as_scalar(Resids(dd)) * as_scalar(Resids(dd)) ;
     }
//     cout <<"\nidhchjhsd\n";
     return(SS);
}

/**
 *Computes the Sum of Squares for given model
 *@param vector of parameters selected from the population
 */
double asvh::compute_fittnes(colvec Par)
{
     double SumSq = 0.0;
     mat One_QRecDATA;
     double So=-99999.9;
     colvec Par_to_rec;

     for(unsigned rc=0;rc<Number_of_Recessions;rc++){
//            cout <<"\nidhchjhsd11111\n";
        One_QRecDATA = get_one_recession(rc);
//     cout <<"\nidhchjhsd1111122222\n";
        So = as_scalar(Par(rc));
        Par_to_rec = params_RESERVOIR(Par);
        SumSq += compute_ss_outflow_1Reccesion(One_QRecDATA, Par_to_rec, So);
        One_QRecDATA.reset();
     }

     return(SumSq);

}

/**
 * Loading the recesion data and counting length of each of reccion
 *
  */
void asvh::get_recessions_data()
{

    recesions_data.load(input_data_file.c_str());
    unsigned int nr,nc, index_last;

    nr = recesions_data.n_rows;
    nc = recesions_data.n_cols;

//    cout << endl << nc << "    " << nr<< endl;
//    cout << recesions_data;

    if(nc!=3) {
      cout << "\nThe wrong number of columns in :"<< input_data_file.c_str() <<" it should be 3.";
      exit(EXIT_FAILURE);
    }

    index_last = static_cast<unsigned int>(as_scalar(recesions_data((nr-1),1)));
 //   cout <<endl << "upsssss\n" << nc << nr << endl<<index_last << endl;
    //cout << endl <<recesions_data ;


    unsigned int help_indd=0;

    for (unsigned int dat=1; dat <nr ;dat++ ){
            if(((as_scalar(recesions_data((dat),1))-as_scalar(recesions_data((dat-1),1)) )>0)&&(dat<(nr-1))){
                ndat_inRec(help_indd) = recesions_data((dat-1),2);
                help_indd++;
//                cout << help_indd << endl;
            }

      if(dat==(nr-1) ){
            ndat_inRec(help_indd) = as_scalar(recesions_data((nr-1),2));
//            cout << help_indd;
      }
      }
    NRECdata = nr;
    sim_obs_recession_data.set_size(NRECdata ,nc+1);
    sim_obs_recession_data.fill(-999999999.9);
//  cout << endl << "uiohyasduifaghf]\n" << ndat_inRec;
    sim_obs_recession_data_ens.set_size(NRECdata , (ensemble));
    sim_obs_recession_data_ens.fill(-99999999.9);
}

/**
 * filters form recession matrrix one recession
 *@param numeber of selected recesion
 */
mat asvh::get_one_recession(unsigned int Ind_rec)
{
    mat data_1recession;
    data_1recession.set_size(ndat_inRec(Ind_rec),3);
    unsigned int help_ind=0;
//    cout << ndat_inRec(Ind_rec) << endl;
    for (unsigned int dd=0;dd<NRECdata ;dd++ ){
        if((static_cast<unsigned int>(as_scalar(recesions_data(dd,1))))==(Ind_rec+1)) {
                data_1recession.row(help_ind) = recesions_data.row(dd);
//            cout << help_ind << "   " << data_1recession.row(help_ind);
//cout << help_ind << endl;
                help_ind++;
        }
    }
//    cout <<"\nidhchjhsd_GETREC\n";

//    cout << endl << ndat_inRec << endl;
    return data_1recession;
}

/**
 * Compute the the Outflow from the one recession
 * @param matrix with ifnormation about 1 recession
 * @param parameters of constituitive eq. obtained from param colvec
  *@param initial storage obtained for parameters colvec
 */
mat asvh::compute_outflow_1Reccesion(mat Qdata, colvec Par, double So)
{

    colvec pars;
    unsigned int nr;
    double S;

    S = So;
    nr = Qdata.n_rows;

    colvec QQ_rec;//recesion outflow
    QQ_rec.set_size(nr);
    QQ_rec.fill(-99999999.9);
    double helpQ=0.0;

    for (unsigned int qdat =0 ;qdat <nr ; qdat++){
      if( S>0) {
            helpQ = res_outflow(Par, S);
            if(helpQ < S) {
                QQ_rec(qdat) = helpQ;
            } else {
                QQ_rec(qdat) = S;
            }
      } else QQ_rec(qdat) = 0.0;
      if((S - as_scalar(QQ_rec(qdat)))>0){
       S = S - as_scalar(QQ_rec(qdat));
      } else {
       S = 0.0;
      }
     }

     mat Qobs_Qsim_data;
     Qobs_Qsim_data.set_size(nr, 4);
     Qobs_Qsim_data.col(0) = Qdata.col(1);
     Qobs_Qsim_data.col(1) = Qdata.col(2);
     Qobs_Qsim_data.col(2) = Qdata.col(0);
     Qobs_Qsim_data.col(3) = QQ_rec;

//     cout <<"\nidhchjhsd\n";
     return(Qobs_Qsim_data);
}

/**
 * Calculates the outcome for model provided the Par vector
 *@param columumn vector of parameters
 *@param file for writing the simulated outflow
 */
void asvh::out_sim_data(colvec Par, string FILET)
{
    mat One_QRecDATA_IN;
    mat One_QRecDATA_OUT;
    colvec Par_to_rec;
    double So=-999999999.9;
    unsigned int u_row, l_row;

    for(unsigned rc=0;rc<Number_of_Recessions;rc++){
//            cout <<"\nidhchjhsd11111\n";
        if(rc ==0){
               l_row = 0;
               u_row = (static_cast<unsigned int >(as_scalar( ndat_inRec(rc) )))-1;
        } else {
            l_row = u_row+1;
            u_row +=  (static_cast<unsigned int >(as_scalar( ndat_inRec(rc) ))) ;
        }
        One_QRecDATA_IN = get_one_recession(rc);
//     cout <<"\nidhchjhsd1111122222\n";
//cout << l_row << "   " <<u_row<< endl;
        So = as_scalar(Par(rc));
        Par_to_rec = params_RESERVOIR(Par);
        One_QRecDATA_OUT = compute_outflow_1Reccesion(One_QRecDATA_IN, Par_to_rec, So);
//        cout << One_QRecDATA_OUT.n_rows << "   " << One_QRecDATA_OUT.n_cols << endl;
        sim_obs_recession_data.submat(l_row, 0,u_row , 3 ) = One_QRecDATA_OUT;
        One_QRecDATA_IN.reset();
        One_QRecDATA_OUT.reset();
      }
      ofstream out_DATA(FILET.c_str(), ios::app);
      if (!out_DATA) {
       cout << "\nIt is impossible to write to file  " << path_out_file;
       exit(EXIT_FAILURE);
     }
     out_DATA << sim_obs_recession_data;
     out_DATA.close();
}

mat asvh::out_sim_data_MAT(colvec Par)
{
    mat One_QRecDATA_IN;
    mat One_QRecDATA_OUT;
    mat sim_obs_recession_data_DAT;
    colvec Par_to_rec;
    double So=-999999999.9;
    unsigned int u_row, l_row;
    for(unsigned rc=0;rc<Number_of_Recessions;rc++){
//            cout <<"\nidhchjhsd11111\n";
        if(rc ==0){
               l_row = 0;
               u_row = (static_cast<unsigned int >(as_scalar( ndat_inRec(rc) )))-1;
        } else {
            l_row = u_row+1;
            u_row +=  (static_cast<unsigned int >(as_scalar( ndat_inRec(rc) ))) ;
        }
        One_QRecDATA_IN = get_one_recession(rc);
        So = as_scalar(Par(rc));
        Par_to_rec = params_RESERVOIR(Par);
        One_QRecDATA_OUT = compute_outflow_1Reccesion(One_QRecDATA_IN, Par_to_rec, So);
        sim_obs_recession_data.submat(l_row, 0,u_row , 3 ) = One_QRecDATA_OUT;
        One_QRecDATA_IN.reset();
        One_QRecDATA_OUT.reset();
      }
     sim_obs_recession_data_DAT =  sim_obs_recession_data;

    return(sim_obs_recession_data_DAT);
}
