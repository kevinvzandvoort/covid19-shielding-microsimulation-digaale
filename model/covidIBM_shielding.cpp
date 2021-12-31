// covidIBM.cpp

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppGSL)]]

#include <iostream>
#include <vector>
#include <cstdlib>
#include "./randomizer.h"
#include <Rcpp.h>
#include <omp.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "./ibm.h"

int main(){}

Rcpp::List RunSimulation(Rcpp::List parameters, int seed){
  
  // load parameter values from parameters list
  Rcpp::DataFrame population_data = parameters["population"];
  int agegroups = parameters["n_agegroups"];
  int Iclasses = 3;
  int n_households = parameters["n_households"];
  int initial_inf = parameters["n_init_infected"];
  int times = parameters["days"];
  int tstep = parameters["tstep_day"];
  double beta_hh = parameters["contact_household"];
  Rcpp::NumericMatrix beta = parameters["contact_matrix"];
  Rcpp::List interventions = parameters["interventions"];
  
  double test_specificity = parameters["test_specificity"];
  //Can be indexed by Individual::Status - 1
  double test_sensitivity[4] = {
    parameters["test_sensitivity_preinfectious"], parameters["test_sensitivity_subclinical"],
    parameters["test_sensitivity_preclinical"], parameters["test_sensitivity_clinical"]
  };
  
  int shield_age = 99;
  
  double contacts_shielded_shielded[2] = { 1.0, 1.0 };
  double contacts_shielded_unshielded = 1.0;
  double contacts_unshielded_shielded = 1.0;
  double contacts_unshielded_unshielded[2] = { 1.0, 1.0 };
  double contacts_shielded_infectious_clinical[2] = { 1.0, 1.0 };
  double contacts_unshielded_infectious_clinical[2] = { 1.0, 1.0 };
  
  int shielding_clinical_exit_delay = 999;
  bool shielding_clinical_exit_all = false;
  
  /*
    allocate memory for household- and individual-pointer vectors
  use pointers to minimize overhead
  */
  std::vector<Household*> households;
  std::vector<Household*> greenzones;
  households.reserve(n_households);
  std::vector<Individual*> individuals;
  int n_pop = population_data.nrows();
  individuals.reserve(n_pop);
  
  // use temporary scope to load data
  {
    Rcpp::NumericVector df_part_id = population_data["participant_id"];
    Rcpp::NumericVector df_household_id = population_data["household_id"];
    Rcpp::NumericVector df_age_group = population_data["age_group"];
    Rcpp::NumericVector df_dE = population_data["dE"];
    Rcpp::NumericVector df_dP = population_data["dP"];
    Rcpp::NumericVector df_dC = population_data["dC"];
    Rcpp::NumericVector df_dS = population_data["dS"];
    Rcpp::NumericVector df_y = population_data["y"];
    Rcpp::NumericVector df_u = population_data["u"];
    Rcpp::NumericVector df_immune = population_data["immune"];
    
    int household_id = -1;
    int h = -1;
    for(int i=0; i<n_pop; i++){
      // add household pointer to household vector
      if(df_household_id[i] != household_id){
        household_id = df_household_id[i];
        h++;
        households.emplace_back(new Household(household_id, 0));
      }
      
      // add individual pointer to household and individual vector
      Individual* indiv = new Individual(df_part_id[i], household_id, df_age_group[i], households[h], df_dE[i], df_dP[i], df_dC[i], df_dS[i], df_y[i], df_u[i], df_immune[i]);
      individuals.emplace_back(indiv);
      households[h]->addMember(indiv);
    }
  }
  
  //allocate three-dimensional array (age-group, I-class, shielded-status)
  int infectious_strata[agegroups][Iclasses][2];
  for(int a=0; a<agegroups; a++){
    for(int i=0; i<Iclasses; i++){
      for(int s=0; s<2; s++){
        infectious_strata[a][i][s] = 0; 
      }
    }
  }
  
  // vector to store infected individuals (pointers)
  std::vector<Individual*> susceptibles;
  std::vector<Individual*> infectious;
  
  // pre-allocate sufficient memory so addresses never change
  susceptibles.reserve(n_pop);
  infectious.reserve(n_pop);
  
  // ensure random number generator works when using multiple threads
  std::vector<Randomizer*> randomizers;
  randomizers.reserve(omp_get_max_threads()+1);
  
  for(int p=0; p <= omp_get_max_threads(); p++){
    Randomizer* rndm = new Randomizer(seed+p);
    randomizers.emplace_back(rndm);
  }
  
  // random number generator used in the main thread
  Randomizer* random_generator = randomizers[omp_get_max_threads()];
  
  // randomly inocculate individual(s)
  while(infectious.size() < initial_inf){
    int rnd_infected = random_generator->UniformInt(0, n_pop-1);
    if(individuals[rnd_infected]->m_status == Individual::susceptible){
      if(random_generator->Binomial(1, individuals[rnd_infected]->m_y) == 1){
        individuals[rnd_infected]->m_status=Individual::infected_preclinical;
      } else {
        individuals[rnd_infected]->m_status=Individual::infected_subclinical;
      }
      individuals[rnd_infected]->m_infected_at = 0.0;
      individuals[rnd_infected]->m_infectious_at = 0.0;
      individuals[rnd_infected]->m_infected_by = 0;
      infectious.emplace_back(individuals[rnd_infected]);
      infectious_strata[individuals[rnd_infected]->m_age][individuals[rnd_infected]->m_status - 2][!individuals[rnd_infected]->m_shielded_status] += 1; 
    }
  }
  
  // loop through timesteps
  for(int t = 1; t <= times*tstep; t++){
    double t_day = t/tstep;
    
    //reset vectors
    susceptibles.clear();
    susceptibles.reserve(n_pop);
    infectious.clear();
    infectious.reserve(n_pop);
    
    //reset total count of infectious strata
    for(int a=0; a<agegroups; a++){
      for(int i=0; i<Iclasses; i++){
        for(int s=0; s<2; s++){
          infectious_strata[a][i][s] = 0; 
        }
      }
    }
    
    int preinfectious = 0;
    //update status of all infected individuals and regenerate susceptible, infectious vectors
    for(int i = 0; i < individuals.size(); ++i){
      Individual* member = individuals[i];
      
      //progress to next state if sufficient time has passed
      bool change_status = member->progressStatus(random_generator, t_day, shielding_clinical_exit_all);
      
      if( member->m_status == Individual::susceptible){
        susceptibles.emplace_back(member);
      } else if(
          member->m_status == Individual::infected_preclinical |
            member->m_status == Individual::infected_clinical |
            member->m_status == Individual::infected_subclinical
      ){
        infectious.emplace_back(member);
        infectious_strata[member->m_age][member->m_status - 2][!member->m_shielded_status] += 1; 
      } else if(member->m_status == Individual::preinfectious){
        preinfectious += 1;
      }
    }
    
    //stop simulation if infection is cleared from population
    if(preinfectious == 0 & infectious.size() == 0){
      break;
    }
    
    //process interventions
    for(int i=0; i < interventions.size(); i++){
      Rcpp::List intervention = interventions[i];
      int iv_processed = intervention["processed"];
      
      //only continue if intervention is not processed before
      if(iv_processed == 1){
        continue; 
      }
      
      //check if intervention should start now
      bool iv_start = false;
      std::string iv_start_event = intervention["start_event"];
      double iv_start_value = intervention["start_value"];
      if(iv_start_event == "time"){
        if(t_day >= iv_start_value){
          iv_start = true;
        }
      } else if(iv_start_event == "prevalence"){
        //prevalence of cases (clinical infectious people)
        int cases = 0;
        for(int a=0; a<agegroups; a++){
          cases += infectious_strata[a][2][0] + infectious_strata[a][2][1];
        }
        
        double prevalence = (double)cases/n_pop;
        if(prevalence >= iv_start_value){
          iv_start = true;
        }
      }
      
      if(!iv_start){
        continue;
      }
      
      intervention["start_time"] = t;
      
      if( intervention.containsElementNamed("shielding_group_size") ){
        int max_size = intervention["shielding_group_size"];
        double shield_probability = intervention["shielding_proportion"];
        int shielding_enter_control = intervention["shielding_enter_control"];
        std::string shielding_mitigate_type = intervention["shielding_mitigate_type"];
        int shielding_mitigate_value = intervention["shielding_mitigate_value"];
        shield_age = intervention["shielding_age_min"];
        
        if(shielding_mitigate_type == "delay"){
          shielding_clinical_exit_delay = shielding_mitigate_value;
        }
        
        if(shielding_mitigate_type == "allexit"){
          bool shielding_clinical_exit_all = true;
        }
        
        int g = 0;
        
        greenzones.reserve(n_households);
        greenzones.emplace_back(new Household(g, 0));
        for(int s=0; s < individuals.size(); s++){
          Individual* member = individuals[s];
          if(member->m_age >= shield_age){
            if( random_generator->Binomial(1, shield_probability) == 1 ){
              if(shielding_enter_control == 2 & member->m_status == Individual::infected_clinical){
                //individuals with clinical symptoms are not shielded
                continue;
              } else {
                if(shielding_enter_control == 3){
                  double test_prob = 0.0;
                  switch(member->m_status){
                  case Individual::preinfectious:
                    test_prob = test_sensitivity[0];
                    break;
                  case Individual::infected_subclinical:
                    test_prob = test_sensitivity[1];
                    break;
                  case Individual::infected_preclinical:
                    test_prob = test_sensitivity[2];
                    break;
                  case Individual::infected_clinical:
                    test_prob = test_sensitivity[3];
                    break;
                  default:
                    test_prob = 1.0 - test_specificity;
                  }
                  if( random_generator->Binomial(1, test_prob) == 1 ){
                    //individuals with a positive test are not shielded
                    continue;
                  }
                }
              }
              member->m_shielded_status = true;
              member->m_duration_shielded = shielding_clinical_exit_delay*tstep;
              member->m_shielded_gz_id = g;
              member->m_shielded_gz_ptr = greenzones[g];
              greenzones[g]->addMember(member);
              
              if(greenzones[g]->m_household_size >= max_size){
                g++;
                greenzones.emplace_back(new Household(g, 0));
              }
            }
          }
        }
      }
      if( intervention.containsElementNamed("contacts_shielded_shielded") ){
        Rcpp::NumericVector iv_contacts_shielded_shielded = intervention["contacts_shielded_shielded"];
        contacts_shielded_shielded[0] *= iv_contacts_shielded_shielded[0];
        contacts_shielded_shielded[1] *= iv_contacts_shielded_shielded[1];
      }
      if( intervention.containsElementNamed("contacts_shielded_unshielded") ){
        double iv_contacts_shielded_unshielded = intervention["contacts_shielded_unshielded"];
        contacts_shielded_unshielded *= iv_contacts_shielded_unshielded;
        contacts_unshielded_shielded *= iv_contacts_shielded_unshielded;
      }
      if( intervention.containsElementNamed("contacts_unshielded_unshielded") ){
        Rcpp::NumericVector iv_contacts_unshielded_unshielded = intervention["contacts_unshielded_unshielded"];
        contacts_unshielded_unshielded[0] *= iv_contacts_unshielded_unshielded[0];
        contacts_unshielded_unshielded[1] *= iv_contacts_unshielded_unshielded[1];
      }
      if( intervention.containsElementNamed("contacts_shielded_infectious_clinical") ){
        Rcpp::NumericVector iv_contacts_shielded_infectious_clinical = intervention["contacts_shielded_infectious_clinical"];
        contacts_shielded_infectious_clinical[0] *= iv_contacts_shielded_infectious_clinical[0];
        contacts_shielded_infectious_clinical[1] *= iv_contacts_shielded_infectious_clinical[1];
      }
      if( intervention.containsElementNamed("contacts_unshielded_infectious_clinical") ){
        Rcpp::NumericVector iv_contacts_unshielded_infectious_clinical = intervention["contacts_unshielded_infectious_clinical"];
        contacts_unshielded_infectious_clinical[0] *= iv_contacts_unshielded_infectious_clinical[0];
        contacts_unshielded_infectious_clinical[1] *= iv_contacts_unshielded_infectious_clinical[1];
      }
      intervention["processed"] = 1;
    }
    
    // check for contact with all susceptible individuals
    #pragma omp parallel for
    for(int s = 0; s < susceptibles.size(); s++){
      Randomizer* random_generator_thread = randomizers[omp_get_thread_num()];
      Individual* member = susceptibles[s];
      
      int infectious_strata_household[agegroups][Iclasses][2];
      for(int a=0; a<agegroups; a++){
        for(int i=0; i<Iclasses; i++){
          for(int z=0; z<2; z++){
            infectious_strata_household[a][i][z] = 0; 
          }
        }
      }
      
      //are any household/greenzone members infectious?
      int hhmem_inf = 0;
      if(member->m_shielded_status){
        for(int h=0; h < member->m_shielded_gz_ptr->members.size(); h++){
          Individual* gzmember = member->m_shielded_gz_ptr->members[h];
          if(gzmember->m_status == Individual::infected_preclinical |
             gzmember->m_status == Individual::infected_clinical |
             gzmember->m_status == Individual::infected_subclinical
          ){
            //all gzmembers will have m_shielded_status == true
            infectious_strata_household[gzmember->m_age][gzmember->m_status - 2][!gzmember->m_shielded_status] += 1;
            hhmem_inf += 1;
          }
        }
      } else {
        for(int h=0; h < member->m_household_ptr->members.size(); h++){
          Individual* hhmember = member->m_household_ptr->members[h];
          if(hhmember->m_status == Individual::infected_preclinical |
             hhmember->m_status == Individual::infected_clinical |
             hhmember->m_status == Individual::infected_subclinical
          ){
            //hhmember m_shielded_status could be true or false
            infectious_strata_household[hhmember->m_age][hhmember->m_status - 2][!hhmember->m_shielded_status] += 1;
            hhmem_inf += 1;
          }
        } 
      }
      
      bool transmission = false;
      for(int a = 0; a < agegroups; a++){
        if(transmission){ break; }
        
        for(int i = 0; i < Iclasses; i++){
          if(transmission){ break; }
          
          if(member->m_shielded_status){
            //transmission with people in same bubble
            int Nhh = infectious_strata_household[a][i][0];
            transmission = member->transmissionQuick(a, i, true, Nhh, beta, beta_hh, t_day, random_generator_thread, contacts_shielded_shielded[0] * (i == 2 ? contacts_shielded_infectious_clinical[0] : 1) );
            if(!transmission){
              //transmission with shielded people in other bubbles (other-gz)
              int Ncom_shielded = infectious_strata[a][i][0];
              transmission = member->transmissionQuick(a, i, false, Ncom_shielded - Nhh, beta, beta_hh, t_day, random_generator_thread, contacts_shielded_shielded[1] * (i == 2 ? contacts_shielded_infectious_clinical[1] : 1)); 
            }
            if(!transmission){
              //transmission with unshielded
              int Ncom_unshielded = infectious_strata[a][i][1];
              transmission = member->transmissionQuick(a, i, false, Ncom_unshielded, beta, beta_hh, t_day, random_generator_thread, contacts_shielded_unshielded * (i == 2 ? contacts_unshielded_infectious_clinical[1] : 1));
            }
          } else {
            //transmission with people in same household (not shielded)
            int Nhh = infectious_strata_household[a][i][1];
            transmission = member->transmissionQuick(a, i, true, Nhh, beta, beta_hh, t_day, random_generator_thread, contacts_unshielded_unshielded[0] * (i == 2 ? contacts_unshielded_infectious_clinical[0] : 1));
            if(!transmission){
              //transmission with shielded people
              int Ncom_shielded = infectious_strata[a][i][0];
              transmission = member->transmissionQuick(a, i, false, Ncom_shielded, beta, beta_hh, t_day, random_generator_thread, contacts_unshielded_shielded * (i == 2 ? contacts_shielded_infectious_clinical[1] : 1)); 
            }
            if(!transmission){
              //transmission with unshielded (non-hh)
              int Ncom_unshielded = infectious_strata[a][i][1];
              transmission = member->transmissionQuick(a, i, false, Ncom_unshielded - Nhh, beta, beta_hh, t_day, random_generator_thread, contacts_unshielded_unshielded[1] * (i == 2 ? contacts_unshielded_infectious_clinical[1] : 1));
            }
          }
        }
      }
    }
  }
  
  //std::cout << "Finished simulation" << std::endl;
  //std::cout << "Data export..." << std::endl;
  
  // Data export
  // set up vectors for R
  std::vector<int> df_part_id;
  df_part_id.reserve(n_pop);
  std::vector<int> df_household_id;
  df_household_id.reserve(n_pop);
  std::vector<int> df_age;
  df_age.reserve(n_pop);
  std::vector<int> df_status;
  df_status.reserve(n_pop);
  std::vector<int> df_shielded_status;
  std::vector<int> df_shielded_gz_id;
  //Rcpp::NumericVector df_infected_by;
  std::vector<int> df_infected_at;
  df_infected_at.reserve(n_pop);
  
  std::vector<int> df_infected_shielded;
  std::vector<int> df_infectious_at;
  df_infectious_at.reserve(n_pop);
  std::vector<int> df_recovered_at;
  df_recovered_at.reserve(n_pop);
  
  // add individual data to vectors
  for(int i=0; i<individuals.size(); i++){
    int status = individuals[i]->m_status;
    
    //do not export susceptible
    if(status > 0){
      df_part_id.emplace_back(individuals[i]->m_part_id);
      df_household_id.emplace_back(individuals[i]->m_household_id);
      df_age.emplace_back(individuals[i]->m_age);
      df_status.emplace_back(status);
      df_shielded_status.emplace_back(individuals[i]->m_shielded_status);
      df_infected_shielded.emplace_back(individuals[i]->m_infected_shielded);
      df_shielded_gz_id.emplace_back(individuals[i]->m_shielded_gz_id);
      
      //df_infected_by.push_back(individuals[i]->m_infected_by);
      df_infected_at.emplace_back(individuals[i]->m_infected_at);
      df_infectious_at.emplace_back(individuals[i]->m_infectious_at);
      if(status == 5){
        df_recovered_at.emplace_back(individuals[i]->m_recovered_at);
      } else {
        df_recovered_at.emplace_back(-1.0);
      }
    }
    
    delete individuals[i];
  }
  
  //std::cout << "Prepare data for R..." << std::endl;
  
  // make sure to copy
  Rcpp::DataFrame population = Rcpp::DataFrame::create(
    Rcpp::Named("part_id") = df_part_id,
    Rcpp::Named("household_id") = df_household_id,
    Rcpp::Named("age") = df_age,
    Rcpp::Named("status") = df_status,
    Rcpp::Named("shielded_status") = df_shielded_status,
    Rcpp::Named("infected_shielded") = df_infected_shielded,
    Rcpp::Named("shielded_gz_id") = df_shielded_gz_id,
    //Rcpp::Named("infected_by") = Rcpp::clone(df_infected_by),
    Rcpp::Named("infected_at") = df_infected_at,
    Rcpp::Named("infectious_at") = df_infectious_at,
    Rcpp::Named("recovered_at") = df_recovered_at
  );
  
  //std::cout << "Process interventions... " << interventions.size() << std::endl;
  
  Rcpp::List intervention_data(interventions.size());
  for(int i=0; i < interventions.size(); i++){
    Rcpp::List intervention = interventions[i];
    Rcpp::List ivlist = Rcpp::List::create(
      Rcpp::Named("processed", intervention["processed"]),
      Rcpp::Named("start_time", intervention["start_time"])
    );
    intervention_data[i] = ivlist;
  }
  
  //std::cout << "Combine list..." << std::endl;
  
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("intervention") = intervention_data,
    Rcpp::Named("population") = population
  );
  
  for(int h=0; h<n_households; h++){
    delete households[h];
  }
  
  //std::cout << "Done" << std::endl;
  
  return(output);
}

// [[Rcpp::export]]
Rcpp::List covidIBM(Rcpp::List parameters, unsigned long int seed = 0)
{
  
  // Initialise parameters for this simulation
  Rcpp::List output;
  
  output = RunSimulation(parameters, seed);
  
  return output;
}
