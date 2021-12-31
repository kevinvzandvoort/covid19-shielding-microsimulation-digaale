#ifdef MODE_DEBUG
#define DEBUG(x) { x }
#else
#define DEBUG(x)
#endif

//forward declarations
class Household;
class Individual;
//Randomizer random_generator;
Rcpp::NumericMatrix beta;
double beta_hh;
Rcpp::NumericVector u;

class Household{
public:
  const unsigned int m_household_id;
  unsigned int m_household_size;
  std::vector<Individual*> members;
public:
  Household(int household_id, int household_size)
    : m_household_id(household_id), m_household_size(household_size)
  {
  }
  Household(const Household& other)
    : m_household_id(other.m_household_id), m_household_size(other.m_household_size)
  {
    members.reserve(m_household_size);
    for(int i=0; i<other.m_household_size; i++)
      members.emplace_back(other.members[i]);
  }
  ~Household(){
  }
  void addMember(Individual* indiv){
    members.emplace_back(indiv);
    m_household_size++;
  }
};

class Individual{
public:
  enum Status{
    susceptible = 0,
    preinfectious = 1,
    infected_subclinical = 2,
    infected_preclinical = 3,
    infected_clinical = 4,
    recovered = 5
  };
  const unsigned int m_part_id;
  const unsigned int m_household_id;
  const unsigned int m_age;
  Status m_status = susceptible;
  bool m_shielded_status;
  Household* m_shielded_gz_ptr;
  int m_shielded_gz_id;
  double m_y;
  double m_u;
  double m_duration_E;
  double m_duration_Is;
  double m_duration_Ip;
  double m_duration_Ic;
  unsigned int m_infected_by;
  double m_infected_at;
  bool m_infected_shielded;
  int m_duration_shielded;
  double m_infectious_at;
  double m_recovered_at;
  Household* m_household_ptr;
  
public:
  //require IDs and age for constructor
  Individual(int part_id, int household_id, int age, Household* household_ptr, double dE, double dP, double dC,
             double dS, double y, double u, int immune) :
             m_part_id(part_id), m_household_id(household_id), m_age(age),
             m_status(immune == 1 ? recovered : susceptible), m_shielded_status(false), m_shielded_gz_id(-1),
             m_household_ptr(household_ptr), m_duration_E(dE), m_duration_Ip(dP), m_duration_Ic(dC), m_duration_Is(dS),
             m_y(y), m_u(u), m_duration_shielded(-1), m_infected_shielded(false)
  {
    if(immune == 1){
      m_infected_at = -1.0;
      m_infectious_at = -1.0;
      m_recovered_at = -1.0;
      m_infected_by = 0;
    }
  }
  //check for copy to optimize code
  Individual(const Individual& other)
    : m_part_id(other.m_part_id), m_household_id(other.m_household_id), m_age(other.m_age), m_status(other.m_status),
      m_infected_by(other.m_infected_by), m_infected_at(other.m_infected_at)
  {
  }
  ~Individual()
  {
  }
  bool progressStatus(Randomizer* random_generator, double time, bool &shielding_clinical_exit_all){
    bool change_status = false;
    switch(m_status){
    case preinfectious:
      if(m_duration_E == 0){
        //update status if sufficient time has passed
        //m_y determines whether individual becomes preclinical or subclinical
        if(random_generator->Binomial(1, m_y) == 1){
          m_status = infected_preclinical;
        } else {
          m_status = infected_subclinical;
        }
        change_status = true;
      } else {
        m_duration_E -= 1;
      }
      break;
    case infected_subclinical:
      if(m_duration_Is == 0){ m_status = recovered; change_status = true; m_recovered_at = time; } else { m_duration_Is -= 1; }
      break;
    case infected_preclinical:
      if(m_duration_Ip == 0){ m_status = infected_clinical; change_status = true; } else { m_duration_Ip -= 1; }
      break;
    case infected_clinical:
      if(m_shielded_status == true){
        //send all shielded home after one develops symptoms
        if(shielding_clinical_exit_all){
          for(int i=0; i < (m_shielded_gz_ptr->m_household_size); i++){
            m_shielded_gz_ptr->members[i]->m_shielded_status = false;
          }
        }
        //remove shielded individuals after delay
        if(m_duration_shielded == 0){ m_shielded_status = false; } else { m_duration_shielded -= 1; }
      }
      if(m_duration_Ic == 0){ m_status = recovered; change_status = true; m_recovered_at = time; } else { m_duration_Ic -= 1; }
      break;
    default:
      change_status = false;
    }
    if(change_status & m_status > 1 & m_status < 4){
      m_infectious_at = time;
    }
    return change_status;
  }
  
  bool transmission(Individual* infected_contact, int hh_index, int gz_index, Rcpp::NumericMatrix &beta, double &beta_hh, double time, Randomizer* random_generator, double contacts_shielded_shielded[2], double contacts_shielded_unshielded[2], double contacts_unshielded_unshielded[2], double infected_contact_adapt_transmission){
    //double t=time;
    if(m_status != susceptible){
      return false;
    } else {
      // subclinical infections are 50% less infectious
      double f = 1.0;
      if(infected_contact->m_status == infected_subclinical){
        f = 0.5;
      }
      //if hh_index or gz_index are 0, contacts are household and/or greenzone members, and homogeneous mixing is used
      bool use_homog_mix = (hh_index+gz_index < 2);
      double lambda = m_u * f * (use_homog_mix ? beta_hh : (beta(infected_contact->m_age,m_age)));
      double member_adapt_transmission = 1.0;
      if(m_shielded_status == true){
        if(infected_contact->m_shielded_status == true){
          //member and contact are shielded
          if(infected_contact->m_shielded_gz_id == m_shielded_gz_id){
            gz_index = 0;
          }
          member_adapt_transmission *= contacts_shielded_shielded[gz_index];
        } else {
          //member is shielded, contact is not shielded
          member_adapt_transmission *= contacts_shielded_unshielded[hh_index];
        }
      } else {
        //member is not shielded
        member_adapt_transmission *= contacts_unshielded_unshielded[hh_index];
      }
      
      lambda = lambda * member_adapt_transmission * infected_contact_adapt_transmission;
      bool is_infected = random_generator->Binomial(1, lambda);
      if(is_infected){
        m_status = preinfectious;
        m_infected_at = time;
        m_infected_by = infected_contact->m_part_id;
        if(m_shielded_status == true){ m_infected_shielded = true; }
      }
      return is_infected;
    }
  }
  
  bool transmissionQuick(int &contacts_age, int &I_class, bool hhmember, int N, Rcpp::NumericMatrix &beta, double &beta_hh, double time, Randomizer* random_generator, double intervention_effect){
    //double t=time;
    if(m_status != susceptible){
      return false;
    } else if(N == 0){
      return false;
    } else {
      // subclinical infections are 50% less infectious
      double f = 1.0;
      if(I_class == 0){
        f = 0.5;
      }
      
      double lambda = m_u * f * (hhmember ? beta_hh : (beta(contacts_age, m_age))) * intervention_effect;
      /*if(time == 1 & N > 1){
        std::cout << "N is " << N << " and lambda is " << lambda << " and lambda_p is " << 1-std::pow((1-lambda),N) << std::endl;
      }*/
      bool is_infected = random_generator->Binomial(1, 1 - std::pow((1 - lambda), N));
      if(is_infected){
        m_status = preinfectious;
        m_infected_at = time;
        //m_infected_by = infected_contact->m_part_id;
        if(m_shielded_status == true){ m_infected_shielded = true; }
      }
      return is_infected;
    }
  }

};
