#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Randomizer{
public:
    Randomizer(unsigned long int seed = 0);
    ~Randomizer();
    void Reset();

    unsigned int Binomial(unsigned int n, double p);
    int UniformInt(int min, int max);
    //inline gsl_rng* GSL_RNG() { return r; }
private:
    unsigned long int seed;
    gsl_rng* r;
};

Randomizer::Randomizer(unsigned long int s)
 : seed(s), r(gsl_rng_alloc(gsl_rng_mt19937))
{
    Reset();
}

Randomizer::~Randomizer(){
    gsl_rng_free(r);
}

void Randomizer::Reset(){
    gsl_rng_set(r, seed);
}

unsigned int Randomizer::Binomial(unsigned int n, double p){
    if (p <= 0) return 0;
    return gsl_ran_binomial(r, p, n);
}

int Randomizer::UniformInt(int min, int max){
    return min + gsl_rng_uniform(r) * (max - min);
}
