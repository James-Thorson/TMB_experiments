
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace R_inla;

  // Options
  DATA_IVECTOR( Options_z );

  // Data
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_FACTOR( j_i );  // Random effect index for observation i

  // SPDE objects
  DATA_STRUCT(spde, spde_t);

  // 2D AR1 objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );

  // Random effects
  PARAMETER_VECTOR( epsilon_j );

  // Objective funcction
  int n_i = c_i.size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();

  // Derived quantities

  // Precision matrix
  Eigen::SparseMatrix<Type> Q_matern = Q_spde(spde, exp(ln_kappa) );
  Eigen::SparseMatrix<Type> Q_ar = (M0*pow(1+exp(ln_kappa*2),2) + M1*(1+exp(ln_kappa*2))*(-exp(ln_kappa)) + M2*exp(ln_kappa*2));
  REPORT( Q_matern );
  REPORT( Q_ar );

  // Probability of random effects
  Type Range;
  Type MargSD;
  if( Options_z(0)==0 ){
    jnll_comp(1) += SCALE( GMRF(Q_matern), 1/exp(ln_tau) )( epsilon_j );
    Range = sqrt(8) / exp( ln_kappa );
    MargSD = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));
  }
  if( Options_z(0)==1 ){
    jnll_comp(1) += SCALE( GMRF(Q_ar), exp(ln_tau) )( epsilon_j );
    Range = log(0.1) / ln_kappa;
    MargSD = exp(ln_tau) / sqrt(1-exp(ln_kappa*2));
  }
  REPORT( Range );
  REPORT( MargSD );

  // Probability of data conditional on random effects
  for( int i=0; i<n_i; i++){
    if( !isNA(c_i(i)) ) jnll_comp(0) -= dpois( c_i(i), exp(beta0 + epsilon_j(j_i(i))), true );
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );

  return jnll;
}
