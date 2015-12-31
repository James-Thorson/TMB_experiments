
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( c_i );

  // Parameters
  PARAMETER( beta );
  PARAMETER( rho );
  PARAMETER( sigma2 );

  // Probability of random effects
  matrix<Type> Q_yy( c_i.size(), c_i.size() );
  Q_yy.setZero();
  for(int y=0; y<c_i.size(); y++) Q_yy(y,y) = (1+pow(rho,2))/sigma2;
  for(int y=1; y<c_i.size(); y++){
    Q_yy(y-1,y) = -1 * rho/sigma2;
    Q_yy(y,y-1) = 1 * rho/sigma2;
  }
  REPORT( Q_yy )

  // EXAMPLE FROM: http://kaskr.github.io/adcomp/matexp_8cpp_source.html
  //using namespace Eigen;
  //Eigen::EigenSolver< matrix > eigensolver;
  //eigensolver.compute( Q_yy );
  //V = eigensolver.eigenvectors();
  //lambda=eigensolver.eigenvalues();
  //vector<Type> eigenvalues_y(n_y);
  //eigenvalues_y = ;compute().eigenvalues();

  // Elementary
  Eigen::SparseMatrix<Type> Qsparse_yy = asSparseMatrix( Q_yy );
  Type Q_yy_determinant = Q_yy.determinant();
  REPORT( Qsparse_yy );
  REPORT( Q_yy_determinant );

  // LDLT operations
  // https://github.com/kaskr/adcomp/blob/master/TMB/inst/include/Eigen/src/Cholesky/LDLT.h
  matrix<Type> L = Q_yy.ldlt().matrixL();
  vector<Type> diagD = Q_yy.ldlt().vectorD();
  vector<Type> V = Q_yy.jacobiSvd().computeV();
  REPORT( L );
  REPORT( diagD );
  REPORT( V );

  // Eigenvalues
  Eigen::EigenSolver< Eigen::SparseMatrix<Type> > eigensolver;
  //Eigen::EigenSolver< matrix<Type> > eigensolver;
  //eigensolver.compute( Qsparse_yy, true );
  //matrix<Type> V = eigensolver.eigenvectors();
  //vector<Type> lambda = eigensolver.eigenvalues();
  vector< std::complex< Type > > eigenvalues_Q_yy = Q_yy.eigenvalues();
  vector< Type > real_eigenvalues_Q_yy = eigenvalues_Q_yy.real();
  vector< Type > imag_eigenvalues_Q_yy = eigenvalues_Q_yy.imag();
  REPORT( real_eigenvalues_Q_yy );
  REPORT( imag_eigenvalues_Q_yy );
  //std::complex<scalartype>

  // Schur decomposition
  // https://github.com/kaskr/adcomp/blob/master/TMB/inst/include/Eigen/src/Eigenvalues/RealSchur.h
  Eigen::RealSchur< matrix<Type> > realschur;
  //realschur.compute( Q_yy, true );

  // SVD experiments
  Q_yy.jacobiSvd();
  Q_yy.jacobiSvd().computeV();

  // Kroenecker
  matrix<Type> QQ_yy = kronecker(Q_yy, Q_yy);
  Eigen::SparseMatrix<Type> QQsparse_yy = asSparseMatrix( QQ_yy );
  REPORT( QQ_yy );
  REPORT( QQsparse_yy );
  //REPORT( QQsparse_yy );

  // Probability of data conditional on random effects
  Type Total_Abundance = 0;
  Type jnll = 0;
  for( int i=0; i<c_i.size(); i++){
    jnll -= dnorm( c_i(i), beta, Type(1.0), true );
  }

  return jnll;
}
