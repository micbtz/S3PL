#define TMB_LIB_INIT R_init_BR3PL
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* data section */
  DATA_MATRIX(y);
  DATA_MATRIX(y1);
  DATA_VECTOR(nodes);
  DATA_VECTOR(weights);  //for GH quadrature
 
  
  /* Parameter section */
  PARAMETER_VECTOR(alpha);  
  PARAMETER_VECTOR(beta); 
  PARAMETER_VECTOR(gamma);
 
  
  using namespace density;

  Type nll=0.0;     // Negative log likelihood function
  int I = y.cols();
  int S = y.rows();
  int Q =  nodes.size();

  
  
  matrix<Type> lprobNOD(I, Q);
  matrix<Type> lprobNOD1(I, Q);
  
  for(int q=0; q<Q; q++)
    {
     Type thetasq = nodes(q) * sqrt(Type(2));
    for(int i=0;i<I;i++)
      {
        Type eta = alpha(i) + thetasq * beta(i);
	Type ci =  invlogit(gamma(i));
	ci = squeeze(ci);
	Type prob =  invlogit(eta) * (Type(1) - ci) + ci;
	prob = squeeze(prob);
	lprobNOD(i, q) = log(prob);
	lprobNOD1(i, q) = log(Type(1) - prob);
        }
    }

  matrix<Type> lprodj(S, Q);
  lprodj = y * lprobNOD + y1 * lprobNOD1; 

  matrix<Type> prodj(S, Q);
  prodj =  exp(lprodj.array()); 

  vector<Type> sqs = prodj * weights;
  vector<Type> lsqs = log(sqs.array());

  nll -= lsqs.sum(); 
  
 
  
  return nll;
}
