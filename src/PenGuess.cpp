// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;



double ridgepenRcpp(NumericVector x)
{
	int J = x.length();
	int j,k;
	double out=0;
	for (j=0;j<(J-1);j++) 
		for (k=(j+1);k<J;k++) 
			out += pow((x[j]-x[k]),2);
	return out;
}




// [[Rcpp::export]]
arma::mat IRTlikRcppA(arma::vec par, arma::mat data, bool pen, arma::vec nodes, 
                      arma::vec weights, double lambda=0)
{
	// Rcout << "par : " << par << "\n";
	int nitems = par.size()/3;
	int nq = nodes.size();
	arma::vec guess=par.subvec(0,nitems-1);
	arma::vec beta0=par.subvec(nitems,2*nitems-1);
	arma::vec beta1=par.subvec(2*nitems,3*nitems-1);
	
	arma::mat uno_nodes(nq,2);
	uno_nodes.ones();
	uno_nodes.col(1) = nodes;
	
	arma::mat beta01 = join_rows(beta0,beta1);//cbind
	// arma::mat beta01(nitems,2);
	// beta01.col(0)=beta0;
	// beta01.col(1)=beta1;
	
	arma::mat prob = (beta01*uno_nodes.t());
	prob = 1/(1+exp(-prob));
	arma::uvec ids = find(prob >= 1-1.4e-07);
	prob.elem(ids).fill(1-1.4e-07);
	ids = find(prob <= 1.4e-07);
	prob.elem(ids).fill(1.4e-07);

	arma::mat guessrep = arma::repmat(1/(exp(-guess)+1), 1, nq);
	//guessrep.print("guessrep");
	prob = guessrep + (1-guessrep) % prob;
	
	//prob.print("prob");
	arma::mat prodj= exp(data*log(prob)+(1-data)*log(1-prob));
	arma::mat logsumq = log(prodj * weights) ;
	arma::mat lik = sum(logsumq);

	double penalty=0;
	NumericVector guess1 = wrap(guess);
	if (pen) penalty = -lambda * ridgepenRcpp(guess1);
	lik += penalty;
	return lik;
}



// [[Rcpp::export]]
arma::mat gradIRTlikRcppA(arma::vec par, arma::mat data, bool pen, 
                          arma::vec nodes, arma::vec weights, double lambda=0)
{
	int nitems = par.size()/3;
	int nq = nodes.size();
	arma::vec beta3=par.subvec(0,nitems-1);
	arma::vec beta0=par.subvec(nitems,2*nitems-1);
	arma::vec beta1=par.subvec(2*nitems,3*nitems-1);
	
	arma::vec uno(nq);
	uno.ones();
	arma::mat uno_nodes = join_rows(uno,nodes);//cbind

	arma::mat beta01 = join_rows(beta0,beta1);//cbind
	
	arma::mat probs = (beta01*uno_nodes.t());
	probs = 1/(1+exp(-probs));
	arma::uvec ids = find(probs >= 1-1.4e-07);
	probs.elem(ids).fill(1-1.4e-07);
	ids = find(probs <= 1.4e-07);
	probs.elem(ids).fill(1.4e-07);

	arma::mat guessrep = arma::repmat(1/(exp(-beta3)+1), 1, nq);
	arma::mat nodesrep = arma::repmat(nodes,1,nitems);
	arma::mat prob = guessrep + (1-guessrep) % probs;
	
	arma::mat dPdbeta1 = (1-guessrep)%probs%(1-probs);
	arma::mat dPdbeta2 = (1-guessrep)%probs%(1-probs)%nodesrep.t();
	arma::mat dPdbeta3 = guessrep%(1-guessrep)%(1-probs);
	
	arma::mat prodj= exp(data*log(prob)+(1-data)*log(1-prob));
	arma::mat sumq = prodj * weights ;
	arma::mat out(nitems*3,1);
	for (int k=0;k<nitems;k++) {
		arma::mat data_k=data;
		data_k.shed_col(k);
		arma::mat prob_k=prob;
		prob_k.shed_row(k);
		arma::mat prod_k = exp(data_k*log(prob_k)+(1.0 - data_k)*log(1.0 - prob_k));
		arma::mat prod_k1 = prod_k % ((2*data.col(k)-1.0) * dPdbeta1.row(k));
		arma::mat prod_k2 = prod_k % ((2*data.col(k)-1.0) * dPdbeta2.row(k));
		arma::mat prod_k3 = prod_k % ((2*data.col(k)-1.0) * dPdbeta3.row(k));
		out(k) = sum(prod_k3 * weights / sumq);
		out(nitems+k) = sum(prod_k1 * weights / sumq);
		out(2*nitems+k) = sum(prod_k2 * weights / sumq);
	}

	arma::mat grpenalty(nitems,1);
	grpenalty.zeros();
	if (pen) { 
		arma::mat D(nitems,nitems);
		D.fill(-1);
		D.diag().fill((double)nitems - 1.0);
		grpenalty = -lambda*2*D*beta3;
	}
	out.rows(0,nitems-1) += grpenalty;
	return out;
}



