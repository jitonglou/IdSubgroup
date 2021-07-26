#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

// --------------------
  // Functions
// --------------------

// [[Rcpp::export]]
arma::mat MDist_C(const int& N, const arma::mat& A, arma::mat& B) {
  arma::mat out = arma::zeros<arma::mat>(N,N);
  for(int i = 0; i < N-1; ++i) {
    for(int j = i+1; j < N; ++j) {
      out.at(i,j) = arma::as_scalar((A.row(i) - A.row(j)) * B * (A.row(i).t() - A.row(j).t()));
    }
  }
  return out;
} // calcute Mahalanobis distance between each subject

// [[Rcpp::export]]
double logSumExp_C(const arma::vec& log_x){
	if( log_x.n_elem == 1 ){
		return log_x.at(0);
	} else {
		double max_val = max(log_x);
		arma::vec log_x_2 = log_x - max_val;
		return log(arma::sum(arma::exp(log_x_2))) + max_val;
	}
} // calcute extremely small values

// [[Rcpp::export]]
double logy_norm_C(const double& y, const double& lf, const double& par2) {
  return -std::pow(y-lf,2.0)/(2.0*par2);
} // calculate log-likelihood of normal distributed variable


// [[Rcpp::export]]
double logy_bern_C(const double& y, const double& lf, const double& par2) {
  double par = std::exp(lf)/(1.0+std::exp(lf));
  return y*std::log(par) + (1.0-y)*std::log(1.0-par);
} // calculate log-likelihood of Bernoulli distributed variable

// [[Rcpp::export]]
double logy_pois_C(const double& y, const double& lf, const double& par2) {
  double par = std::exp(lf);
  return y*std::log(par) - par - std::lgamma(y+1.0);
} // calculate log-likelihood of Poisson distributed variable






