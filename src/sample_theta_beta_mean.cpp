#include <RcppArmadillo.h>
#include <math.h>
#include <shrinkTVP.h>
using namespace Rcpp;

void sample_theta_beta_mean(arma::vec& theta,
                            arma::vec& beta_mean,
                            const arma::mat& beta,
                            double tau2,
                            double s0,
                            double S0){
  // Get dimension
  int d = theta.n_elem;
  int N = beta.n_cols - 1; //-1 to correct for 0th state

  arma::vec beta_mean_new(d);
  arma::vec theta_new(d);


  for (int j = 0; j < d; j++){
    double sigma2_beta_mean = 1.0/(1.0/tau2 + 1.0/(theta(j)));
    double mu_beta_mean = beta(j, 0) * tau2/(tau2 + theta(j));
    beta_mean_new(j) = R::rnorm(mu_beta_mean, std::sqrt(sigma2_beta_mean));
  }

  for (int j = 0; j < d; j++){
    double s0_post = (N + 2.0*s0 + 1.0)*0.5;
    double S0_post = .5*(std::pow(beta(j, 0) - beta_mean(j), 2) +
                         arma::as_scalar(arma::accu(arma::pow(arma::diff(beta.row(j)), 2)))) +
                         S0;
    theta_new(j) = S0_post/R::rgamma(s0_post, 1);
  }

  std::for_each(theta_new.begin(), theta_new.end(), shrinkTVP::res_protector);
  std::for_each(beta_mean_new.begin(), beta_mean_new.end(), shrinkTVP::res_protector);

  beta_mean = beta_mean_new;
  theta = theta_new;
}
