#include <RcppArmadillo.h>
#include <stochvol.h>
#include <shrinkTVP.h>
#include <progress.hpp>
#include <math.h>
#include "sample_theta_beta_mean.h"
using namespace Rcpp;


//[[Rcpp::export]]
List BVSTVP_cpp(arma::vec y,
                arma::mat x,
                int niter,
                int nburn,
                int nthin,
                double c0,
                double g0,
                double G0,
                double tau2,
                double s0,
                double S0,
                bool display_progress,
                bool sv,
                double Bsigma_sv,
                double a0_sv,
                double b0_sv,
                double bmu,
                double Bmu,
                Rcpp::List starting_vals) {

  // Progress bar setup
  arma::vec prog_rep_points = arma::round(arma::linspace(0, niter, 50));
  Progress p(50, display_progress);

  // Define some necessary dimensions
  int N = y.n_elem;
  int d = x.n_cols;
  int nsave = std::floor((niter - nburn)/nthin);

  // Storage objects (returned to R at the end)
  arma::cube beta_save(N+1, d, nsave, arma::fill::zeros);
  arma::cube sigma2_save(N,1, nsave, arma::fill::zeros);
  arma::mat theta_save(d, nsave, arma::fill::zeros);
  arma::mat beta_mean_save(d, nsave, arma::fill::zeros);

  arma::vec C0_save;
  arma::vec sv_mu_save;
  arma::vec sv_phi_save;
  arma::vec sv_sigma2_save;

  if (sv == false) {
    C0_save = arma::vec(nsave, arma::fill::zeros);
  } else {
    sv_mu_save = arma::vec(nsave, arma::fill::zeros);
    sv_phi_save = arma::vec(nsave, arma::fill::zeros);
    sv_sigma2_save = arma::vec(nsave, arma::fill::zeros);
  }

  // Initial values and objects holding current values of samples
  // Parameters where learning can be toggled are also overwritten if learning is not activated
  arma::mat beta_samp(d, N+1, arma::fill::zeros);
  arma::mat beta_nc_samp(d, N+1, arma::fill::zeros);

  arma::vec beta_mean_samp(d);
  beta_mean_samp = as<arma::vec>(starting_vals["beta_mean_st"]);

  arma::vec theta_samp(d);
  arma::vec theta_sr_samp(d);
  theta_samp = as<arma::vec>(starting_vals["theta_st"]);
  theta_sr_samp = arma::sqrt(theta_samp);

  arma::vec sigma2_samp;
  double C0_samp;
  arma::vec sv_para(3);
  if (sv == true) {
    sigma2_samp = as<arma::vec>(starting_vals["sigma2_st"]);
    sv_para = {as<double>(starting_vals["sv_mu_st"]),
               as<double>(starting_vals["sv_phi_st"]),
               as<double>(starting_vals["sv_sigma2_st"])};
  } else {
    sigma2_samp.copy_size(y);
    sigma2_samp.fill(as<double>(starting_vals["sigma2_st"]));
    C0_samp = as<double>(starting_vals["C0_st"]);
  }

  // Objects required for stochvol to work
  arma::uvec r(N); r.fill(5);
  double h0_samp = as<double>(starting_vals["h0_st"]);
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {  // prior specification object for the update_*_sv functions
    PriorSpec::Latent0(),  // stationary prior distribution on priorlatent0
    PriorSpec::Mu(PriorSpec::Normal(bmu, std::sqrt(Bmu))),  // normal prior on mu
    PriorSpec::Phi(PriorSpec::Beta(a0_sv, b0_sv)),  // stretched beta prior on phi
    PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma_sv))  // normal(0, Bsigma) prior on sigma
  };  // heavy-tailed, leverage, regression turned off
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert {  // very expert settings for the Kastner, Fruehwirth-Schnatter (2014) sampler
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    2,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };


  // Values for LPDS calculation
  arma::cube m_N_save(d, 1, nsave);
  arma::cube chol_C_N_inv_save(d, d, nsave);
  arma::vec m_N_samp;
  arma::mat chol_C_N_inv_samp;

  // Values to check if the sampler failed or not
  bool succesful = true;
  std::string fail;
  int fail_iter;


  // Introduce additional index post_j that is used to calculate accurate storage positions in case of thinning
  int post_j = 1;
  // Begin Gibbs loop
  for (int j = 0; j < niter; j++) {
    // sample time varying beta.tilde parameters (NC parametrization)
    try {

      // Weave back into non-centered parameterization
      shrinkTVP::to_NCP(beta_nc_samp,
                        beta_samp,
                        theta_sr_samp,
                        beta_mean_samp);
      shrinkTVP::FFBS(beta_nc_samp,
                      y,
                      x,
                      theta_sr_samp,
                      beta_mean_samp,
                      sigma2_samp,
                      m_N_samp,
                      chol_C_N_inv_samp);
      // Weave into centered parameterization
      shrinkTVP::to_CP(beta_samp,
                       beta_nc_samp,
                       theta_sr_samp,
                       beta_mean_samp);
    } catch (...) {
      beta_nc_samp.fill(nanl(""));
      Rcout << theta_sr_samp << "\n" << sigma2_samp << "\n" << beta_mean_samp << "\n";
      if (succesful == true) {
        fail = "sample beta_nc";
        fail_iter = j + 1;
        succesful = false;
      }
    }



    // sample beta_mean and theta
    try {
      sample_theta_beta_mean(theta_samp,
                             beta_mean_samp,
                             beta_samp,
                             tau2,
                             s0,
                             S0);
      theta_sr_samp = arma::sqrt(theta_samp);
    } catch (...) {
      theta_samp.fill(nanl(""));
      beta_mean_samp.fill(nanl(""));
      if (succesful == true) {
        fail = "sample theta/beta";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    // sample sigma2 from homoscedastic or SV case
    try {
      if (sv) {
        arma::vec datastand = 2 * arma::log(arma::abs(y - x * beta_mean_samp - (x % beta_nc_samp.cols(1,N).t()) * theta_sr_samp));
        std::for_each(datastand.begin(), datastand.end(), shrinkTVP::res_protector);

        // update_sv needs sigma and not sigma^2
        double mu = sv_para(0);
        double phi = sv_para(1);
        double sigma = std::sqrt(sv_para(2));

        arma::vec cur_h = arma::log(sigma2_samp);
        stochvol::update_fast_sv(datastand, mu, phi, sigma, h0_samp, cur_h, r, prior_spec, expert);

        // Write back into sample object
        sigma2_samp = arma::exp(cur_h);

        // change back to sigma^2
        sv_para = {mu, phi, std::pow(sigma, 2)};

        std::for_each(sigma2_samp.begin(), sigma2_samp.end(), shrinkTVP::res_protector);

      } else {
        shrinkTVP::sample_sigma2(sigma2_samp,
                                 y,
                                 x,
                                 beta_nc_samp,
                                 beta_mean_samp,
                                 theta_sr_samp,
                                 c0,
                                 C0_samp);
      }
    } catch(...) {
      sigma2_samp.fill(nanl(""));
      if (succesful == true) {
        fail = "sample sigma2";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    if(sv == false) {
      try {
        C0_samp =  shrinkTVP::sample_C0(sigma2_samp,
                                        g0,
                                        c0,
                                        G0);
      } catch(...) {
        C0_samp = nanl("");
        if (succesful == true) {
          fail = "sample C0";
          fail_iter = j + 1;
          succesful = false;
        }
      }
    }

    // Increment index post_j if burn-in period is over
    if (j > nburn) {
      post_j++;
    }

    // Store everything (skip the steps when the sampler is used in the standalone case)
    if (!(niter == 1 && nburn == 0)){
      if ((post_j % nthin == 0) && (j >= nburn)) {
        // shrinkTVP::to_CP(beta_samp,
        //                  beta_nc_samp,
        //                  theta_sr_samp,
        //                  beta_mean_samp);
        sigma2_save.slice((post_j-1)/nthin) = sigma2_samp;
        theta_save.col((post_j-1)/nthin) = theta_samp;
        beta_mean_save.col((post_j-1)/nthin) = beta_mean_samp;
        beta_save.slice((post_j-1)/nthin) = beta_samp.t();
        m_N_save.slice((post_j-1)/nthin) = m_N_samp;
        chol_C_N_inv_save.slice((post_j-1)/nthin) = chol_C_N_inv_samp;

        if (sv == false) {
          C0_save((post_j-1)/nthin) = C0_samp;
        } else {
          sv_mu_save((post_j-1)/nthin) = sv_para(0);
          sv_phi_save((post_j-1)/nthin) = sv_para(1);
          sv_sigma2_save((post_j-1)/nthin) = sv_para(2);
        }
      }
    }

    // Increment progress bar
    if (arma::any(prog_rep_points == j)) {
      p.increment();
    }

    // Check for user interrupts
    if (j % 500 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Break loop if succesful is false
    if (!succesful) {
      break;
    }
  }

  // return everything as a nested list (due to size restrictions on Rcpp::Lists)

  return List::create(_["beta"] = beta_save,
                      _["beta_mean"] = beta_mean_save.t(),
                      _["theta"] = theta_save.t(),
                      _["sigma2"] = sigma2_save,
                      _["C0"] = C0_save,
                      _["sv_mu"] = sv_mu_save,
                      _["sv_phi"] = sv_phi_save,
                      _["sv_sigma2"] = sv_sigma2_save,
                      _["internals"] = List::create(
                        _["m_N"] = m_N_save,
                        _["chol_C_N_inv"] = chol_C_N_inv_save,
                        _["success_vals"] = List::create(
                          _["success"] = succesful,
                          _["fail"] = fail,
                          _["fail_iter"] = fail_iter),
                          _["x"] = x));


}


