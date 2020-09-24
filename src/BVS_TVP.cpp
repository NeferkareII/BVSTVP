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
  arma::mat mixprob(10, N);
  arma::vec mixprob_vec(mixprob.begin(), mixprob.n_elem, false);
  arma::ivec r(N);
  double h0_samp = as<double>(starting_vals["h0_st"]);
  double B011inv         = 1e-8;
  double B022inv         = 1e-12;
  bool Gammaprior        = true;
  double MHcontrol       = -1;
  int parameterization   = 3;
  bool centered_baseline = parameterization % 2; // 1 for C, 0 for NC baseline
  int MHsteps = 2;
  bool dontupdatemu = 0;
  double cT = N/2.0;
  double C0_sv = 1.5*Bsigma_sv;
  bool truncnormal = false;
  double priorlatent0 = -1;

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

      shrinkTVP::sample_beta_McCausland(beta_nc_samp,
                                        y,
                                        x,
                                        theta_sr_samp,
                                        sigma2_samp,
                                        beta_mean_samp,
                                        m_N_samp,
                                        chol_C_N_inv_samp);
    } catch (...) {
      beta_nc_samp.fill(nanl(""));
      if (succesful == true) {
        fail = "sample beta_nc";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    // Weave into centered parameterization
    shrinkTVP::to_CP(beta_samp,
                     beta_nc_samp,
                     theta_sr_samp,
                     beta_mean_samp);

    // sample beta_mean and theta
    try {
      sample_theta_beta_mean(theta_samp,
                             beta_mean_samp,
                             beta_samp,
                             tau2,
                             s0,
                             S0);

      // Weave back into NCP
      theta_sr_samp = arma::sqrt(theta_samp);
      // shrinkTVP::to_NCP(beta_nc_samp,
      //                   beta_samp,
      //                   theta_sr_samp,
      //                   beta_mean_samp);

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

        // shrinkTVP::to_CP(beta_samp,
        //                  beta_nc_samp,
        //                  theta_sr_samp,
        //                  beta_mean_samp);
        arma::vec datastand = 2 * arma::log(arma::abs(y - arma::sum(x % beta_samp.cols(1,N).t(), 1)));

        // update_sv needs sigma and not sigma^2
        sv_para(2) = std::sqrt(sv_para(2));

        arma::vec cur_h = arma::log(sigma2_samp);
        stochvol::update_sv(datastand,
                            sv_para,
                            cur_h,
                            h0_samp,
                            mixprob_vec,
                            r,
                            centered_baseline,
                            C0_sv,
                            cT,
                            Bsigma_sv,
                            a0_sv,
                            b0_sv,
                            bmu,
                            Bmu,
                            B011inv,
                            B022inv,
                            Gammaprior,
                            truncnormal,
                            MHcontrol,
                            MHsteps,
                            parameterization,
                            dontupdatemu,
                            priorlatent0);

        // Write back into sample object
        sigma2_samp = arma::exp(cur_h);

        std::for_each(sigma2_samp.begin(), sigma2_samp.end(), shrinkTVP::res_protector);


        // change back to sigma^2
        sv_para(2) = std::pow(sv_para(2), 2);

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
                          _["fail_iter"] = fail_iter)));


}


