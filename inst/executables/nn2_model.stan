functions {

  // Function to estimate K at temperature 'temp'from K600.  Wanninkhof 1992.Temps 4-35 C
  real Kcor(real temp, real K600) {
    return K600/(600/(1615-(temp*92.15)+(2.349*temp^2)-(0.0240*temp^3)))^-0.5;
  }
}
data {
  // Model specifications
  int<lower = 1> mod;
  int<lower = 1> nobs;
  int<lower = 1> ncol;
  int<lower = 1> nparam;

  // Data
  real data_obj[nobs, ncol];

  // For models 3 and 4 data_obj is the following
  // -- Column 1 = tempup (celcius)
  // -- Column 2 = tempdown (celcius)
  // -- Column 3 = n2equilup (g m^-3)
  // -- Column 4 = n2equildown (g m^-3)
  // -- Column 5 = n2up (g m^-3)
  // -- Column 6 = n2down (g m^-3)

  // For models 1 and 2 data_obj is the following
  // -- Column 1 = temp (celcius)
  // -- Column 2 = n2equil (g m^-3)
  // -- Column 3 = n2 (g m^-3)
  // -- Column 4 = delta t (days^-1)


  // Fixed values
  real z; // Water depth (m)
  real tt; // Travel time between stations
  int<lower = 0> lag; // Time lag (if mod < 3 it is 0)

  real PPFD[nobs + lag]; // Photosynthetic photon flux density (PPFD/light) (mol m-2 s-1)
  real PPFDtotal; // Daily total of PPFD (mol m-2 s-1 d-1)

  // Prior specifications
  real Kmean;
  real Ksd;
}
parameters {
  // params[1] = DN = Denitrification rate (g m-2 d-1)
  // params[2] = K600 = Gas exchange rates (d-1)
  // params[3] = Nfix = N2 fixation rate (g m-2 d-1)
  real params[nparam];

  real<lower=0> sigma; // Variance
}
model {
  // Objects
  vector[nobs] n2hat;

  // Priors
  params[1] ~ normal(1, 5);
  params[2] ~ normal(Kmean, Ksd);
  sigma ~ normal(0,10);

  // Additional prior for model 2 and 4
  if(mod == 2){
    params[3] ~ normal(-1, 5);
  }
  if(mod == 4){
    params[3] ~ normal(-1, 5);
  }

  // -- One-station models (Models 1 and 1)
  if(mod < 3){
    // Initialize
    n2hat[1] = data_obj[1, 3];

    // Model loop
    for(i in 2:nobs){
      n2hat[i] = n2hat[i - 1] + params[1] * data_obj[i, 4] / z ;
      n2hat[i] +=  Kcor(data_obj[i, 1], params[2]) * data_obj[i, 4] * (data_obj[i, 2] + data_obj[i - 1, 2] - n2hat[i - 1])/2;

      // Add N consumption (DN + Nconsume) for model 2
      if(mod == 2){
         n2hat[i] += params[3] / z * PPFD[i] / PPFDtotal;
      }

      n2hat[i] /=  (1 + Kcor(data_obj[i, 1], params[2]) * data_obj[i, 4] / 2);
    }

    // Likelihood
    data_obj[2:,3] ~ normal( n2hat[2:], sigma);
  }


  // -- Two-station models (Models 3 and 4)
  if(mod > 2){
    // Model loop
    for(i in 1:nobs){
      n2hat[i] = data_obj[i, 5] + params[1] * tt / z ;
      n2hat[i] +=  Kcor(data_obj[i, 1], params[2]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2;

      // Add N consumption (DN + Nconsume) for model 4
      if(mod == 4){
         n2hat[i] += params[3] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal;
      }

      n2hat[i] /=  (1 + Kcor(data_obj[i, 1], params[2]) * tt / 2);
    }

    // Likelihood
    data_obj[,6] ~ normal( n2hat, sigma);
  }

}
generated quantities{
  // Objects
  vector[nobs] n2hat;
  vector[nobs] n2pred;

      // -- One-station models (Models 1 and 1)
  if(mod < 3){
    // Initialize
    n2hat[1] = data_obj[1, 3];

    // Model loop
    for(i in 2:nobs){
      n2hat[i] = n2hat[i - 1] + params[1] * data_obj[i, 4] / z ;
      n2hat[i] +=  Kcor(data_obj[i, 1], params[2]) * data_obj[i, 4] * (data_obj[i, 2] + data_obj[i - 1, 2] - n2hat[i - 1])/2;

      // Add N consumption (DN + Nconsume) for model 2
      if(mod == 2){
         n2hat[i] += params[3] / z * PPFD[i] / PPFDtotal;
      }

      n2hat[i] /=  (1 + Kcor(data_obj[i, 1], params[2]) * data_obj[i, 4] / 2);
    }

    // Posterior predictive
    n2pred[1] =  n2hat[1];
    for(i in 2:nobs){
      n2pred[i] = normal_rng( n2hat[i], sigma);
    }
  }


  // -- Two-station models (Models 3 and 4)
  if(mod > 2){
    // Model loop
    for(i in 1:nobs){
      n2hat[i] = data_obj[i, 5] + params[1] * tt / z ;
      n2hat[i] +=  Kcor(data_obj[i, 1], params[2]) * tt * (data_obj[i, 3] - data_obj[i, 5] + data_obj[i, 4])/2;

      // Add N consumption (DN + Nconsume) for model 4
      if(mod == 4){
         n2hat[i] += params[3] / z * sum(PPFD[i:(i+lag)]) / PPFDtotal;
      }

      n2hat[i] /=  (1 + Kcor(data_obj[i, 1], params[2]) * tt / 2);
    }

    // Posterior predictive
    for(i in 1:nobs){
      n2pred[i] = normal_rng( n2hat[i], sigma);
    }
  }
}
