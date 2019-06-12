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
  int<lower = 1> nup;
  int<lower = 1> ndown;
  int<lower = 1> nparam;

  // Data
  real tempup[nup];
  real tempdown[ndown];

  real n2equilup[nup];
  real n2equildown[ndown];

  real n2up[nup];
  real n2down[ndown];

  // Fixed values
  real z; // Water depth (m)
  real tt; // Travel time between stations
  int<lower = 1> lag; // Time lag

  real light[ndown + lag];

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
  vector[nup] n2hat;

  // Global priors
  params[1] ~ normal(1, 5);
  params[2] ~ normal(Kmean, Ksd);
  sigma ~ normal(0,10);


  // -- MODEL 3
  // Two-station model without N consumption (DN base model)
  if(mod == 3){

    // Model loop
    for(i in 1:nup){
      n2hat[i] = n2up[i] + params[1] * tt / z + Kcor(tempup[i], params[2]) * tt * (n2equilup[i] - n2up[i] + n2equildown[i])/2;
      n2hat[i] /=  (1 + Kcor(tempup[i], params[2]) * tt / 2);
    }
  }

  // -- MODEL 4
  // Two station model with N consumption (DN + Nconsume)
  if(mod == 4){
    // Additional prior
    params[3] ~ normal(-1, 5);

    // Model loop
    for(i in 1:nup){
      n2hat[i] = n2up[i] + params[3] / z * sum(light[i:(i+lag)]) / sum(light) + params[1] * tt / z ;
      n2hat[i] +=  Kcor(tempup[i], params[2]) * tt * (n2equilup[i] - n2up[i] + n2equildown[i])/2;
      n2hat[i] /=  (1 + Kcor(tempup[i], params[2]) * tt / 2);
    }
  }

  // Likelihood
  n2down ~ normal( n2hat, sigma);

}
generated quantities{
    // Objects
  vector[nup] n2hat;
  vector[nup] n2pred;

// -- MODEL 3
  // Two-station model without N consumption (DN base model)
  if(mod == 3){

    // Model loop
    for(i in 1:nup){
      n2hat[i] = n2up[i] + params[1] * tt / z + Kcor(tempup[i], params[2]) * tt * (n2equilup[i] - n2up[i] + n2equildown[i])/2;
      n2hat[i] /=  (1 + Kcor(tempup[i], params[2]) * tt / 2);
    }
  }

  // -- MODEL 4
  // Two station model with N consumption (DN + Nconsume)
  if(mod == 4){
    // Model loop
    for(i in 1:nup){
      n2hat[i] = n2up[i] + params[3] / z * sum(light[i:(i+lag)]) / sum(light) + params[1] * tt / z ;
      n2hat[i] +=  Kcor(tempup[i], params[2]) * tt * (n2equilup[i] - n2up[i] + n2equildown[i])/2;
      n2hat[i] /=  (1 + Kcor(tempup[i], params[2]) * tt / 2);
    }
  }

 // Posterior predictive
 for(i in 1:nup){
  n2pred[i] = normal_rng( n2hat[i], sigma);
 }

}
