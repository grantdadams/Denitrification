functions {

// Function to estimate K at temperature 'temp'from K600.  Wanninkhof 1992.Temps 4-35 C
  real Kcor(real temp, real k600) {
    return K600/(600/(1615-(temp*92.15)+(2.349*temp^2)-(0.0240*temp^3)))^-0.5;
  }
}
data {
// Model specifications
  int<lower = 1> model;
  int<lower = 1> nobs;
    int<lower = 1> nup;
      int<lower = 1> ndown;

  // Data
  real tempup[n];
  real tempdown[n];

  real n2equilup[n];
  real n2equildown[n];

  real light[n];

  // Fixed values
  real z; // Water depth (m)
  real tt; // Travel time between stations

  // Prior specifications
  real Kmean;
  real Ksd;
}
parameters {
  real DN; // Denitrification rate (g m-2 d-1), 
  real K; // Gas exchange rates (d-1) 
  real<lower=0> sigma; // Variance
}
model {
  // Priors
  DN ~ normal(1, 5);  
  K ~ normal(Kmean, Ksd);
  sigma ~ normal(0,10);


// Two station model without N consumption
  if(model == 3){
for(t in 1:n){
n2hat[i] = (n2up[i] + DN * tt / z + ( Kcor(tempup[i], K)) * tt * (n2equilup[i] - n2up[i] + n2equildown[i])/2)
n2hat[i] /=  (1 + Kcor(tempup[i], K) * tt / 2)
}
  }
  
  // Likelihood
  for(i in 1:n)
  weight[i] ~ normal( a * length[i]^b, sigma);
}
