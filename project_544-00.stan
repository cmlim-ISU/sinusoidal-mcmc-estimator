
data {
  int<lower=1> nobs;              // number of data points
  real y[nobs];                     // observation                           
  real<lower=0, upper=1> t[nobs]; //time step  
  //real<lower=-pi, upper=+pi> a[nobs]; //range for omega and gamma
}
parameters{
  real<lower=0> omega;    // phase angle
  real<lower=0> alpha;                   // amplitude   
  real<lower=-pi(), upper=pi()> gamma;                    // phase shift
  real<lower=0> sigma;                    // initial stv
  real<lower=0> epsilon;
}
transformed parameters{
  real yt[nobs];                         
  for (i in 1:nobs)
    yt[i] = alpha * sin(2*pi()*omega * t[i]+gamma);
                               
}
model {
  // priors and likelihood
  epsilon ~ normal(0, sigma);
  alpha ~ uniform(0,10); //the value should be 5 to 7 works dun touch
  sigma ~ normal(0, 1); // Prior for the standard deviation
  //omega ~ normal(5, 2); //works well should be 5~6
  omega ~ uniform(0, 10);
  
  //gamma ~ uniform(0, 6); #might be 0.5
  // Likelihood
  for (i in 1:nobs) {
    y[i] ~ normal(yt[i], epsilon);
  }
  
}

generated quantities {
  real y_rep[nobs];
  for (i in 1:nobs) {
    y_rep[i] = normal_rng(yt[i], epsilon);
  }
}


