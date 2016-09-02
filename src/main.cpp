// This file contains three C++ functions for hindcasting CR data
// Note that emergence and recruitment are interchangeable

// double llRecruit(NumericVector params)
// double PrNotDetect(NumericVector params)
// NumericVector ExpectedCapture(NumericVector params)
// NumericVector ExpectedFirstCapture(NumericVector params)

#include <Rcpp.h>
using namespace Rcpp;

// Log-likelihood of the model that includes age-dependent capture rates
// [[Rcpp::export]]
double llRecruit(NumericVector params) {
  double tbar   = params[0];
  double sigma  = params[1];
  double alpha0 = params[2];
  double alpha1 = params[3];
  double beta0  = params[4];
  double beta1  = params[5];

  // access the global environment
  Environment env = Environment::global_env();

  Rcpp::IntegerMatrix y = env["y"]; // capture history for all individuals
  Rcpp::IntegerVector F = env["f"]; // index of first capture event
  Rcpp::IntegerVector L = env["l"]; // index of last caoture event
  Rcpp::NumericVector E = env["E"]; // effort across capture events
  Rcpp::IntegerVector T = env["T"]; // day of capture events
  int TF = env["T.F"]; // earliest day of emergence
  int TL = env["T.L"]; // latest day of emergence

  int I = F.size(); // number of individuals caught
  int J = E.size(); // number of sampling days

  // copy the global vectors f and l (don't want to change them)
  int f[I];
  int l[I];
  for (int i = 0; i < I; ++i) {
    f[i] = F(i) - 1; // decrement index so 0 is the base index
    l[i] = L(i) - 1; // decrement index so 0 is the base index
  }

  // make sure TL is not less than last day sampled
  if (TL < T(J-1)) {
    TL = T(J-1);
  }

  int EmergeDays = TL - TF + 1; // possible number of emergence days

  // calculate next sampling time
  int n[EmergeDays]; // next sampling event after emergence
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    n[d2] = 0.0;
  }
  for (int d = TF; d <= T(0); ++d) {
    n[d-TF] = 0; // first sampling time
  }
  for (int j = 1; j < J; ++j) {
    for (int d = T(j-1) + 1; d <= T(j); ++d) {
      n[d-TF] = j;
    }
  }

  // create the vector of emergence probabilities
  double u[EmergeDays];
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    u[d2] = 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = exp(-0.5*pow((i + TF - tbar)/sigma,2));
    sum += u[i];
  }
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = u[i]/sum; // normalise so u is a probability
  }

  // calculate probabilities of catching butterflies
  double c[J][EmergeDays]; // prob. catch on sampling day j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      c[j][d2] = 0.0;
    }
  }
  int d;
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling day
      if (T[j] >= d) { // a valid sampling day
        c[j][d2] = 1 - exp(-alpha0*exp(alpha1*(T[j] - d)) * E(j));
      } else {
        c[j][d2] = 0.0;
      }
    }
  }

  // calculate survival to next sampling time
  double s[J][EmergeDays]; // prob. survive to sample j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      s[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling event
      if (T[j] >= d) { // a valid sampling day
        if (beta1 != 0.0) {
          s[j][d2] = exp(-(beta0 / beta1) * (exp(beta1 * (T[j] - d)) - 1.0));
        } else {
          s[j][d2] = exp(-beta0*(T[j] - d));
        }
      } else {
        s[j][d2] = 0.0;
      }
    }
  }

  // calculate survivals between sampling times (from j-1 to j)
  double z[J][EmergeDays]; // prob survive between sampling times
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      z[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = n[d2] + 1; j < J; ++j) {
      z[j][d2] = s[j][d2]/s[j-1][d2];
    }
  }

  double qprime = 0.0; // probability an individual will not be seen
  double prod1 = 1.0;
  double prod2 = 1.0;
  double sum1 = 0.0;
  double sum2 = 0.0;

  // calculate the probability an individual is not seen at all = q
  double q[EmergeDays]; // prob not seen at all if emerged on day
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    q[d2] = 0.0;
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    prod1 = 1.0;
    for (int j = n[d2] + 1; j < J; ++j) {
      prod1 = prod1 * z[j][d2] * (1 - c[j][d2]);
    }
    sum1 = prod1;
    for (int j = n[d2] + 1; j < J; ++j) {
      prod1 = 1.0;
      for (int k = n[d2] + 1; k <= j-1; ++k) {
        prod1 = prod1 * z[k][d2] * (1 - c[k][d2]);
      }
      sum1 += (1-z[j][d2])*prod1;
    }
    q[d2] = 1.0 - s[n[d2]][d2] + s[n[d2]][d2] * (1.0 - c[n[d2]][d2]) * sum1;
    qprime += u[d2] * q[d2];
  }

  // calculate the probabilities of first detection
  double v[J][EmergeDays]; // prob. of first detection
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      v[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // d2 is day of emergence
    for (int f2 = n[d2]; f2 < J; ++f2) { // f2 is the first observation
      prod1 = s[n[d2]][d2]*c[f2][d2];
      for (int j = n[d2] + 1; j <= f2; ++j) {
        prod1 *= (1.0 - c[j-1][d2]) * z[j][d2];
      }
      v[f2][d2] = prod1;
    }
  }

  // calculate probability of capture history
  double w[I][EmergeDays]; // store for each individual
  for (int i = 0; i < I; ++i) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      w[i][d2] = 0.0;
    }
  }

  // calculate log-likelihood across all indivduals
  for (int i = 0; i < I; ++i) {
    for (int d2 = 0; d2 <= T(f[i]) - TF; ++d2) { // valid emergence days
      prod2 = 1.0;
      for (int j = f[i] + 1; j <= l[i]; ++j) { // survived periods
        prod2 *= z[j][d2]*(y(i,j)*c[j][d2] + (1.0-y(i,j))*(1.0-c[j][d2]));
      }

      prod1 = 1.0;
      for (int j = l[i]+1; j < J; ++j) {
        prod1 *= z[j][d2]*(1.0 - c[j][d2]);
      }
      sum1 = prod1;

      for (int j = l[i]+1; j < J; ++j) {
        prod1 = 1.0 - z[j][d2];
        for(int k = l[i] + 1; k <= j-1; ++k) {
          prod1 *= z[k][d2] * (1.0 - c[k][d2]);
        }
        sum1 += prod1;
      }

      w[i][d2] = prod2 * sum1;
    }
  }

  sum1 = -I*log(1.0 - qprime);
  int f2; // sample when individual is first observed
  for (int i = 0; i < I; ++i) {
    f2 = f[i];
    sum2 = 0.0;
    for (int d2 = 0; d2 <= T(f2) - TF; ++d2) { // valid emergence days
      sum2 += u[d2]*v[f2][d2]*w[i][d2];
    }
    sum1 += log(sum2);
  }

  return(sum1); // log-likelihood of probability of all data given model
}

// Probability an individual is not detected during sampling
// [[Rcpp::export]]
double PrNotDetect(NumericVector params) {
  double tbar   = params[0];
  double sigma  = params[1];
  double alpha0 = params[2];
  double alpha1 = params[3];
  double beta0  = params[4];
  double beta1  = params[5];

  // access the global environment
  Environment env = Environment::global_env();

  Rcpp::NumericVector E = env["E"]; // effort across capture events
  Rcpp::IntegerVector T = env["T"]; // day of capture events
  int TF = env["T.F"]; // earliest day of emergence
  int TL = env["T.L"]; // latest day of emergence

  int J = E.size(); // number of sampling days

  int EmergeDays = TL - TF + 1; // possible number of emergence days

  // calculate next sampling time
  int n[EmergeDays]; // next sampling event after emergence
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    n[d2] = 0.0;
  }
  for (int d = TF; d <= T(0); ++d) {
    n[d-TF] = 0; // first sampling time
  }
  for (int j = 1; j < J; ++j) {
    for (int d = T(j-1) + 1; d <= T(j); ++d) {
      n[d-TF] = j;
    }
  }

  // create the vector of emergence probs
  double u[EmergeDays];
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    u[d2] = 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = exp(-0.5*pow((i + TF - tbar)/sigma,2));
    sum += u[i];
  }
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = u[i]/sum; // normalise so u is a probability
  }

  // calculate probabilities of catching butterflies
  double c[J][EmergeDays]; // prob. catch on sampling day j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      c[j][d2] = 0.0;
    }
  }
  int d;
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling day
      if (T[j] >= d) { // a valid sampling day
        c[j][d2] = 1 - exp(-alpha0*exp(alpha1*(T[j] - d)) * E(j));
      } else {
        c[j][d2] = 0.0;
      }
    }
  }

  // calculate survival to next sampling time
  double s[J][EmergeDays]; // prob. survive to sampling event j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      s[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling event
      if (T[j] >= d) { // a valid sampling day
        if (beta1 != 0.0) {
          s[j][d2] = exp(-(beta0 / beta1) * (exp(beta1 * (T[j] - d)) - 1.0));
        } else {
          s[j][d2] = exp(-beta0*(T[j] - d));
        }
      } else {
        s[j][d2] = 0.0;
      }
    }
  }

  // calculate survivals between sampling times (from j-1 to j)
  double z[J][EmergeDays]; // prob survive between sampling times
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      z[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = n[d2] + 1; j < J; ++j) {
      z[j][d2] = s[j][d2]/s[j-1][d2];
    }
  }

  double qprime = 0.0; // probability an individual will not be seen
  double prod1 = 1.0;
  double prod2 = 1.0;
  double sum1 = 0.0;
  double sum2 = 0.0;

  // calculate the probability an individual is not seen at all = q
  double q[EmergeDays]; // prob not seen at all if emerged on day
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    q[d2] = 0.0;
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    prod1 = 1.0;
    for (int j = n[d2] + 1; j < J; ++j) {
      prod1 = prod1 * z[j][d2] * (1 - c[j][d2]);
    }
    sum1 = prod1;
    for (int j = n[d2] + 1; j < J; ++j) {
      prod1 = 1.0;
      for (int k = n[d2] + 1; k <= j-1; ++k) {
        prod1 = prod1 * z[k][d2] * (1 - c[k][d2]);
      }
      sum1 += (1-z[j][d2])*prod1;
    }
    q[d2] = 1.0 - s[n[d2]][d2] + s[n[d2]][d2] * (1.0 - c[n[d2]][d2]) * sum1;
    qprime += u[d2] * q[d2];
  }

  return(qprime); // log-likelihood of probability of all data given model
}

// Return a vector of expected capture fractions per capture event
// [[Rcpp::export]]
NumericVector ExpectedCapture(NumericVector params) {
  double tbar   = params[0];
  double sigma  = params[1];
  double alpha0 = params[2];
  double alpha1 = params[3];
  double beta0  = params[4];
  double beta1  = params[5];

  // access the global environment
  Environment env = Environment::global_env();

  Rcpp::NumericVector E = env["E"]; // effort across capture events
  Rcpp::IntegerVector T = env["T"]; // day of capture events
  int TF = env["T.F"]; // earliest day of emergence
  int TL = env["T.L"]; // latest day of emergence

  int J = E.size(); // number of sampling days

  int EmergeDays = TL - TF + 1; // possible number of emergence days

  // create the vector of emergence probs
  double u[EmergeDays];
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    u[d2] = 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = exp(-0.5*pow((i + TF - tbar)/sigma,2));
    sum += u[i];
  }
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = u[i]/sum; // normalise so u is a probability
  }

  // calculate probabilities of catching butterflies
  double c[J][EmergeDays]; // prob. catch on sampling day j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      c[j][d2] = 0.0;
    }
  }
  int d;
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling day
      if (T[j] >= d) { // a valid sampling day
        c[j][d2] = 1 - exp(-alpha0*exp(alpha1*(T[j] - d)) * E(j));
      } else {
        c[j][d2] = 0.0;
      }
    }
  }

  // calculate survival to next sampling time
  double s[J][EmergeDays]; // prob. survive to sampling event j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      s[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling event
      if (T[j] >= d) { // a valid sampling day
        if (beta1 != 0.0) {
          s[j][d2] = exp(-(beta0 / beta1) * (exp(beta1 * (T[j] - d)) - 1.0));
        } else {
          s[j][d2] = exp(-beta0*(T[j] - d));
        }
      } else {
        s[j][d2] = 0.0;
      }
    }
  }

  NumericVector g(J); // prob not seen at all if emerged on day
  double sum1 = 0.0;

  for (int j = 0; j < J; ++j) {
    sum1 = 0.0;
    for (int d2 = 0; d2 <= T(j) - TF; ++d2) {
      sum1 += u[d2]*s[j][d2]*c[j][d2];
    }
    g[j] = sum1;
  }

  return(g); // probability vector
}

// Return a vector of expected fractions of first detections per capture event
// [[Rcpp::export]]
NumericVector ExpectedFirstCapture(NumericVector params) {
  double tbar   = params[0];
  double sigma  = params[1];
  double alpha0 = params[2];
  double alpha1 = params[3];
  double beta0  = params[4];
  double beta1  = params[5];

  // access the global environment
  Environment env = Environment::global_env();

  Rcpp::NumericVector E = env["E"]; // effort across capture events
  Rcpp::IntegerVector T = env["T"]; // day of capture events
  int TF = env["T.F"]; // earliest day of emergence
  int TL = env["T.L"]; // latest day of emergence

  int J = E.size(); // number of sampling days

  // make sure TL is not less than last day sampled
  if (TL < T(J-1)) {
    TL = T(J-1);
  }

  int EmergeDays = TL - TF + 1; // possible number of emergence days

  // calculate next sampling time
  int n[EmergeDays]; // next sampling event after emergence
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    n[d2] = 0.0;
  }
  for (int d = TF; d <= T(0); ++d) {
    n[d-TF] = 0; // first sampling time
  }
  for (int j = 1; j < J; ++j) {
    for (int d = T(j-1) + 1; d <= T(j); ++d) {
      n[d-TF] = j;
    }
  }

  // create the vector of emergence probs
  double u[EmergeDays];
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    u[d2] = 0.0;
  }
  double sum = 0.0;
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = exp(-0.5*pow((i + TF - tbar)/sigma,2));
    sum += u[i];
  }
  for (int i = 0; i < EmergeDays; ++i) {
    u[i] = u[i]/sum; // normalise so u is a probability
  }

  // calculate probabilities of catching butterflies
  double c[J][EmergeDays]; // prob. catch on sampling day j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      c[j][d2] = 0.0;
    }
  }
  int d;
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling day
      if (T[j] >= d) { // a valid sampling day
        c[j][d2] = 1 - exp(-alpha0*exp(alpha1*(T[j] - d)) * E(j));
      } else {
        c[j][d2] = 0.0;
      }
    }
  }

  // calculate survival to next sampling time
  double s[J][EmergeDays]; // prob. survive to sampling event j if emerged on day d
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      s[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = 0; j < J; ++j) { // sampling event
      if (T[j] >= d) { // a valid sampling day
        if (beta1 != 0.0) {
          s[j][d2] = exp(-(beta0 / beta1) * (exp(beta1 * (T[j] - d)) - 1.0));
        } else {
          s[j][d2] = exp(-beta0*(T[j] - d));
        }
      } else {
        s[j][d2] = 0.0;
      }
    }
  }

  // calculate survivals between sampling times (from j-1 to j)
  double z[J][EmergeDays]; // prob survive between sampling times
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      z[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    d = TF + d2; // d = day of emergence
    for (int j = n[d2] + 1; j < J; ++j) {
      z[j][d2] = s[j][d2]/s[j-1][d2];
    }
  }

  double qprime = 0.0; // probability an individual will not be seen
  double prod1 = 1.0;
  double prod2 = 1.0;
  double sum1 = 0.0;
  double sum2 = 0.0;

  // calculate the probability an individual is not seen at all = q
  double q[EmergeDays]; // prob not seen at all if emerged on day
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
    q[d2] = 0.0;
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) {
    prod1 = 1.0;
    for (int j = n[d2] + 1; j < J; ++j) {
      prod1 = prod1 * z[j][d2] * (1 - c[j][d2]);
    }
    sum1 = prod1;
    for (int j = n[d2] + 1; j < J; ++j) {
      prod1 = 1.0;
      for (int k = n[d2] + 1; k <= j-1; ++k) {
        prod1 = prod1 * z[k][d2] * (1 - c[k][d2]);
      }
      sum1 += (1-z[j][d2])*prod1;
    }
    q[d2] = 1.0 - s[n[d2]][d2] + s[n[d2]][d2] * (1.0 - c[n[d2]][d2]) * sum1;
    qprime += u[d2] * q[d2];
  }

  // calculate the probabilities of first detection
  double v[J][EmergeDays]; // prob. of first detection
  for (int j = 0; j < J; ++j) {
    for (int d2 = 0; d2 < EmergeDays; ++d2) { // valid emergence days
      v[j][d2] = 0.0;
    }
  }
  for (int d2 = 0; d2 < EmergeDays; ++d2) { // d2 is day of emergence
    for (int f2 = n[d2]; f2 < J; ++f2) { // f2 is the first observation
      prod1 = s[n[d2]][d2]*c[f2][d2];
      for (int j = n[d2] + 1; j <= f2; ++j) {
        prod1 *= (1.0 - c[j-1][d2]) * z[j][d2];
      }
      v[f2][d2] = prod1;
    }
  }

  NumericVector g(J); // prob seen for the first time if emerged on day
  sum1 = 0.0;

  for (int j = 0; j < J; ++j) {
    sum1 = 0.0;
    for (int d2 = 0; d2 <= T(j) - TF; ++d2) {
      sum1 += u[d2]*v[j][d2];
    }
    g[j] = sum1;
  }

  return(g); // probability vector
}
