#include <mathieu.h>
#include <iostream>

using namespace mathieu;

Mathieu::Mathieu(double a, double q, unsigned int n_max) :
  a(a), q(q), n_max(n_max), wronskian(-1),
  coeff_p(n_max+1),
  coeff_n(n_max+1)
{
  calculate_coeff();
  calculate_misc();
}


std::array<double, 4> Mathieu::eval(double x) const {
  std::array<double, 4> r;

  double arg = mu*x;
  double c = cos(arg) * coeff_p[0];
  double s = sin(arg) * coeff_n[0];
  double cp = -sin(arg) * mu * coeff_p[0];
  double sp = cos(arg) * mu * coeff_n[0];

  for (int i=1; i<n_max; i++) {
    double arg_p = mu - i*2.0;
    double arg_n = mu + i*2.0;

    double cos_p = cos(arg_p*x) * coeff_p[i];
    double cos_n = cos(arg_n*x) * coeff_n[i]; 

    double sin_p = sin(arg_p*x) * coeff_p[i];
    double sin_n = sin(arg_n*x) * coeff_n[i]; 

    c += (cos_p + cos_n);
    s += (sin_p + sin_n);

    cp -= (sin_p * arg_p + sin_n * arg_n);
    sp += (cos_p * arg_p + cos_n * arg_n);
  }

  r[0] = c;
  r[1] = s;
  r[2] = cp;
  r[3] = sp;

  return r;
}

void Mathieu::calculate_coeff() {
  std::vector<double> gamma(n_max);
  std::vector<double> d(n_max);
  std::vector<double> alpha(n_max);

  std::vector<double> ratio_n(n_max);
  std::vector<double> ratio_p(n_max);
  std::vector<double> ksi_mu_p(n_max);
  std::vector<double> ksi_mu_n(n_max);

  gamma[0] = -q/a;

  // Handle the singularity caused by a << 0

  for (int i=1; i<n_max; i++) {
    gamma[i] = q/(4*i*i-a);
    alpha[i] = gamma[i] * gamma[i-1];
  }

  d[0] = 1;
  d[1] = 1-2*alpha[1];
  d[2] = (2*alpha[1]+alpha[2]-1)*(alpha[2]-1);

  for (int i=3; i<n_max; i++) 
    d[i] = (1-alpha[i])*d[i-1] - alpha[i]*(1-alpha[i])*d[i-2] + alpha[i]*alpha[i-1]*alpha[i-1]*d[i-3];

    
  double (*f)(double) = a >= 0 ? (double (*)(double))&cos : (double (*)(double))&cosh;
  double aa = fabs(a);

  mu = acos( 1 - d[n_max-1] * ( 1 - f(M_PI* sqrt(aa)))) / M_PI;
  //std::cout << 1 - d[n_max-1] * ( 1 - f(M_PI* sqrt(aa))) << "\n";
  //std::cout << d[n_max-1] << "\n";
  //std::cout << "mu=  " << mu;

  if (mu != mu) mu = 0;

  for (int i=0; i<n_max; i++) {
    double arg_p = (-2*i - mu);
    double arg_n = ( 2*i - mu);
    ksi_mu_n[i] = (arg_p*arg_p-a)/q;
    ksi_mu_p[i] = (arg_n*arg_n-a)/q;
  }

  // Boundary condition, assume ratio becomes infinitely small
  ratio_n[n_max-1] = 0;
  ratio_p[n_max-1] = 0;

  for (int i=n_max-1;i>0;i--) {
    ratio_n[i-1] = -1/(ratio_n[i] + ksi_mu_n[i]);
    ratio_p[i-1] = -1/(ratio_p[i] + ksi_mu_p[i]);
  }

  coeff_p[0] = 1;
  coeff_n[0] = 1;

  double nf = 1;
  for (int i=1; i<n_max ; i++ ) {
    coeff_p[i] = coeff_p[i-1] * ratio_p[i-1];
    coeff_n[i] = coeff_n[i-1] * ratio_n[i-1];
    nf += (coeff_n[i]*coeff_n[i]+coeff_p[i]*coeff_p[i]);
  }
  nf = sqrt(nf);

  for (int i=0; i<n_max;i++) {
    coeff_p[i] /= nf;
    coeff_n[i] /= nf;
  }
}

void Mathieu::calculate_misc() {
  std::array<double, 4> r = eval(0);
  wronskian = r[0] * r[3];
}
