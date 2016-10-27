#include <stdio.h>
#include <iostream>
#include <cmath>
using namespace std;

#define PARAMS_FROM_KEYBOARD
bool DEBUG = false;

double* u;
double* u_old;

double* u2;
double* u2_old;
double* kappa;
double* omega;

double L;
int N;
double h;
double tau;
int J;
int step;
double t = 0;

void print_layer();

void init_parameters()
{
#ifdef PARAMS_FROM_KEYBOARD
  L = 1;
  h = 0.1;
  J = 2;
  cout << "L: " << L << endl;
  cout << "h: " << h << endl;
  cout << "J: " << J << endl;
#else
  cout << "L: "; cin >> L;
  cout << "h: "; cin >> h;
  cout << "J: "; cin >> J;
#endif
  N = (int)(L / h);
  cout << "N: " << N << endl;
}

double a(double x, double t)
{
  return 1;
}

double psi1(double x)
{
  return x;
}

double psi2(double t)
{
  return t*131;
}

double psi3(double t)
{
  return t*81 + L;
}

double phi(double x, double t)
{
  return x + t;
}

double alpha0(double t)
{
  return 0;
}

double betta0(double t)
{
  return 1;
}

double alpha1(double t)
{
  return 0;
}

double betta1(double t)
{
  return 1;
}

bool check_conditions()
{
  return (psi1(0) == psi2(0)) && (psi1(L) == psi3(0));
}

void init_containers()
{
  u = new double[N+1];
  u_old = new double[N+1];
  u2 = new double[N+1];
  u2_old = new double[N+1];
  omega = new double[N+1];
  kappa = new double[N+1];
}

void fill_first_layer()
{
  for(int i = 0; i <= N; i++) {
    u[i] = psi1(i*h);
    u_old[i] = psi1(i*h);
    u2[i] = psi1(i*h);
    u2_old[i] = psi1(i*h);
  }
}

void fill_bolders()
{
  u[0] = 1.0/(betta0(t) - alpha0(t)/h)*psi2(t) - alpha0(t)/(betta0(t)*h - alpha0(t))*u_old[1];
  u[N] = 1.0/(alpha1(t)/h + betta1(t))*psi3(t) + alpha1(t)/(alpha1(t) + betta1(t)*h)*u_old[N-1];
  u2[0] = 1.0/(betta0(t) - alpha0(t)/h)*psi2(t) - alpha0(t)/(betta0(t)*h - alpha0(t))*u_old[1];
  u2[N] = 1.0/(alpha1(t)/h + betta1(t))*psi3(t) + alpha1(t)/(alpha1(t) + betta1(t)*h)*u_old[N-1];
}

void calculate_tau() {
  double u_max = 0;

  for(int i = 0; i <= N; i++) {
    double val = a(i*h, t)*a(i*h, t)*(1+u[i]*u[i]);
    if (val > u_max) {
      u_max = val;
    }
  }

  tau = h*h/(2*u_max);
}

void swap(double** first, double** second)
{
  double** tmp = first;
  *first = *second;
  *second = *tmp;
}

void calculate_explicit()
{
  swap(&u, &u_old);
  for(int i = 1; i < N; i++) {
    double ath = (tau/(h*h))*a(i*h, t)*a(i*h, t)*(1+u_old[i]*u_old[i]);

    u[i] = (1-2*ath)*u_old[i] +
           ath *     u_old[i+1] +
           ath *     u_old[i-1] +
           tau *     phi(i*h, t);
  }
}

bool check_resolution(double a, double c, double b)
{
  if (abs(c) >= abs(a) + abs(b)) {
    return true;
  } else {
    printf("WARNING: |%lf| < |%lf| + |%lf|\n", c, a, b);
    return false;
  }
}

void calculate_implicit()
{
  swap(&u2, &u2_old);
  double ath = (tau/(h*h))*a(0, t)*a(0, t)*(1+u2_old[0]*u2_old[0]);
  double A;
  double C = 1-2*ath;
  double B = -ath;
  double F = tau*phi(h, t+tau) + u2_old[1]; // phi on j+1

  kappa[1] = -B/C;
  omega[1] = F/C;
  printf("0 %lf %lf %lf\n", C, B, F);
  check_resolution(0, C, B);

  for(int i = 1; i < N; i++) {
    ath = (tau/(h*h))*a(i*h, t)*a(i*h, t)*(1+u2_old[i]*u2_old[i]);
    A = -ath;
    C = 1+2*ath;
    B = -ath;
    F = tau*phi(i*h, t+tau) + u2_old[i]; // phi on j+1

    printf("%lf %lf %lf %lf\n", A, C, B, F);
    check_resolution(A, C, B);

    kappa[i+1] = -B/(A*kappa[i] + C);
    omega[i+1] = (F - A*omega[i])/(A*kappa[i] + C);
  }

  ath = (tau/(h*h))*a(N*h, t)*a(N*h, t)*(1+u2_old[N]*u2_old[N]);
  A = -ath;
  C = 1-2*ath;
  F = tau*phi(N*h, t+tau) + u2_old[N]; // phi on j+1
  printf("%lf %lf %lf %lf\n", A, C, B, F);
  check_resolution(A, C, B);

  u2[N] = (F - A*omega[N])/(C + A*kappa[N]);
  printf("%lf\n", u2[N]);
  for(int i = N-1; i >= 0; i--) {
  	u2[i] = kappa[i+1]*u2[i+1] + omega[i+1];
  }
}

void calculate()
{
  if (DEBUG) {
    printf("t: %0.5lf\nTau: %0.5lf\n", t, tau);
    print_layer();
  }
  for(step = 0; step <= J; step++) {
    calculate_tau();
    if (DEBUG) {
      printf("t: %0.5lf\nTau: %0.5lf\n", t, tau);
    }
    calculate_explicit();
    calculate_implicit();
    t += tau;
    if (DEBUG) {
      print_layer();
    }
    fill_bolders();
  }
}

void print_results()
{
  cout << "x explicit" << endl;
  for(int i = 0; i <= N; i++) {
    printf("%0.2lf %0.2lf %0.2lf %0.2lf\n", i*h, u[i], u2[i], abs(u2[i]-u[i]));
  }
}

void print_layer()
{
  for(int i = 0; i <= N; i++) {
    printf(" %0.2lf", u[i]);
  }

  printf("\n");
  for(int i = 0; i <= N; i++) {
    printf(" %0.2lf", u2[i]);
  }

  printf("\n");
}

int main()
{
  init_parameters();
  if(!check_conditions) {
    printf("Wrong input functions");
    return 0;
  }

  init_containers();
  fill_first_layer();
  calculate();
  print_results();

  return 0;
}
