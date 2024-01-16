#include "Constraints.h"
#include "DecayTable.h"
#include "HBHS.h"
#include "SM.h"
#include "THDM.h"
#include <iostream>
#include <math.h>
#include <getopt.h>
#include <cstdlib>
#include <vector>

using namespace std;

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;
    double delta = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        double value = start + i * delta;
        result.push_back(value);
    }
    return result;
}

int main(int argc, char *argv[]) {

  double mh = 125.;
  double mH = 200.;
  double mA = 150.;
  double mC = mH;
  double sba = 1.0;
  double lambda_6 = 0.;
  double lambda_7 = 0.;
  double tb = 2.0;
  int yt = 1;

  int opt;
  static struct option long_options[] = {
      {"mh", required_argument, 0, 'h'},
      {"mH", required_argument, 0, 'H'},
      {"mA", required_argument, 0, 'A'},
      {"mC", required_argument, 0, 'C'},
      {"sba", required_argument, 0, 's'},
      {"lambda_6", required_argument, 0, '6'},
      {"lambda_7", required_argument, 0, '7'},
      {"tb", required_argument, 0, 't'},
      {"yt", required_argument, 0, 'y'},
      {0, 0, 0, 0}
  };
  int option_index = 0;
  while ((opt = getopt_long(argc, argv, "h:H:A:C:s:6:7:t:y:m:", long_options, &option_index)) != -1) {
      switch (opt) {
          case 'h':
              mh = strtod(optarg, nullptr);
              break;
          case 'H':
              mH = strtod(optarg, nullptr);
              break;
      	  case 'A':
              mA = strtod(optarg, nullptr);
              break;
          case 'C':
              mC = strtod(optarg, nullptr);
              break;
          case 's':
              sba = strtod(optarg, nullptr);
              break;
      	  case '6':
              lambda_6 = strtod(optarg, nullptr);
              break;
          case '7':
              lambda_7 = strtod(optarg, nullptr);
              break;
       	  case 't':
              tb = strtod(optarg, nullptr);
              break;
          case 'y':
              yt = strtod(optarg, nullptr);
              break;
          case '?':
              std::cerr << "Unknown option: " << optopt << std::endl;
              return 1;
          default:
              std::cerr << "Unexpected option: " << opt << std::endl;
              return 1;
      }
  }


  double m12_2 = mH*mH*cos(atan(tb))*sin(atan(tb));
  double interval = 100.0;
  bool positive = true;

  // Reference SM Higgs mass for EW precision observables
  double mh_ref = 125.;

  // Create SM and set parameters
  SM sm;
  sm.set_qmass_pole(6, 172.5);
  sm.set_qmass_pole(5, 4.75);
  sm.set_qmass_pole(4, 1.42);
  sm.set_lmass_pole(3, 1.77684);
  sm.set_alpha(1. / 127.934);
  sm.set_alpha0(1. / 137.0359997);
  sm.set_alpha_s(0.119);
  sm.set_MZ(91.15349);
  sm.set_MW(80.36951);
  sm.set_gamma_Z(2.49581);
  sm.set_gamma_W(2.08856);
  sm.set_GF(1.16637E-5);

  // Create 2HDM and set SM parameters
  THDM model;
  model.set_SM(sm);

  bool valid = false;
  for (int i = 0; i < 100; i++) {

    bool pset =
        model.set_param_phys(mh, mH, mA, mC, sba, lambda_6, lambda_7, m12_2, tb);

    if (!pset) {
      cerr << "The specified parameters are not valid\n";
      return -1;
    }

    // Set Yukawa couplings to type II
    model.set_yukawas_type(yt);

    // Print the parameters in different parametrizations
    model.print_param_phys();
    //model.print_param_gen();
    //model.print_param_higgs();
    //model.print_param_hybrid();

    // Prepare to calculate observables
    Constraints constr(model);

    double S, T, U, V, W, X;

    constr.oblique_param(mh_ref, S, T, U, V, W, X);

    printf("\nConstraints:\n");
    printf("  Potential stability: %s\n",
           (constr.check_stability() ? "OK" : "Not OK"));
    printf(" Tree-level unitarity: %s\n",
           (constr.check_unitarity() ? "OK" : "Not OK"));
    printf("       Perturbativity: %s\n",
           (constr.check_perturbativity() ? "OK" : "Not OK"));

    valid = (constr.check_perturbativity() && constr.check_unitarity() && constr.check_stability());
    if (valid) {
      break;
    }

    if (constr.check_stability()) {
      if (!positive) {
        interval = interval/2.0;
      }
      m12_2 = m12_2 + interval;
      positive = true;
    } else {
      if (positive) {
        interval = interval/2.0;
      }
      m12_2 = m12_2 - interval;
      positive = false;
    }

    std::cout << m12_2 << " " << interval << std::endl;
  }

  if (valid) {
    std::cout << "Valid" << std::endl;
  } else {
    std::cout << "Not valid" << std::endl;
  }

}
