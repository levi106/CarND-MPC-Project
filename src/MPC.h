#ifndef MPC_H
#define MPC_H

#include <vector>
#include <fstream>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

  static const size_t N;
  static const double dt;
  static const double Lf;
private:
  std::ofstream f_;
};

#endif /* MPC_H */
