#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;


  //  the estimation vector size should not be zero
  //  the estimation vector size should equal ground truth vector size
  assert(estimations.size() > 0);
  assert(estimations.size() == ground_truth.size());
  
  VectorXd c;
  // accumulate squared residuals
  for (int i=0; i < estimations.size(); ++i) {
    c = estimations[i] - ground_truth[i];
    rmse = rmse + VectorXd(c.array()*c.array());
  }

  // calculate the mean
  rmse = rmse/estimations.size();
  // calculate the squared root
  rmse = rmse.array().sqrt();
  
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);


  // check division by zero
  if ((px < 0.00001) && (py < 0.00001)){
      std::cout << "Divide By Zero" << std::endl;
      
  }
  else{
    float rho = std::sqrt(px*px + py*py);
      Hj << px/rho, py/rho, 0.0, 0.0,
            -py/(px*px+py*py), px/(px*px + py*py), 0.0, 0.0,
            py*(vx*py-vy*px)/std::pow(px*px + py*py, 1.5), px*(vy*px-vx*py)/std::pow(px*px+py*py,1.5),px/rho, py/rho; 
  }

  return Hj;
}