#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;


/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  ekf_.x_ = VectorXd(4);
  ekf_.P_ = MatrixXd(4,4);
  ekf_.F_ = MatrixXd(4,4);
  ekf_.Q_ = MatrixXd(4,4);
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    float x, y;
    
    // first measurement (Don't know what this is for)
    cout << "EKF: " << endl;
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
	  y = sin(measurement_pack.raw_measurements_(1))*measurement_pack.raw_measurements_(0);
      x = cos(measurement_pack.raw_measurements_(1))*measurement_pack.raw_measurements_(0);
      
      // Propogation of uncertainty from Radar to x and y based on equations from: https://en.wikipedia.org/wiki/Propagation_of_uncertainty
      float sin_cov = cos(measurement_pack.raw_measurements_(1))*cos(measurement_pack.raw_measurements_(1))*0.0009;
      float y_cov = sin(measurement_pack.raw_measurements_(1))*sin(measurement_pack.raw_measurements_(1))*0.09 + 				measurement_pack.raw_measurements_(0)*measurement_pack.raw_measurements_(0)*sin_cov + sin_cov*0.09;
      float cos_cov = sin(measurement_pack.raw_measurements_(1))*sin(measurement_pack.raw_measurements_(1))*0.0009;
      float x_cov = cos(measurement_pack.raw_measurements_(1))*cos(measurement_pack.raw_measurements_(1))*0.09 + 				measurement_pack.raw_measurements_(0)*measurement_pack.raw_measurements_(0)*cos_cov + cos_cov*0.09;
                                                                  
                                                                  
      ekf_.P_ << x_cov, 0, 0, 0,
            0, y_cov, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1.;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      x = measurement_pack.raw_measurements_(0);
      y = measurement_pack.raw_measurements_(1);
      ekf_.P_ << .0225, 0, 0, 0,
            0, .0225, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
    }
	ekf_.x_ << x, y, 0.0, 0.0;
    
    
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
 
  

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  
  
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // convert to seconds or already in seconds?
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_  << 1.0, 0.0, dt, 0.0,
            0.0, 1.0, 0.0, dt,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0; 
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  ekf_.Q_ << std::pow(dt,4.0)*noise_ax/4.0, 0.0,  std::pow(dt,3.0)*noise_ax/2.0, 0.0,
            0.0, std::pow(dt,4.0)*noise_ay/4.0, 0.0, std::pow(dt,3.0)*noise_ay/2.0,
            std::pow(dt,3.0)*noise_ax/2.0, 0.0, std::pow(dt,2.0)*noise_ax, 0.0,
            0.0, std::pow(dt,3.0)*noise_ay/2.0, 0.0, std::pow(dt,2.0)*noise_ay;
  
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
    //Check for dividing by zero
    if (fabs(ekf_.x_(0)*ekf_.x_(0) + ekf_.x_(1)*ekf_.x_(1)) < .01){
      return; // skip thisw radar measurement
    }
    // TODO: Radar updates
    ekf_.H_ =  tools.CalculateJacobian(ekf_.x_);
    
    ekf_.R_ = R_radar_;
    
    VectorXd rad_meas = VectorXd(3);
    rad_meas << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
    
	ekf_.UpdateEKF(rad_meas);
  } 
  else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    VectorXd laser_meas = VectorXd(2);
    laser_meas << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
    ekf_.Update(laser_meas);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
