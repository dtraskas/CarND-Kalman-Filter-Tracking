#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  previous_timestamp_ = 0;
  is_initialized_ = false;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0, 0, 0.0225;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  noise_ax = 9;
  noise_ay = 9;

  VectorXd x_ = VectorXd(4);
  x_ << 0, 0, 0, 0;

  //state covariance matrix P
  MatrixXd P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

  MatrixXd F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  MatrixXd Q_ = MatrixXd(4, 4);
  Q_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        1, 0, 1, 0,
        0, 1, 0, 1;

  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);  
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    double px = 0.0;
    double py = 0.0;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);      
    
      px = rho * cos(phi);
      py = rho * sin(phi);      
      
      if (px == 0 && py == 0)
        return;
      
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      px = measurement_pack.raw_measurements_(0);
      py = measurement_pack.raw_measurements_(1);      

      if (px == 0 && py == 0)
        return;
          
    }    

    ekf_.x_ << px, py, 0, 0;

    // time starts from zero seconds ant not from Jan 1,1970
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  if (dt > 0.001){
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;

    //set the process covariance matrix Q
    ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
               0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
               dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
               0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

    ekf_.Predict();
  }
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}
