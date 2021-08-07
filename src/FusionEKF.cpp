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
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0;

  Hj_ << 1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1, 0;

  ekf_.x_ = VectorXd(4);
  ekf_.Q_ = MatrixXd(4, 4);
  
  noise_ax = 9; // sax2
  noise_ay = 9; // say2
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  
  cout << ">>> FusionEKF::ProcessMeasurement" << endl;
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      
      cout << "Init RADAR" << endl;
      
      // TODO: Convert radar from polar to cartesian coordinates and initialize state.
      float ro = measurement_pack.raw_measurements_(0); 	// meas_rho
      float phi = measurement_pack.raw_measurements_(1);	// meas_phi
      float ro_dot = measurement_pack.raw_measurements_(2);	// meas_rho_dot
      
      ekf_.x_ << ro * cos(phi), ro * sin(phi), ro_dot * cos(phi), ro_dot * sin(phi);					// px, py, vx, vy
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Init LASER" << endl;
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0; 	// px, py, vx, vy
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    
    
    // P is the state covariance matrix, which contains information about the uncertainty of the object's position and velocity
    P_ = MatrixXd(4, 4);
    P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

    // The State Transition Matrix
    // x'= Fx + noise
    F_ = MatrixXd(4, 4);
    F_ << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

    ekf_.Init(ekf_.x_, /*x_in*/
              P_, /*P_in*/
              F_, /*F_in*/
              H_laser_, 
              Hj_, 
              R_laser_, 
              R_radar_, 
              ekf_.Q_); /*Q_in*/
    
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
  //set the acceleration noise components
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Process Covariance Matrix Q
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0, 
  			0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay, 
  			dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0, 
  			0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;

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
    // TODO: Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  
  cout << "<<< FusionEKF::ProcessMeasurement" << endl;
}
