#include "kalman_filter.h"

#include <iostream>

using namespace std;

#define PI 3.14159265

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in_laser, MatrixXd &H_in_radar, MatrixXd &R_in_laser,MatrixXd &R_in_radar, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser_ = H_in_laser;
  H_radar_ = H_in_radar;
  R_laser_ = R_in_laser;
  R_radar_ = R_in_radar;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  
  cout << ">>> KalmanFilter::Predict" << endl;
  
  // print the state
  cout << "KalmanFilter::Predict start state" << endl;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  cout << "F_ = " << F_ << endl;
  cout << "Q_ = " << Q_ << endl;
  
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;			// P′= FPFT + Q
  
  // print the state
  cout << "KalmanFilter::Predict end state" << endl;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  
  cout << "<<< KalmanFilter::Predict" << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  cout << ">>> KalmanFilter::Update" << endl;
  
  // print the state
  cout << "KalmanFilter::Update start state" << endl;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;						// y = z − Hx′
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
  
  // print the state
  cout << "KalmanFilter::Update end state" << endl;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  
  cout << "<<< KalmanFilter::Update" << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  cout << ">>> KalmanFilter::UpdateEKF" << endl;
  
  // print the state
  cout << "KalmanFilter::UpdateEKF start state" << endl;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */  
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  // If rho == 0, skip the update step to avoid dividing by zero.
  // This is crude but should be fairly robust on our data set.
  if( px == 0. && py == 0. )
    return;
  
  H_radar_ = tools.CalculateJacobian( x_ );
  VectorXd hofx(3);
  float rho = sqrt( px * px + py * py );
  hofx << rho, atan2( py, px ), ( px * vx + py * vy ) / rho;
  
  // Update the state using Extended Kalman Filter equations
  VectorXd y = z - hofx;
  if( y[1] > PI )
    y[1] -= 2.f * PI;
  if( y[1] < -PI )
    y[1] += 2.f * PI;
  
  MatrixXd Hj_transpose = H_radar_.transpose();
  MatrixXd S = H_radar_ * P_ * Hj_transpose + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Hj_transpose * Si;

  // Compute new state
  x_ = x_ + ( K * y );
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = ( I - K * H_radar_ ) * P_;  
  
  // print the state
  cout << "KalmanFilter::UpdateEKF end state" << endl;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  
  cout << "<<< KalmanFilter::UpdateEKF" << endl;
}
