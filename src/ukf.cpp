#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .30; // Looks off

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .030; // Looks off

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;
  n_sig_ = 2 * n_aug_ + 1;
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  weights_.tail(2*n_aug_).fill(0.5/(lambda_ + n_aug_));
}


/**
* Destructor.
*/
UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    // first measurement
    float x_est;
    float y_est;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad


    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
            /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro_meas = meas_package.raw_measurements_(0);
      float theta_meas = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);
      x_est = ro_meas*cos(theta_meas);
      y_est =  ro_meas*sin(theta_meas);
      x_ << x_est, y_est, 0, 0, 0;
      //Errors adding in quadrature:
      double var_x_init = pow(std_radr_*cos(theta_meas), 2) + 
                    pow(ro_meas*sin(theta_meas)*std_radphi_, 2);
      double var_y_init = pow(std_radr_*sin(theta_meas), 2) + 
                    pow(ro_meas*cos(theta_meas)*std_radphi_, 2);
      P_ << var_x_init, 0, 0, 0, 0,
        0, var_y_init, 0, 0, 0,
        0, 0, 1000, 0, 0,
        0, 0, 0, 1000, 0,
        0, 0, 0, 0, 1000;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      //
      x_est = meas_package.raw_measurements_[0];
      y_est = meas_package.raw_measurements_[1];
      x_ << x_est, y_est, 0, 0, 0;
      P_ << (std_laspx_ * std_laspx_), 0, 0, 0, 0, 
            0, (std_laspy_ * std_laspy_), 0, 0, 0,
            0, 0, 1000, 0, 0,
            0, 0, 0, 1000, 0,
            0, 0, 0, 0, 1000;
    }
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;



  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! 
  Estimate the object's location. 
  Modify the state vector, x_. 
  Predict sigma points, 
  the state, 
  and the state covariance matrix.
  */

  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);

  x_aug.head(n_x_) = x_;

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;


  MatrixXd A = P_aug.llt().matrixL();
  double A_prefactor = sqrt(lambda_+n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.colwise() += x_aug;
  Xsig_aug.block(0,1,n_aug_,n_aug_) += A*A_prefactor;
  Xsig_aug.block(0,n_aug_+1,n_aug_,n_aug_) -= A*A_prefactor;

  VectorXd x_old = VectorXd(n_aug_);
  VectorXd dx = VectorXd(n_aug_);
  double v_k;
  double psi_k;
  double psi_dot_k;
  double nu_a_k;
  double nu_psi_k;
  //predict sigma points
  for (int i = 0; i < n_sig_; i++){
    x_old = Xsig_aug.col(i);
    v_k = x_old(2);
    psi_k = x_old(3);
    psi_dot_k = x_old(4);
    nu_a_k = x_old(5);
    nu_psi_k = x_old(6);
    //avoid division by zero
    if (std::abs(psi_dot_k) < 0.00001) {
      dx(0) = v_k * delta_t * cos(psi_k);
      dx(1) = v_k * delta_t * sin(psi_k);
    }
    else{
      dx(0) = v_k * (sin(psi_k + psi_dot_k*delta_t) - sin(psi_k)) / psi_dot_k;
      dx(1) = v_k * (-cos(psi_k + psi_dot_k*delta_t) + cos(psi_k)) / psi_dot_k;
    }
    double delta_t_2 = delta_t*delta_t;
    dx(0) += 0.5 * delta_t_2 * cos(psi_k) * nu_a_k;
    dx(1) += 0.5 * delta_t_2 * sin(psi_k) * nu_a_k;

    dx(2) = delta_t * nu_a_k;
    dx(3) = psi_dot_k * delta_t + 0.5 * delta_t_2 * nu_psi_k;
    dx(4) = delta_t * nu_psi_k;

    dx(5) = 0;
    dx(6) = 0;

    //write predicted sigma points into right column

    Xsig_pred_.col(i) = x_old.head(n_x_) + dx.head(n_x_);
  
  }

  for (int i = 0; i < n_x_; i++){
    x_(i) = Xsig_pred_.row(i).dot(weights_);
  }

  P_.fill(0.0);
  VectorXd resid = VectorXd(n_x_);
  for (int j = 0; j < 2*n_aug_+1; j++){
    resid = Xsig_pred_.col(j) - x_;
    P_ += weights_(j)*( resid * resid.transpose());
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z = meas_package.raw_measurements_;

  int n_z = 2;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  VectorXd x_old = VectorXd(n_x_);
  //transform sigma points into measurement space
  double rho_i;
  double phi_i;
  double rho_dot_i;
  
  double p_x_i;
  double p_y_i;
  double v_i;
  double psi_i;
  double psi_dot_i;

  for (int i = 0; i < n_sig_; i++){
    x_old = Xsig_pred_.col(i);

    p_x_i = x_old(0);
    p_y_i = x_old(1);
    
    Zsig.col(i) << p_x_i, p_y_i ;
  }
  //calculate mean predicted measurement
  for (int j = 0; j < n_z; j++){
    z_pred(j) = Zsig.row(j).dot(weights_);
  }

  //calculate measurement covariance matrix S
  S.fill(0.0);
  VectorXd resid = VectorXd(n_x_);
  for (int i = 0; i < n_sig_; i++){
    resid = Zsig.col(i) - z_pred;
    S += weights_(i) * ( resid * resid.transpose());
  }
  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;


//calculate cross correlation matrix
  VectorXd z_resid;
  VectorXd x_resid;
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++){
    x_resid = (Xsig_pred_.col(i) - x_);
    z_resid = (Zsig.col(i) - z_pred).transpose();
    x_resid(3) = fmod(x_resid(3) + M_PI, 2*M_PI) - M_PI;
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  z_resid = z - z_pred;
  x_ += K*z_resid;
  P_ += -K*S*K.transpose(); 

  //Normalized innovation squared
  double epsilon = z_resid.transpose() * S.inverse() * z_resid;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. 
  Modify the state vector, x_, and 
  covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;

  int n_z = 3;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  VectorXd x_old = VectorXd(n_x_);
  //transform sigma points into measurement space
  double rho_i;
  double phi_i;
  double rho_dot_i;
  
  double p_x_i;
  double p_y_i;
  double v_i;
  double psi_i;
  double psi_dot_i;

  for (int i = 0; i < n_sig_; i++){
    x_old = Xsig_pred_.col(i);
    p_x_i = x_old(0);
    p_y_i = x_old(1);
    v_i = x_old(2);
    psi_i = x_old(3);
    psi_dot_i = x_old(4);

    rho_i = sqrt((p_x_i*p_x_i) + (p_y_i*p_y_i));
    if (abs(rho_i) > 0.00001){
      phi_i = atan2(p_y_i, p_x_i);
      rho_dot_i = (p_x_i*cos(psi_i)*v_i + p_y_i*sin(psi_i)*v_i)/rho_i;
    } else {
      phi_i = 0.0; 
      rho_dot_i = 0.0; 
    }
    
    Zsig.col(i) << rho_i, phi_i, rho_dot_i ;
  }
  //calculate mean predicted measurement
  for (int j = 0; j < n_z; j++){
    z_pred(j) = Zsig.row(j).dot(weights_);
  }

  //calculate measurement covariance matrix S
  S.fill(0.0);
  VectorXd resid = VectorXd(n_x_);
  for (int i = 0; i < n_sig_; i++){
    resid = Zsig.col(i) - z_pred;
    resid(1) = fmod(resid(1)+M_PI, 2*M_PI) - M_PI;
    
    S += weights_(i) * ( resid * resid.transpose());
  }
  S(0,0) += std_radr_ * std_radr_;
  S(1,1) += std_radphi_ * std_radphi_;
  S(2,2) += std_radrd_ * std_radrd_;

//calculate cross correlation matrix
  VectorXd z_resid;
  VectorXd x_resid;
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(1.0);
  for (int i = 0; i < n_sig_; i++){
    x_resid = (Xsig_pred_.col(i) - x_);
    z_resid = (Zsig.col(i) - z_pred).transpose();
    x_resid(3) = fmod(x_resid(3) + M_PI, 2*M_PI) - M_PI;
    z_resid(1) = fmod(z_resid(1) + M_PI, 2*M_PI) - M_PI;
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  z_resid = z - z_pred;
  z_resid(1) = fmod(z_resid(1) + M_PI, 2*M_PI) - M_PI;
  x_ += K*z_resid;
  P_ += -K*S*K.transpose(); 

  //Normalized innovation squared
  double epsilon = z_resid.transpose() * S.inverse() * z_resid;

}
