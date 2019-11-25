#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  
  
  ekf_.P_ = MatrixXd(4,4);

  ekf_.P_ << 1, 0, 0,    0,
	     0, 1, 0,    0,
	     0, 0, 1000, 0,
	     0, 0, 0,    1000;
  
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  noise_ax = 9;
  noise_ay = 9;


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
       
    cout<< "Initializing EKF......................."<<endl;
    ekf_.x_ = VectorXd(4);
    
    
    //first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {


      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */      
      float rho = measurement_pack.raw_measurements_[0];
      float phi =  measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];     
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      
      /**
      Initialize state.
      */ 
      ekf_.x_ << px, py, vx, vy;
      
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */      
      
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1),0.0,0.0;
    }   
    
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;

    // first measurement is used to initialize. No update is done
    return;
    
    }




  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  if( dt > 0.001 ) 
  { 
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0,  dt,
             0, 0, 1,  0,
             0, 0, 0,  1;
  
  float dt2 = dt*dt;
  float dt3 = dt2*dt;
  float dt4 = dt3*dt;      
  
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ <<  (dt4/4)*noise_ax, 0, (dt3/2)*noise_ax, 0,
              0, (dt4/4)*noise_ay, 0, (dt3/2)*noise_ay,
              (dt3/2)*noise_ax, 0, dt2*noise_ax, 0,
              0, (dt3/2)*noise_ay, 0, dt2*noise_ay;
  

  // make prediction  
  ekf_.Predict();

  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);    
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);   


  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;

}
