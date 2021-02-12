#pragma once

#include "BaseQuadEstimator.h"
#include "matrix/math.hpp"
#include "Math/Quaternion.h"

using matrix::Vector;
using matrix::Matrix;
using matrix::SquareMatrix;

#include "Eigen/Dense"
#include "Eigen/SVD"
using Eigen::MatrixXf;
using Eigen::VectorXf;

struct PredictExtArg
{
    float dt;
    float rollEst;
    float pitchEst;
    V3F accel;
    V3F gyro;
};

class QuadEstimatorEKF : public BaseQuadEstimator
{
public:
    ///*Constructor.
    QuadEstimatorEKF(string config, string name);

    ///*/Destructor.
    virtual ~QuadEstimatorEKF();

    ///*Initialization of kalman filter data.
    virtual void Init();

    // Prediction Predicts the state and the state covariance.
    // using the process model.
    // @param dt Time between k and k+1 in s.
    // @param accel Accelerometer readings.
    // @param gyro Gayeroscope readings.
    virtual void Predict(float dt, V3F accel, V3F gyro);
    virtual void UpdateFromIMU(V3F accel, V3F gyro);
    virtual void UpdateFromGPS(V3F pos, V3F vel);
    virtual void UpdateFromBaro(float z) {};
    virtual void UpdateFromMag(float magYaw);

    //Number of Augemented states
    static const int QUAD_EKF_NUM_STATES = 7;

    // process covariance
    MatrixXf Q;

    // GPS measurement covariance
    MatrixXf R_GPS;

    // Magnetometer measurement covariance
    MatrixXf R_Mag;

    // Extra Input for Prediction step
    struct PredictExtArg predict_inputs;

    // attitude filter state
    float pitchEst, rollEst;
    float accelPitch, accelRoll; // raw pitch/roll angles as calculated from last accelerometer.. purely for graphing.
    V3F accelG;
    V3F lastGyro;

    // EKF state and covariance
    VectorXf ekfState;
    MatrixXf ekfCov;

    // params
    float attitudeTau;
    float dtIMU;

    // Access functions for graphing variables
    virtual bool GetData(const string& name, float& ret) const;
    virtual vector<string> GetFields() const;
    string _name;

    // error vs ground truth (trueError = estimated-actual)
    virtual void UpdateTrueError(V3F truePos, V3F trueVel, SLR::Quaternion<float> trueAtt);
    VectorXf trueError;
    float pitchErr, rollErr, maxEuler;

    float posErrorMag, velErrorMag;

    // EstimatedPosition estimated position getter.
    virtual V3F EstimatedPosition() 
    {
        return V3F(ekfState(0), ekfState(1), ekfState(2));
    }

    // EstimatedVelocity estimated velocity getter.
    virtual V3F EstimatedVelocity()
    {
        return V3F(ekfState(3), ekfState(4), ekfState(5));
    }

    // EstimatedAttitude estimated Attitude getter.
    virtual Quaternion<float> EstimatedAttitude()
    {
        return Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, ekfState(6));
    }

    // EstimatedOmega estimated Omega getter.
    virtual V3F EstimatedOmega()
    {
        return lastGyro;
    }

    // Calculate the condition number of the EKF ovariance matrix
    float CovConditionNumber() const;


    // calculates the mean state vector based dynamic model.
    //@param x The state vector.
    //@param[in] kd an object contains all kalman data {KalmanData}.
    //@return the mean state vector.
    static Eigen::VectorXf gFun(const VectorXf &mean, const void *p_args=NULL);

    // calculates the derivative of g_function.
    //@param x The state vector.
    //@param[in] kd an object contains all kalman data {KalmanData}.
    //@return the state transition matrix.
    static Eigen::MatrixXf g_prime(const Eigen::VectorXf &mean, const void *p_args=NULL);
};
