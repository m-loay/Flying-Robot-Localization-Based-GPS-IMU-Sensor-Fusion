#include "Common.h"
#include "QuadEstimatorEKF.h"
#include "Utility/SimpleConfig.h"
#include "Utility/StringUtils.h"
#include "Math/Quaternion.h"
#include "UT.h"
#include "kalmanFilter.h"
#include "Math/Angles.h"
using namespace SLR;
static V g_(const V &mean, const void *p_args);
static M g_prime_ (const V &mean, const void *p_args);
const int QuadEstimatorEKF::QUAD_EKF_NUM_STATES;

QuadEstimatorEKF::QuadEstimatorEKF(string config, string name)
: BaseQuadEstimator(config),
    Q(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES),
    R_GPS(6, 6),
    R_Mag(1, 1),
    ekfState(QUAD_EKF_NUM_STATES),
    ekfCov(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES),
    trueError(QUAD_EKF_NUM_STATES)
{
    _name = name;
    Init();
}

QuadEstimatorEKF::~QuadEstimatorEKF()
{

}

void QuadEstimatorEKF::Init()
{
    ParamsHandle paramSys = SimpleConfig::GetInstance();

    paramSys->GetFloatVector(_config + ".InitState", ekfState);

    VectorXf initStdDevs(QUAD_EKF_NUM_STATES);
    paramSys->GetFloatVector(_config + ".InitStdDevs", initStdDevs);
    ekfCov.setIdentity();
    for (int i = 0; i < QUAD_EKF_NUM_STATES; i++)
    {
        ekfCov(i, i) = initStdDevs(i) * initStdDevs(i);
    }

    // complementary filter params
    attitudeTau = paramSys->Get(_config + ".AttitudeTau", .1f);
    dtIMU = paramSys->Get(_config + ".dtIMU", .002f);

    pitchEst = 0;
    rollEst = 0;

    // GPS measurement model covariance
    R_GPS.setZero();
    R_GPS(0, 0) = R_GPS(1, 1) = powf(paramSys->Get(_config + ".GPSPosXYStd", 0), 2);
    R_GPS(2, 2) = powf(paramSys->Get(_config + ".GPSPosZStd", 0), 2);
    R_GPS(3, 3) = R_GPS(4, 4) = powf(paramSys->Get(_config + ".GPSVelXYStd", 0), 2);
    R_GPS(5, 5) = powf(paramSys->Get(_config + ".GPSVelZStd", 0), 2);

    // magnetometer measurement model covariance
    R_Mag.setZero();
    R_Mag(0, 0) = powf(paramSys->Get(_config + ".MagYawStd", 0), 2);

    // load the transition model covariance
    Q.setZero();
    Q(0, 0) = Q(1, 1) = powf(paramSys->Get(_config + ".QPosXYStd", 0), 2);
    Q(2, 2) = powf(paramSys->Get(_config + ".QPosZStd", 0), 2);
    Q(3, 3) = Q(4, 4) = powf(paramSys->Get(_config + ".QVelXYStd", 0), 2);
    Q(5, 5) = powf(paramSys->Get(_config + ".QVelZStd", 0), 2);
    Q(6, 6) = powf(paramSys->Get(_config + ".QYawStd", 0), 2);
    Q *= dtIMU;

    rollErr = pitchErr = maxEuler = 0;
    posErrorMag = velErrorMag = 0;
}

void QuadEstimatorEKF::UpdateFromIMU(V3F accel, V3F gyro)
{
    // Improve a complementary filter-type attitude filter
    // 
    // Currently a small-angle approximation integration method is implemented
    // The integrated (predicted) value is then updated in a complementary filter style with attitude information from accelerometers
    // 
    // Implement a better integration method that uses the current attitude estimate (rollEst, pitchEst and ekfState(6))
    // to integrate the body rates into new Euler angles.
    //
    // HINTS:
    //  - there are several ways to go about this, including:
    //    1) create a rotation matrix based on your current Euler angles, integrate that, convert back to Euler angles
    //    OR 
    //    2) use the Quaternion<float> class, which has a handy FromEuler123_RPY function for creating a quaternion from Euler Roll/PitchYaw
    //       (Quaternion<float> also has a IntegrateBodyRate function, though this uses quaternions, not Euler angles)

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    // SMALL ANGLE GYRO INTEGRATION:
    // (replace the code below)
    // make sure you comment it out when you add your own code -- otherwise e.g. you might integrate yaw twice

    Quaternion<float> qt;
    qt = qt.FromEuler123_RPY(rollEst, pitchEst, ekfState(6));

    Quaternion<float> dq;
    dq = dq.IntegrateBodyRate(gyro, dtIMU);

    Quaternion<float> qt_bar (dq*qt);

    float predictedPitch = qt_bar.Pitch();
    float predictedRoll = qt_bar.Roll();
    ekfState(6) = qt_bar.Yaw();

    AngleNormF(ekfState(6));
    /////////////////////////////// END STUDENT CODE ////////////////////////////

    // CALCULATE UPDATE
    accelRoll = atan2f(accel.y, accel.z);
    accelPitch = atan2f(-accel.x, 9.81f);

    // FUSE INTEGRATION AND UPDATE
    rollEst = attitudeTau / (attitudeTau + dtIMU) * (predictedRoll)+dtIMU / (attitudeTau + dtIMU) * accelRoll;
    pitchEst = attitudeTau / (attitudeTau + dtIMU) * (predictedPitch)+dtIMU / (attitudeTau + dtIMU) * accelPitch;

    lastGyro = gyro;
}

void QuadEstimatorEKF::UpdateTrueError(V3F truePos, V3F trueVel, Quaternion<float> trueAtt)
{
    VectorXf trueState(QUAD_EKF_NUM_STATES);
    trueState(0) = truePos.x;
    trueState(1) = truePos.y;
    trueState(2) = truePos.z;
    trueState(3) = trueVel.x;
    trueState(4) = trueVel.y;
    trueState(5) = trueVel.z;
    trueState(6) = trueAtt.Yaw();

    trueError = ekfState - trueState;
    if (trueError(6) > F_PI) trueError(6) -= 2.f*F_PI;
    if (trueError(6) < -F_PI) trueError(6) += 2.f*F_PI;

    pitchErr = pitchEst - trueAtt.Pitch();
    rollErr = rollEst - trueAtt.Roll();
    maxEuler = MAX(fabs(pitchErr), MAX(fabs(rollErr), fabs(trueError(6))));

    posErrorMag = truePos.dist(V3F(ekfState(0), ekfState(1), ekfState(2)));
    velErrorMag = trueVel.dist(V3F(ekfState(3), ekfState(4), ekfState(5)));
}

void QuadEstimatorEKF::Predict(float dt, V3F accel, V3F gyro)
{
    //get the extra inputs
    predict_inputs.dt = dt;
    predict_inputs.accel = accel;
    predict_inputs.gyro = gyro;
    predict_inputs.pitchEst = pitchEst;
    predict_inputs.rollEst = rollEst;

    // Predict the current covariance forward by dt using the current accelerations and body rates as input.
    // INPUTS: 
    //   dt: time step to predict forward by [s]
    //   accel: acceleration of the vehicle, in body frame, *not including gravity* [m/s2]
    //   gyro: body rates of the vehicle, in body frame [rad/s]
    //   state (member variable): current state (state at the beginning of this prediction)
    //   
    // OUTPUT:
    //   update the member variable cov to the predicted covariance

    // HINTS
    // - update the covariance matrix cov according to the EKF equation.
    // 
    // - you may find the current estimated attitude in variables rollEst, pitchEst, state(6).
    //
    // - use the class MatrixXf for matrices. To create a 3x5 matrix A, use MatrixXf A(3,5).
    //
    // - the transition model covariance, Q, is loaded up from a parameter file in member variable Q
    // 
    // - This is unfortunately a messy step. Try to split this up into clear, manageable steps:
    //   1) Calculate the necessary helper matrices, building up the transition jacobian
    //   2) Once all the matrices are there, write the equation to update cov.
    //
    // - if you want to transpose a matrix in-place, use A.transposeInPlace(), not A = A.transpose()
    // 

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    kalmanFilter::predict(ekfState, ekfCov, Q, g_, g_prime_, &predict_inputs);


    /////////////////////////////// END STUDENT CODE ////////////////////////////
}

void QuadEstimatorEKF::UpdateFromGPS(V3F pos, V3F vel)
{
    VectorXf z(6), zFromX(6),Y(6);
    z(0) = pos.x;
    z(1) = pos.y;
    z(2) = pos.z;
    z(3) = vel.x;
    z(4) = vel.y;
    z(5) = vel.z;

    zFromX(0) = ekfState(0);
    zFromX(1) = ekfState(1);
    zFromX(2) = ekfState(2);
    zFromX(3) = ekfState(3);
    zFromX(4) = ekfState(4);
    zFromX(5) = ekfState(5);

    /**********************************************************************
     *    set output matrix H
     *********************************************************************/
    MatrixXf hPrime(6, QUAD_EKF_NUM_STATES);
    hPrime.setZero();
    hPrime.topLeftCorner(6 , 6) = M::Identity(6, 6);

    /**********************************************************************
     *    Calculate Innovation
     *********************************************************************/
    Y = z - zFromX;

    /**********************************************************************
     *    Calculate Kalman Gain
     *********************************************************************/
    M K = kalmanFilter::CalculateKalmanGain(ekfCov, hPrime, R_GPS);

    /**********************************************************************
     *    Update Linear
     *********************************************************************/
    kalmanFilter::update(ekfState, ekfCov, Y, hPrime, K);

}

void QuadEstimatorEKF::UpdateFromMag(float magYaw)
{
    Eigen::VectorXf Y(1);

    // MAGNETOMETER UPDATE
    // Hints: 
    //  - Your current estimated yaw can be found in the state vector: ekfState(6)
    //  - Make sure to normalize the difference between your measured and estimated yaw
    //    (you don't want to update your yaw the long way around the circle)
    //  - The magnetomer measurement covariance is available in member variable R_Mag
    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    /**********************************************************************
     *    set output matrix H
     *********************************************************************/
    Eigen::MatrixXf hPrime(1, 7);
    hPrime.setZero();
    hPrime(0,6) = 1;

    /**********************************************************************
     *    Calculate Innovation
     *********************************************************************/
    float diff_yaw = magYaw - ekfState(6);
    diff_yaw = AngleNormF(diff_yaw);
    Y(0) = diff_yaw;

    /**********************************************************************
     *    Calculate Kalman Gain
     *********************************************************************/
    M K = kalmanFilter::CalculateKalmanGain(ekfCov, hPrime, R_Mag);

    /**********************************************************************
     *    Update Linear
     *********************************************************************/
    kalmanFilter::update(ekfState, ekfCov, Y, hPrime, K);

    /////////////////////////////// END STUDENT CODE ////////////////////////////
}

// Execute an EKF update step
// z: measurement
// H: Jacobian of observation function evaluated at the current estimated state
// R: observation error model covariance 
// zFromX: measurement prediction based on current state
void QuadEstimatorEKF::Update(VectorXf& z, MatrixXf& H, MatrixXf& R, VectorXf& zFromX)
{
    assert(z.size() == H.rows());
    assert(QUAD_EKF_NUM_STATES == H.cols());
    assert(z.size() == R.rows());
    assert(z.size() == R.cols());
    assert(z.size() == zFromX.size());

    MatrixXf toInvert(z.size(), z.size());
    toInvert = H*ekfCov*H.transpose() + R;
    MatrixXf K = ekfCov * H.transpose() * toInvert.inverse();

    ekfState = ekfState + K*(z - zFromX);

    MatrixXf eye(QUAD_EKF_NUM_STATES, QUAD_EKF_NUM_STATES);
    eye.setIdentity();

    ekfCov = (eye - K*H)*ekfCov;
}

    // Calculate the condition number of the EKF ovariance matrix (useful for numerical diagnostics)
    // The condition number provides a measure of how similar the magnitudes of the error metric beliefs 
    // about the different states are. If the magnitudes are very far apart, numerical issues will start to come up.
    float QuadEstimatorEKF::CovConditionNumber() const
{
    MatrixXf m(7, 7);
    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            m(i, j) = ekfCov(i, j);
        }
    }

    Eigen::JacobiSVD<MatrixXf> svd(m);
    float cond = svd.singularValues()(0)
    / svd.singularValues()(svd.singularValues().size() - 1);
    return cond;
}

    // Access functions for graphing variables
    bool QuadEstimatorEKF::GetData(const string& name, float& ret) const
    {
    if (name.find_first_of(".") == string::npos) return false;
    string leftPart = LeftOf(name, '.');
    string rightPart = RightOf(name, '.');

    if (ToUpper(leftPart) == ToUpper(_name))
    {
        #define GETTER_HELPER(A,B) if (SLR::ToUpper(rightPart) == SLR::ToUpper(A)){ ret=(B); return true; }
        GETTER_HELPER("Est.roll", rollEst);
        GETTER_HELPER("Est.pitch", pitchEst);

        GETTER_HELPER("Est.x", ekfState(0));
        GETTER_HELPER("Est.y", ekfState(1));
        GETTER_HELPER("Est.z", ekfState(2));
        GETTER_HELPER("Est.vx", ekfState(3));
        GETTER_HELPER("Est.vy", ekfState(4));
        GETTER_HELPER("Est.vz", ekfState(5));
        GETTER_HELPER("Est.yaw", ekfState(6));

        GETTER_HELPER("Est.S.x", sqrtf(ekfCov(0, 0)));
        GETTER_HELPER("Est.S.y", sqrtf(ekfCov(1, 1)));
        GETTER_HELPER("Est.S.z", sqrtf(ekfCov(2, 2)));
        GETTER_HELPER("Est.S.vx", sqrtf(ekfCov(3, 3)));
        GETTER_HELPER("Est.S.vy", sqrtf(ekfCov(4, 4)));
        GETTER_HELPER("Est.S.vz", sqrtf(ekfCov(5, 5)));
        GETTER_HELPER("Est.S.yaw", sqrtf(ekfCov(6, 6)));

        // diagnostic variables
        GETTER_HELPER("Est.D.AccelPitch", accelPitch);
        GETTER_HELPER("Est.D.AccelRoll", accelRoll);

        GETTER_HELPER("Est.D.ax_g", accelG[0]);
        GETTER_HELPER("Est.D.ay_g", accelG[1]);
        GETTER_HELPER("Est.D.az_g", accelG[2]);

        GETTER_HELPER("Est.E.x", trueError(0));
        GETTER_HELPER("Est.E.y", trueError(1));
        GETTER_HELPER("Est.E.z", trueError(2));
        GETTER_HELPER("Est.E.vx", trueError(3));
        GETTER_HELPER("Est.E.vy", trueError(4));
        GETTER_HELPER("Est.E.vz", trueError(5));
        GETTER_HELPER("Est.E.yaw", trueError(6));
        GETTER_HELPER("Est.E.pitch", pitchErr);
        GETTER_HELPER("Est.E.roll", rollErr);
        GETTER_HELPER("Est.E.MaxEuler", maxEuler);

        GETTER_HELPER("Est.E.pos", posErrorMag);
        GETTER_HELPER("Est.E.vel", velErrorMag);

        GETTER_HELPER("Est.D.covCond", CovConditionNumber());
        #undef GETTER_HELPER
    }
    return false;
};

vector<string> QuadEstimatorEKF::GetFields() const
{
    vector<string> ret = BaseQuadEstimator::GetFields();
    ret.push_back(_name + ".Est.roll");
    ret.push_back(_name + ".Est.pitch");

    ret.push_back(_name + ".Est.x");
    ret.push_back(_name + ".Est.y");
    ret.push_back(_name + ".Est.z");
    ret.push_back(_name + ".Est.vx");
    ret.push_back(_name + ".Est.vy");
    ret.push_back(_name + ".Est.vz");
    ret.push_back(_name + ".Est.yaw");

    ret.push_back(_name + ".Est.S.x");
    ret.push_back(_name + ".Est.S.y");
    ret.push_back(_name + ".Est.S.z");
    ret.push_back(_name + ".Est.S.vx");
    ret.push_back(_name + ".Est.S.vy");
    ret.push_back(_name + ".Est.S.vz");
    ret.push_back(_name + ".Est.S.yaw");

    ret.push_back(_name + ".Est.E.x");
    ret.push_back(_name + ".Est.E.y");
    ret.push_back(_name + ".Est.E.z");
    ret.push_back(_name + ".Est.E.vx");
    ret.push_back(_name + ".Est.E.vy");
    ret.push_back(_name + ".Est.E.vz");
    ret.push_back(_name + ".Est.E.yaw");
    ret.push_back(_name + ".Est.E.pitch");
    ret.push_back(_name + ".Est.E.roll");

    ret.push_back(_name + ".Est.E.pos");
    ret.push_back(_name + ".Est.E.vel");

    ret.push_back(_name + ".Est.E.maxEuler");

    ret.push_back(_name + ".Est.D.covCond");

    // diagnostic variables
    ret.push_back(_name + ".Est.D.AccelPitch");
    ret.push_back(_name + ".Est.D.AccelRoll");
    ret.push_back(_name + ".Est.D.ax_g");
    ret.push_back(_name + ".Est.D.ay_g");
    ret.push_back(_name + ".Est.D.az_g");
    return ret;
};

/**
 * @brief g Function 
 *  which calculates the mean state vector based dynamic model.
 *
 * @param[in] mean
 *  the state vector {VectorXd&}.
 * 
 * @param[in] p_args
 *  Extra arguments {const void *}.
 * 
 * @return F.x
 *  the mean state vector {{VectorXd}}.
 */
static V g_(const V &mean, const void *p_args)
{
    assert(mean.size() == 7);
    V predictedState = V(mean.size());
    const struct PredictExtArg *p_predict = (struct PredictExtArg*)p_args;
    //getting the extra inputs
    const float dt(p_predict->dt);
    const float rollEst (p_predict->rollEst);
    const float pitchEst (p_predict->pitchEst);
    const V3F accel(p_predict->accel);
    const V3F gyro(p_predict->gyro);

    // Predict the current state forward by time dt using current accelerations and body rates as input
    // INPUTS: 
    //   curState: starting state
    //   dt: time step to predict forward by [s]
    //   accel: acceleration of the vehicle, in body frame, *not including gravity* [m/s2]
    //   gyro: body rates of the vehicle, in body frame [rad/s]
    //   
    // OUTPUT:
    //   return the predicted state as a vector

    // HINTS 
    // - dt is the time duration for which you should predict. It will be very short (on the order of 1ms)
    //   so simplistic integration methods are fine here
    // - we've created an Attitude Quaternion for you from the current state. Use 
    //   attitude.Rotate_BtoI(<V3F>) to rotate a vector from body frame to inertial frame
    // - the yaw integral is already done in the IMU update. Be sure not to integrate it again here

    Quaternion<float> attitude = Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, mean(6));

    ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

    predictedState(0) = mean(0) + mean(3)* dt;
    predictedState(1) = mean(1) + mean(4)* dt;
    predictedState(2) = mean(2) + mean(5)* dt;
    predictedState(3) = mean(3) + attitude.Rotate_BtoI(accel).x*dt;
    predictedState(4) = mean(4) + attitude.Rotate_BtoI(accel).y*dt;
    predictedState(5) = mean(5) - 9.81*dt + attitude.Rotate_BtoI(accel).z*dt;
    predictedState(6) = mean(6);
    /////////////////////////////// END STUDENT CODE ////////////////////////////

    return predictedState;
}

/**
 * @brief g_prime the derivative of g_function.
 *  In linear case it shall return the state transition Matrix.
 *  In non-linear it shall return the jacobians. 
 *
 * @param[in] mean
 *  the state vector {VectorXd&}.
 * 
 * @param[in] p_args
 *  Extra arguments {const void *}.
 * 
 * @return F 
 *  the state transition matrix {MatrixXd}.
 */
static M g_prime_ (const V &mean, const void *p_args)
{
    assert(mean.size() == 7);
    const struct PredictExtArg *p_predict = (struct PredictExtArg*)p_args;
    //getting the extra inputs
    const float dt(p_predict->dt);
    const float theta (p_predict->rollEst);
    const float phi (p_predict->pitchEst);
    const float psi (mean(6));
    const V3F accel(p_predict->accel);
    const V3F gyro(p_predict->gyro);

    // first, figure out the Rbg_prime
    M RbgPrime(3, 3);
    RbgPrime.setZero();
    RbgPrime(0,0) = -cos(theta)*sin(psi);
    RbgPrime(0,1) = -sin(phi)*sin(theta)*sin(psi) - cos(phi)*cos(psi);
    RbgPrime(0,2) = -cos(phi)*sin(theta)*sin(psi) + sin(phi)*cos(psi);
    
    RbgPrime(1,0) =  cos(theta)*cos(psi);
    RbgPrime(1,1) =  sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);
    RbgPrime(1,2) =  cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);

    //second get the new covar
    M gPrime(7, 7);
    gPrime.setIdentity();
    gPrime(0,3) = dt;
    gPrime(1,4) = dt;
    gPrime(2,5) = dt;

    gPrime(3,6) = dt*(accel.x * RbgPrime(0,0) + accel.y * RbgPrime(0,1) + accel.z * RbgPrime(0,2));
    gPrime(4,6) = dt*(accel.x * RbgPrime(1,0) + accel.y * RbgPrime(1,1) + accel.z * RbgPrime(1,2));
    gPrime(5,6) = dt*(accel.x * RbgPrime(2,0) + accel.y * RbgPrime(2,1) + accel.z * RbgPrime(2,2));   
    
    return gPrime;
}