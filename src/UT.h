/*
 * kf_lib.h
 *
 *  Created on: Apr 18, 2019
 *      Author: mody
 */
/** @file kf_lib.h
 *  @ingroup kf
 *  @brief kf class
 */

/**
 *  @addtogroup kf
 *  @{
 */

#ifndef UT_H_
#define UT_H_
#include <Eigen/Dense>
#include <vector>
#include <iostream>
using M = Eigen::MatrixXf;
using V = Eigen::VectorXf;
class UT
{
public:
    /**
     * @brief CalculateSigmaPoints
     * calculates the sigma points based on x,P.
     *
     * @param[in] mean 
     * The mean of stat vector  {VectorXd}.
     * 
     * @param[in] covariance
     * The mean of stat vector  {MatrixXd}.
     * 
	 * @return sigmaPoints 
     * the sigma points matrix  {MatrixXd}.
     * 
     */
    static M CalculateSigmaPoints(const V &mean,
                                  const M &covariance)
    {
        //extract sizes
        const int numState = mean.size();
        const int m_numSigmaPoints = 2*numState +1;
        const int lambda = 3 - numState;

        //create sigma points
        M sigmaPoints = M(numState,m_numSigmaPoints);

        // create square root matrix
        M L = covariance.llt().matrixL();

        // create augmented sigma points
        sigmaPoints.col(0)  = mean;
        for (int i = 0; i< numState; ++i) 
        {
            sigmaPoints.col(i+1)           = mean + sqrt(lambda+numState) * L.col(i);
            sigmaPoints.col(i+1+numState) = mean - sqrt(lambda+numState) * L.col(i);
        }
        return sigmaPoints;
    }

    /**
     * @brief CalculateSigmaPoints
     * calculates the sigma points based on x,P.
     *
     * @param[in] mean 
     * The mean of stat vector  {VectorXd}.
     * 
     * @param[in] covariance
     * The mean of stat vector  {MatrixXd}.
     * 
     * @param[in] Q 
     * the process noise matrix  {MatrixXd}.
     * 
	 * @return sigmaPoints 
     * the sigma points matrix  {MatrixXd}.
     * 
     */
    static M CalculateSigmaPoints(const V &mean,
                                  const M &covariance,
                                  const M &Q)
    {
        //extract sizes
        const int numState = mean.size();
        const int numAugState = numState + Q.rows();
        const int numSigmaPoints = 2*numAugState +1;
        const int lambda = 3 - numAugState;
        const int diff = numAugState - numState;

        // create augmented mean vector
        V x_aug = V(numAugState);
        x_aug.fill(0.0);
        x_aug.head(numState) = mean;

        // create augmented state covariance
        M P_aug = M(numAugState, numAugState);
        P_aug.fill(0.0);
        P_aug.topLeftCorner(numState, numState) = covariance;
        P_aug.bottomRightCorner(diff , diff) = Q;

        //create sigma points
        M sigmaPoints = M(numAugState, numSigmaPoints);

        // create square root matrix
        M L = P_aug.llt().matrixL();

        // create augmented sigma points
        sigmaPoints.col(0)  = x_aug;
        for (int i = 0; i< numAugState; ++i) 
        {
            sigmaPoints.col(i+1)             = x_aug + sqrt(lambda + numAugState) * L.col(i);
            sigmaPoints.col(i+1+numAugState) = x_aug - sqrt(lambda + numAugState) * L.col(i);
        }
        return sigmaPoints;
    }

    /**
     * @brief PredictSigmaPoints
     * Propagate sigma points into model(motion model/Sensor Model).
     *
     * @param[in] sigmaPoints 
     * is kalman data object  {MatrixXd}.
     * 
     * @param[in] model
     * The propagation model for sigma points  {V(KalmanDataUT &,V&)}.
     * 
	 * @param[in] sigmaPoints 
     * the kalman gain matrix  {MatrixXd}.
     * 
     * @param[in] p_args
     *  Extra arguments {const void *}.
     * 
	 * @param[in] na_aug 
     * number of augmented state {int}.
     * 
	 * @return predictedSigmaPoints 
     * the predicted sigma points {MatrixXd}.
     *
     */
    static M PredictSigmaPoints(const M &sigmaPoints, 
                                std::function<V(const V &,const void *)>model,
                                const void *p_args = NULL,
                                const int na_aug = 0)
    {
        //extract sizes
        const int numSigmaPoints = sigmaPoints.cols();
        const int numState = sigmaPoints.rows() - na_aug;
        const int colSize = sigmaPoints.rows();

        //create predicted sigma points matrix
        M predictedSigmaPoints = M(numState, numSigmaPoints);

        //create column vector to get sigma point col.
        V col(colSize);

        //propagate sigma pointion in dynamic model
        for (int i = 0; i < numSigmaPoints ; i++)
        {
            col = sigmaPoints.col(i);
            predictedSigmaPoints.col(i) = model(col,p_args);
        }
        return predictedSigmaPoints;
    }

    /**
     * @brief CalculateWeigts
     * It calculates the weights of sigma points to revive mean and covariance.
     *
	 * @param[in] numSigmaPoints 
     * number of sigma points {int}.
     * 
	 * @param[in] numAug 
     * number of augmented state {int}.
     * 
	 * @return weights 
     * the weights of sigma points {VectorXd}.
     *
     */
    static V CalculateWeigts(const int numSigmaPoints,
                             const int numAug = 0)
    {
        const int lambda = 3 - numAug;
        V weights = V(numSigmaPoints);
        weights.fill(1.0 / static_cast<double>(2.0 * static_cast<double>(lambda + numAug)));
        weights(0) = static_cast<double>(lambda) /static_cast<double>(lambda + numAug);
        return weights;
    }

    /**
     * @brief PredictMean
     * Revive the mean from the predicted sigma points.
     *
     * @param[in] predictedSigmaPoints 
     * the predicted sigma points  {MatrixXd}.
     * 
     * @param[in] weights 
     * the weights of sigma points  {VectorXd}.
     * 
	 * @return mean 
     * The mean of state vector {VectorXd}.
     *
     */
    static V PredictMean(const M &predictedSigmaPoints,
                         const V &weights)
    {
        //extract sizes
        const int numSigmaPoints = predictedSigmaPoints.cols();
        const int numState = predictedSigmaPoints.rows();

        //Initialize mean matrix
        V mean = V(numState);
        mean.fill(0.0);

        //calculate the mean
        for(int i=0 ; i<numSigmaPoints; i++)
        {
            mean = mean + weights(i) * predictedSigmaPoints.col(i);
        }
        return mean;
    }

    /**
     * @brief PredictCovariance
     * Revive the Covariance from the predicted sigma points.
     *
	 * @param[in] mean 
     * The mean of state vector {VectorXd}.
     * 
     * @param[in] predictedSigmaPoints 
     * the predicted sigma points  {MatrixXd}.
     * 
     * @param[in] weights 
     * the weights of sigma points  {VectorXd}.
     * 
     * @param[in] helperFunc
     * The helperFunc to perfrom specific operation on difference vector  {V(V&)}.
     *
     */
    static M PredictCovariance(const V &mean,
                               const M &predictedSigmaPoints,
                               const V &weights,
                               std::function<V(const V&)>helperFunc = NULL)
    {
        //extract sizes
        const int numSigmaPoints = predictedSigmaPoints.cols();
        const int numState = predictedSigmaPoints.rows();

        //Initialize covariance matrix
        M covariance = M(numState,numState);
        covariance.fill(0.0);
        
        // state difference
        V x_diff = V(numState);
        
        for(int i=0 ; i< numSigmaPoints; i++)
        {
            x_diff = predictedSigmaPoints.col(i);
            x_diff = x_diff - mean;

            if(helperFunc != NULL)
            {
                x_diff = helperFunc(x_diff);
            }
            else
            {
                /*do nothing*/
            }
            covariance = covariance + weights(i) * x_diff * x_diff.transpose();
        }
        return covariance;
    }

    /**
     * @brief TransformPredictedSigmaMeasurement
     * Transform predicted sigma points from state space to measurement space.
     *
     * @param[in] predictedSigmaPoints 
     * the predicted sigma points  {MatrixXd}.
     * 
	 * @param[in] meas_size 
     * number of measurement size {int}.
     * 
     * @param[in] model
     * The helperFunc to perfrom specific operation on difference vector  {V(V&)}.
     * 
     * @param[in] args
     *  Extra arguments {const void *}.
     * 
     * @return Zsig
     * The tranformed sigma points.
     *
     */
    static M TransformPredictedSigmaToMeasurement(const M &predictedSigmaPoints,
                                                  const int meas_size,
                                                  std::function<V(const V&,const void *args )>model,
                                                  const void *args = NULL)
    {
        //extract sizes
        const int numSigmaPoints = predictedSigmaPoints.cols();

        //create Predicted measurement matrix
        M Zsig = M(meas_size, numSigmaPoints);

        //create vectore measurement
        V col(meas_size);

        for (int i = 0; i < numSigmaPoints ; i++)
        {
            col = predictedSigmaPoints.col(i);
            Zsig.col(i) = model(col,args);
        }
        return Zsig;
    }

    /**
     * @brief CalculateKalmanGainUT
     * The calculation of kalman gain.
     *
     * @param[in,out] mean 
     * The mean of state vector  {VectorXd}.
     * 
     * @param[in,out] zpred 
     * The mean of measurement state vector  {VectorXd}.
     * 
     * @param[in] weights 
     * the weights of sigma points  {VectorXd}.
     * 
     * @param[in] xPredictedSigmaPoints 
     * the predicted sigma points  {MatrixXd}.
     * 
     * @param[in] zPredictedSigmaPoints 
     * the predicted sigma points of measurements  {MatrixXd}.
     * 
     * @param[in] S 
     * the sensor noise matrix  {MatrixXd}.
     * 
    * @param[in] xfun
     * The helperFunc to perfrom specific operation on difference vector  {V(V&)}.
     *
     * @param[in] zfun
     * The helperFunc to perfrom specific operation on difference vector for measurement  {V(V&)}.
     * 
	 * @return K 
     * the kalman gain matrix  {MatrixXd}.
     * 
     */
    static M CalculateKalmanGainUT(const V &mean,
                                   const V &zpred,
                                   const V &weights,
                                   const M &xPredictedSigmaPoints,
                                   const M &zPredictedSigmaPoints,
                                   const M &S,
                                   std::function<V(const V&)>xfun,
                                   std::function<V(const V&)>zfun)
    {
        //extract sizes
        const int numSigmaPoints = xPredictedSigmaPoints.cols();
        const int numState = mean.rows();
        const int numMeas = zpred.rows();

        // create matrix for cross correlation Tc
        M Tc = M(numState, numMeas);
        Tc.fill(0.0);
        
        // state difference
        V x_diff = V(numState);
        V z_diff = V(numMeas);
        
        for(int i=0 ; i<numSigmaPoints; i++)
        {
            x_diff = xPredictedSigmaPoints.col(i) - mean;

            if(xfun != NULL)
            {
                x_diff = xfun(x_diff);
            }
            else
            {
                /*do nothing*/
            }

            z_diff = zPredictedSigmaPoints.col(i) - zpred;

            if(zfun != NULL)
            {
                z_diff = zfun(z_diff);
            }
            else
            {
                /*do nothing*/
            }
            Tc = Tc + weights(i) * x_diff * z_diff.transpose();
        }
        // Kalman gain K;
        M K = Tc * S.inverse();
        return K;        
    }

    /**
     * @brief updateUT
     * Update step.
     *
     * @param[in,out] mean 
     * The mean of state vector  {VectorXd}.
     * 
     * @param[in,out] covariance
     * The mean of covariance matrix  {MatrixXd}.
     * 
     * @param[in] S 
     * the sensor noise matrix  {MatrixXd}.
     * 
	 * @param[in] K 
     * the kalman gain matrix  {MatrixXd}.
     */
    static void updateUT(V & mean,
                         M & covariance,
                         const V &Innovation,
                         const M &S,
                         const M &K)
    {
        // update state mean and covariance
        mean = mean + K * Innovation;
        covariance = covariance - K * S * K.transpose();        
    }
};

#endif /* UT_H_ */
/**
 *  @}
 */
