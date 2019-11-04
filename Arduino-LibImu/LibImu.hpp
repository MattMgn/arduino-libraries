#ifndef LibImu_hpp
#define LibImu_hpp

#include <string.h>
#include <LinearAlgebraSvd.hpp>

#define SPHERICAL_MAX_SAMPLE_NB     2000

/*
 * ************************************************************************************************
 * Circular calibration of a single component sensor
 * ************************************************************************************************
 */

struct CircularCalib {
    float max[3];
    float min[3];
    float timer;
    float timeout;
    int sample_nb;
};

void LibImuCircularCalibInit(struct CircularCalib *calib, float timeout);
void LibImuCircularCalibProcess(struct CircularCalib *calib, const float *input, float dt);
void LibImuCircularCalibComputeBiasScale(struct CircularCalib *calib, float *bias, float *scale);

/*
 * ************************************************************************************************
 * Spherical calibration of a three components sensor
 * ************************************************************************************************
 */

struct SphericalCalib {
    float x_vec[SPHERICAL_MAX_SAMPLE_NB];
    float y_vec[SPHERICAL_MAX_SAMPLE_NB];
    float z_vec[SPHERICAL_MAX_SAMPLE_NB];
    float s_vec[SPHERICAL_MAX_SAMPLE_NB];
    float A_matrix[SPHERICAL_MAX_SAMPLE_NB][4];
    float timeout;
    int sample_nb;
};

void LibImuSphericalCalibFeedVectors(struct SphericalCalib *calib, float x, float y, float z);


#endif