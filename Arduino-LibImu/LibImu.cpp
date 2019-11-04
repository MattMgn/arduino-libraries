#include <LibImu.hpp>

/*
 * ************************************************************************************************
 * Circular calibration of a single component sensor
 * ************************************************************************************************
 */

void LibImuCircularCalibInit(struct CircularCalib *calib, float timeout)
{
    for (int i = 0; i < 3; i++) {
        calib->min[i] = 0.0;
        calib->max[i] = 0.0;
    }
    calib->timer = 0.0;
    calib->timeout = timeout;
    calib->sample_nb = 0;
}

void LibImuCircularCalibProcess(struct CircularCalib *calib, const float *input, float dt)
{
    float tmp[3];

    /* copy values to avoid overwritting */
    memcpy(tmp, input, sizeof(tmp));

    if (calib->sample_nb == 0) {
        /* initialize values */
        for (int j = 0; j < 3; j++) {
            calib->max[j] = tmp[j];
            calib->min[j] = tmp[j];
        }
    } else {
        /* keep extremum values */
        for (int j = 0; j < 3; j++) {
            if(tmp[j] > calib->max[j])
                calib->max[j] = tmp[j];
            if(tmp[j] < calib->min[j])
                calib->min[j] = tmp[j];
        }
    }

    calib->timer += dt;
    calib->sample_nb++;
}

void LibImuCircularCalibComputeBiasScale(struct CircularCalib *calib, float *bias, float *scale)
{
    float avg_rad;

    for (int i = 0; i < 3; i++) {
        /* Get hard iron correction */
        bias[i]  = (calib->max[i] + calib->min[i]) / 2.0;

        /* Get soft iron correction estimate */
        scale[i]  = (calib->max[i] - calib->min[i]) / 2.0;
    }

    avg_rad = (scale[0] + scale[1] + scale[2]) / 3.0;

    scale[0] /= avg_rad;
    scale[1] /= avg_rad;
    scale[2] /= avg_rad;
}


/*
 * ************************************************************************************************
 * Spherical calibration of a three components sensor
 * ************************************************************************************************
 */



void LibImuSphericalCalibFeedVectors(struct SphericalCalib *calib, float x, float y, float z)
{

    if (calib->sample_nb >= SPHERICAL_MAX_SAMPLE_NB - 1)
        return;

    calib->sample_nb++;

    for (int i = 0; i < 3; i++) {
        /* Simple vector */
        calib->x_vec[calib->sample_nb] = x;
        calib->y_vec[calib->sample_nb] = y;
        calib->z_vec[calib->sample_nb] = z;

        calib->s_vec[calib->sample_nb] = x*x + y*y +z*z;
        calib->A_matrix[i][0] = 2.0 * x;
        calib->A_matrix[i][1] = 2.0 * y;
        calib->A_matrix[i][2] = 2.0 * z;
        calib->A_matrix[i][3] = 1.0;

    }
}

void LibImuSphericalCalibInit(struct SphericalCalib *calib)
{
    int err;
    int nrows = calib->sample_nb;
    int ncols = 4;

    double* dummy_array;
    double U[nrows][ncols];
    double V[ncols][ncols];
    double singular_values[ncols];
    double Astar[ncols][nrows];

    dummy_array = (double*) malloc(ncols * sizeof(double)); 

    Singular_Value_Decomposition((double*) calib->A_matrix, nrows, ncols, (double*)U,
                                  singular_values, (double*)V, dummy_array);

    Singular_Value_Decomposition_Inverse((double*)U, singular_values, (double*)V,
                                         singular_values[0] * DBL_EPSILON * ncols, nrows, ncols,
                                         (double*) Astar);

}