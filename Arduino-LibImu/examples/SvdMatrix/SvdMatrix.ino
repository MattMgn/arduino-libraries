#include <LibImu.hpp>

#define FLOAT_PRECISION 6

#define M   2
#define N   3

double A[M][N] = {2.0, -1.0, 0.0,
                  4.0, 3.0, -2.0};

int nrows = M;
int ncols = N;

double U[M][N];
double V[N][N];
double singular_values[N];
double* dummy_array;
double Astar[N][M];

double data[3000] = {0.0};

void setup() {
    Serial.begin(115200);
    Serial.println("SVD test");

    Serial.println("A = ");
    for (int i = 0; i < nrows; ++i ) {
        for (int j = 0; j < ncols; ++j) {
            Serial.print(A[i][j], FLOAT_PRECISION);
            Serial.print(" | ");
        }
        Serial.println();
    }

    /* Compute SVD of A */
    dummy_array = (double*) malloc(N * sizeof(double)); 
    int err = Singular_Value_Decomposition((double*) A, nrows, ncols, (double*)U, singular_values, (double*)V, dummy_array);

    Serial.print("err = "); Serial.println(err);

    Serial.println("U = ");
    for (int i = 0; i < nrows; ++i ) {
        for (int j = 0; j < nrows; ++j) {
            Serial.print(U[i][j], FLOAT_PRECISION);
            Serial.print(" | ");
        }
        Serial.println();
    }

    Serial.println("V = ");
    for (int i = 0; i < ncols; ++i ) {
        for (int j = 0; j < ncols; ++j) {
            Serial.print(V[i][j], FLOAT_PRECISION);
            Serial.print(" | ");
        }
        Serial.println();
    }

    Serial.println("singular_values = ");
    for (int i = 0; i < ncols; ++i ) {
            Serial.print(singular_values[i], FLOAT_PRECISION);
            Serial.print(" | ");
    } Serial.println();

    /* Compute Pseudo Inverse of A called Astar */
    Singular_Value_Decomposition_Inverse((double*)U, singular_values, (double*)V, singular_values[0] * DBL_EPSILON * ncols, nrows, ncols,(double*) Astar);

    Serial.println("Astar = ");
    for (int i = 0; i < ncols; ++i ) {
        for (int j = 0; j < nrows; ++j) {
            Serial.print(Astar[i][j], FLOAT_PRECISION);
            Serial.print(" | ");
        }
        Serial.println();
    }

}

void loop() {

}
