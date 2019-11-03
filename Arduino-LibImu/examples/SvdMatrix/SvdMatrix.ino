#include <LibImu.hpp>

#define FLOAT_PRECISION 6

double A[12] ={1.0, 2.0, 3.0,
               4.0, 5.0, 6.0,
               7.0, 8.0, 9.0,
               10.0, 11.0, 12.0};

int nrows = 4;
int ncols = 3;

double U[12];
double D[12];
double V[12];
double dummy_array[12];

double data[3000] = {0};

double Astar[12];

void setup() {
    Serial.begin(115200);
    Serial.println("SVD test");

    /* Compute SVD of A */
    Singular_Value_Decomposition(A, nrows, ncols, U, D, V, dummy_array);

    Serial.println("U = ");
    for (int i = 0; i < nrows; ++i ) {
        for (int j = 0; j < nrows; ++j) {
            Serial.print(U[nrows * i  + j], FLOAT_PRECISION);
            Serial.print(" | ");
        }
        Serial.println();
    }

    Serial.println("V = ");
    for (int i = 0; i < ncols; ++i ) {
        for (int j = 0; j < ncols; ++j) {
            Serial.print(V[ncols * i  + j], FLOAT_PRECISION);
            Serial.print(" | ");
        }
        Serial.println();
    }

    Serial.println("D = ");
    for (int i = 0; i < nrows; ++i ) {
        for (int j = 0; j < ncols; ++j) {
            Serial.print(D[ncols * i  + j], FLOAT_PRECISION);
            Serial.print(" | ");
        }
        Serial.println();
    }

    /* Compute Pseudo Inverse of A called Astar */
    Singular_Value_Decomposition_Inverse(U, D, V, D[0] * DBL_EPSILON * ncols, nrows, ncols, Astar);

    Serial.println("Astar = ");
    for (int i = 0; i < nrows; ++i ) {
        for (int j = 0; j < ncols; ++j) {
            Serial.print(Astar[ncols * i  + j]);
            Serial.print(" | ");
        }
        Serial.println();
    }


}

void loop() {

}
