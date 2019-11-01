#include <LinearAlgebraSvd.hpp>

#define FLOAT_PRECISION 6

double A[12] ={1.0, 2.0, 3.0,
               4.0, 5.0, 6.0,
               7.0, 8.0, 9.0,
               10.0, 11.0, 12.0};

int m = 4;
int n = 3;

double U[12];
double singular_values[12];
double V[12];
double dummy_array[12];

double data[3000] = {0};

void setup() {
    Serial.begin(115200);
    Serial.println("SVD test");

    // put your setup code here, to run once:
    Singular_Value_Decomposition(A, m, n, U, singular_values, V, dummy_array);

    Serial.println("U = ");
    for (int i = 0; i < 12; i++)
        Serial.println(U[i], FLOAT_PRECISION);

    Serial.println("V = ");
    for (int i = 0; i < 12; i++)
    for (int i = 0; i < 12; i++)
        Serial.println(V[i], FLOAT_PRECISION);

    Serial.println("D = ");
    for (int i = 0; i < 12; i++)
    for (int i = 0; i < 12; i++)
        Serial.println(singular_values[i], FLOAT_PRECISION);

}

void loop() {

}
