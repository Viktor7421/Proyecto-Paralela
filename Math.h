#pragma once
#include "./Globals.h"

//-------------------------------------------------------------------------
// Return (a*s + c) MOD m; a, s, c and m must be < 2^35
//
double MultModM(double a, double s, double c, double m) {
  double v;
  long a1;

  v = a * s + c;

  if (v >= two53 || v <= -two53) {
    a1 = static_cast<long>(a / two17);
    a -= a1 * two17;
    v = a1 * s;
    a1 = static_cast<long>(v / m);
    v -= a1 * m;
    v = v * two17 + a * s + c;
  }

  a1 = static_cast<long>(v / m);
  /* in case v < 0)*/
  if ((v -= a1 * m) < 0.0)
    return v += m;
  else
    return v;
}


//-------------------------------------------------------------------------
// Compute the vector v = A*s MOD m. Assume that -m < s[i] < m.
// Works also when v = s.
//
void MatVecModM(const double A[3][3], const double s[3], double v[3],
                double m) {
  int i;
  double x[3]; // Necessary if v = s

  for (i = 0; i < 3; ++i) {
    x[i] = MultModM(A[i][0], s[0], 0.0, m);
    x[i] = MultModM(A[i][1], s[1], x[i], m);
    x[i] = MultModM(A[i][2], s[2], x[i], m);
  }
  for (i = 0; i < 3; ++i)
    v[i] = x[i];
}


//-------------------------------------------------------------------------
// Compute the matrix C = A*B MOD m. Assume that -m < s[i] < m.
// Note: works also if A = C or B = C or A = B = C.
//
void MatMatModM(const double A[3][3], const double B[3][3], double C[3][3],
                double m) {
  int i, j;
  double V[3], W[3][3];

  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j)
      V[j] = B[j][i];
    MatVecModM(A, V, V, m);
    for (j = 0; j < 3; ++j)
      W[j][i] = V[j];
  }
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      C[i][j] = W[i][j];
}


//-------------------------------------------------------------------------
// Compute the matrix B = (A^(2^e) Mod m);  works also if A = B.
//
void MatTwoPowModM(const double A[3][3], double B[3][3], double m, long e) {
  int i, j;

  /* initialize: B = A */
  if (A != B) {
    for (i = 0; i < 3; ++i)
      for (j = 0; j < 3; ++j)
        B[i][j] = A[i][j];
  }
  /* Compute B = A^(2^e) mod m */
  for (i = 0; i < e; i++)
    MatMatModM(B, B, B, m);
}


//-------------------------------------------------------------------------
// Compute the matrix B = (A^n Mod m);  works even if A = B.
//
void MatPowModM(const double A[3][3], double B[3][3], double m, long n) {
  int i, j;
  double W[3][3];

  /* initialize: W = A; B = I */
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j) {
      W[i][j] = A[i][j];
      B[i][j] = 0.0;
    }
  for (j = 0; j < 3; ++j)
    B[j][j] = 1.0;

  /* Compute B = A^n mod m using the binary decomposition of n */
  while (n > 0) {
    if (n % 2)
      MatMatModM(W, B, B, m);
    MatMatModM(W, W, W, m);
    n /= 2;
  }
}
