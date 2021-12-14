#pragma once
#include "./Globals.h"
#include "./Math.h"
#include <cstdlib>
#include <iostream>
#include <string>

class RngStream {
  using seed_t = unsigned long[6];

private:
  double Cg[6], Bg[6], Ig[6];
  bool anti, incPrec;
  std::string name;
  static double nextSeed[6];

  double U01() {
    long k;
    double p1, p2, u;

    p1 = a12 * Cg[1] - a13n * Cg[0];
    k = static_cast<long>(p1 / m1);
    p1 -= k * m1;
    if (p1 < 0.0)
      p1 += m1;
    Cg[0] = Cg[1];
    Cg[1] = Cg[2];
    Cg[2] = p1;

    p2 = a21 * Cg[5] - a23n * Cg[3];
    k = static_cast<long>(p2 / m2);
    p2 -= k * m2;
    if (p2 < 0.0)
      p2 += m2;
    Cg[3] = Cg[4];
    Cg[4] = Cg[5];
    Cg[5] = p2;

    /* Combination */
    u = ((p1 > p2) ? (p1 - p2) * norm : (p1 - p2 + m1) * norm);

    return (anti == false) ? u : (1 - u);
  }

  double U01d() {
    double u;
    u = U01();
    if (anti) {
      // Don't forget that U01() returns 1 - u in the antithetic case
      u += (U01() - 1.0) * fact;
      return (u < 0.0) ? u + 1.0 : u;
    } else {
      u += U01() * fact;
      return (u < 1.0) ? u : (u - 1.0);
    }
  }

public:
  RngStream(const char *name = "") : name(name) {
    anti = false;
    incPrec = false;

    for (int i = 0; i < 6; i++) {
      Bg[i] = Cg[i] = Ig[i] = nextSeed[i];
    }

    MatVecModM(A1p127, nextSeed, nextSeed, m1);
    MatVecModM(A2p127, &nextSeed[3], &nextSeed[3], m2);
  }

  double RandU01() {
    if (incPrec)
      return U01d();
    else
      return U01();
  }

  int RandInt(int low, int high) {
    return low + static_cast<int>((high - low + 1.0) * RandU01());
  }

  void GetState(seed_t seed) const {
    for (int i = 0; i < 6; ++i)
      seed[i] = static_cast<unsigned long>(Cg[i]);
  }

  void WriteState() {
    std::cout << "The current state of the Rngstream";
    if (name.size() > 0)
      std::cout << " " << name;
    std::cout << ":\n   Cg = { ";

    for (int i = 0; i < 5; i++) {
      std::cout << static_cast<unsigned long>(Cg[i]) << ", ";
    }
    std::cout << static_cast<unsigned long>(Cg[5]) << " }\n\n";
  }

  void WriteStateFull() {
    int i;

    std::cout << "The RngStream";
    if (name.size() > 0)
      std::cout << " " << name;
    std::cout << ":\n   anti = " << (anti ? "true" : "false") << "\n";
    std::cout << "   incPrec = " << (incPrec ? "true" : "false") << "\n";

    std::cout << "   Ig = { ";
    for (i = 0; i < 5; i++) {
      std::cout << static_cast<unsigned long>(Ig[i]) << ", ";
    }
    std::cout << static_cast<unsigned long>(Ig[5]) << " }\n";

    std::cout << "   Bg = { ";
    for (i = 0; i < 5; i++) {
      std::cout << static_cast<unsigned long>(Bg[i]) << ", ";
    }
    std::cout << static_cast<unsigned long>(Bg[5]) << " }\n";

    std::cout << "   Cg = { ";
    for (i = 0; i < 5; i++) {
      std::cout << static_cast<unsigned long>(Cg[i]) << ", ";
    }
    std::cout << static_cast<unsigned long>(Cg[5]) << " }\n\n";
  }

  void ResetStartStream() {
    for (int i = 0; i < 6; i++)
      Cg[i] = Bg[i] = Ig[i];
  }

  void ResetStartSubstream() {
    for (int i = 0; i < 6; i++)
      Cg[i] = Bg[i];
  }

  void ResetNextSubstream() {
    MatVecModM(A1p76, Bg, Bg, m1);
    MatVecModM(A2p76, &Bg[3], &Bg[3], m2);
    for (int i = 0; i < 6; i++)
      Cg[i] = Bg[i];
  }

  void AdvanceState(long e, long c) {
    double B1[3][3], C1[3][3], B2[3][3], C2[3][3];

    if (e > 0) {
      MatTwoPowModM(A1p0, B1, m1, e);
      MatTwoPowModM(A2p0, B2, m2, e);
    } else if (e < 0) {
      MatTwoPowModM(InvA1, B1, m1, -e);
      MatTwoPowModM(InvA2, B2, m2, -e);
    }

    if (c >= 0) {
      MatPowModM(A1p0, C1, m1, c);
      MatPowModM(A2p0, C2, m2, c);
    } else {
      MatPowModM(InvA1, C1, m1, -c);
      MatPowModM(InvA2, C2, m2, -c);
    }

    if (e) {
      MatMatModM(B1, C1, C1, m1);
      MatMatModM(B2, C2, C2, m2);
    }

    MatVecModM(C1, Cg, Cg, m1);
    MatVecModM(C2, &Cg[3], &Cg[3], m2);
  }

  bool SetSeed(const seed_t seed) {
    if (CheckSeed(seed)) {
      return false;
    }

    for (int i = 0; i < 6; i++)
      Cg[i] = Bg[i] = Ig[i] = seed[i];

    return true;
  }

  void IncreasedPresic(bool incp) { incPrec = incp; }

  void SetAntithetic(bool a) { anti = a; }

public:
  static bool SetPackageSeed(const seed_t seed) {
    if (CheckSeed(seed)) {
      return false;
    }

    for (int i = 0; i < 6; i++)
      nextSeed[i] = seed[i];

    return true;
  }

private:
  static int CheckSeed(const seed_t seed) {
    for (int i = 0; i < 3; i++) {
      if (seed[i] >= m1) {
        std::cerr << "****************************************\n"
                  << "ERROR: Seed[" << i << "] >= 4294967087, Seed is not set."
                  << "\n****************************************\n\n";
        return (-1);
      }
    }

    for (int i = 3; i < 6; i++) {
      if (seed[i] >= m1) {
        std::cerr << "****************************************\n"
                  << "ERROR: Seed[" << i << "] >= 4294967087, Seed is not set."
                  << "\n****************************************\n\n";
        return (-1);
      }
    }

    if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0) {
      std::cerr << "****************************\n"
                << "ERROR: First 3 seeds = 0.\n"
                << "****************************\n\n";
      return (-1);
    }

    if (seed[3] == 0 && seed[4] == 0 && seed[5] == 0) {
      std::cerr << "****************************\n"
                << "ERROR: Last 3 seeds = 0.\n"
                << "****************************\n\n";
      return (-1);
    }

    return 0;
  }
};

double RngStream::nextSeed[6] = {12345.0, 12345.0, 12345.0,
                                 12345.0, 12345.0, 12345.0};
