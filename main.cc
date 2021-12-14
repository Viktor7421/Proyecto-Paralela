#include "./RngStream.h"
#include <iostream>
using namespace std;

auto main() -> int {
  unsigned long seed[6] = {327612383, 317095578,  14704821,
                           884064067, 1017894425, 16401881};
  RngStream::SetPackageSeed(seed);

  RngStream R;

  R.IncreasedPresic(true);
  R.SetAntithetic(true);

  cout << "Valor: " << R.RandU01() << endl;
  cout << "Valor: " << R.RandU01() << endl;
  cout << "Valor: " << R.RandU01() << endl;

  R.ResetNextSubstream();
  // R.AdvanceState(10, 6);
  cout << endl;

  cout << "Valor: " << R.RandU01() << endl;
  cout << "Valor: " << R.RandU01() << endl;
  cout << "Valor: " << R.RandU01() << endl;

  cout << endl;

  // R.WriteStateFull();
  return 0;
}
