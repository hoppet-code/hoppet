#include <iostream>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;


int main()
{
using namespace std;
using namespace GiNaC;
ex x1 = numeric(8,3);
ex x2 = numeric(1,5);
ex i1 = 1;
ex i2 = 1;
cout << "Li_{1,1}(8/3,1/5) = "
<< evalf(zeta({3,1,3,1})) << endl;
return 0;
}