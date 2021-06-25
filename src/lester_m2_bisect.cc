#include "SVJ/Production/interface/lester_mt2_bisect.h"

namespace Lester {
bool ellipsesAreDisjoint(const EllipseParams & e1, const EllipseParams & e2) {
  /* We want to construct the polynomial "Det(lamdba A + B)" where A and B are the 3x3 matrices associated with e1 and e2, and we want to get that
  polynomial in the form lambda^3 + a lambda^2 + b lambda + c.


  Note that by default we will not have unity as the coefficient of the lambda^3 term, however the redundancy in the parametrisation of A and B allows us to scale the whole ply until the first term does have a unit coefficient.
  */

  if (e1==e2) {
    return false; // Probably won't catch many cases, but may as well have it here.
  }

  // first get unscaled terms:
  const double coeffLamPow3 = e1.det; // Note that this is the determinant of the symmetric matrix associated with e1.
  const double coeffLamPow2 = e1.lesterFactor(e2);
  const double coeffLamPow1 = e2.lesterFactor(e1);
  const double coeffLamPow0 = e2.det; // Note that this is the determinant of the symmetric matrix associated with e2.

  // Since question is "symmetric" and since we need to dovide by coeffLamPow3 ... do this the way round that involves dividing by the largest number:

  if (fabs(coeffLamPow3) >= fabs(coeffLamPow0)) {
    return __private_ellipsesAreDisjoint(coeffLamPow3, coeffLamPow2, coeffLamPow1, coeffLamPow0); // normal order
  } else {
    return __private_ellipsesAreDisjoint(coeffLamPow0, coeffLamPow1, coeffLamPow2, coeffLamPow3); // reversed order
  }
}
bool __private_ellipsesAreDisjoint(const double coeffLamPow3, const double coeffLamPow2, const double coeffLamPow1, const double coeffLamPow0) {

  // precondition of being called:
  //assert(fabs(coeffLamPow3)>=fabs(coeffLamPow0));

  if(coeffLamPow3==0) {
    // The ellipses were singular in some way.
    // Cannot determine whether they are overlapping or not.
    throw 1;
  }

  // now scale terms to monomial form:
  const double a = coeffLamPow2 / coeffLamPow3;
  const double b = coeffLamPow1 / coeffLamPow3;
  const double c = coeffLamPow0 / coeffLamPow3;

#ifdef LESTER_DEEP_FIDDLE
  {
    const double thing1 = -3.0*b + a*a;
    const double thing2 = -27.0*c*c + 18.0*c*a*b + a*a*b*b - 4.0*a*a*a*c - 4.0*b*b*b;
    std::cout 
      << (thing1>0) << " && " << (thing2>0) << " && [[ " << (a>=0) << " " << (3.0*a*c + b*a*a - 4.0*b*b<0)  << " ] or "
      << "[ " << (a< 0)   << " ] =("<< ((a >= 0 /*&& thing1 > 0*/ && 3.0*a*c + b*a*a - 4.0*b*b< 0 /*&& thing2 > 0*/) ||
                                 (a <  0 /*&& thing1 > 0*/                                 /*&& thing2 > 0*/)) << ")] " << (
          ( (a >= 0 && thing1 > 0 && 3.0*a*c + b*a*a - 4.0*b*b< 0 && thing2 > 0) ||
                                 (a <  0 && thing1 > 0                                 && thing2 > 0))
          
          ) << std::endl;
  }
#endif

  // Use the main result of the above paper:
  const double thing1 = -3.0*b + a*a;
  if (thing1<=0) return false;
  const double thing2 = -27.0*c*c + 18.0*c*a*b + a*a*b*b - 4.0*a*a*a*c - 4.0*b*b*b;
  if (thing2<=0) return false;

  // ans true means ellipses are disjoint:
  const bool ans = ( (a >= 0 /*&& thing1 > 0*/ && 3.0*a*c + b*a*a - 4.0*b*b< 0 /*&& thing2 > 0*/) ||
                     (a <  0 /*&& thing1 > 0*/                                 /*&& thing2 > 0*/));
  return ans;

}

}

void myversion(){

  std::cout << "Version is : 2014_11_13" << std::endl;

}

