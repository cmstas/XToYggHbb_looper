#include <math.h>

float deltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}
float deltaR2(float eta1, float phi1, float eta2, float phi2) {
  float deta = std::abs(eta1-eta2);
  float dphi = deltaPhi(phi1,phi2);
  return deta*deta + dphi*dphi;
}
float deltaR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}


