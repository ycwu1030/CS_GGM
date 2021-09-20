#ifndef SM_GAC_H
#define SM_GAC_H

#include <cmath>
#include <vector>

#ifndef alpha1
#define alpha1 127.9
#endif

#ifndef MZ
#define MZ 91.1876
#endif

#ifndef GF
#define GF 0.0000116637
#endif

// #ifndef Mt
// #define Mt 173.5
// #endif

// #ifndef Mh
// #define Mh 125.0
// #endif

#ifndef Pi
#define Pi 3.14159265358979312
#endif

double ScaleRange(double input, double min, double max);

class SM {
 public:
  SM();
  ~SM(){};

  double MW;
  double thetaW;
  double alpha;
  double vev;       // (sqrt(2)GF)^(-0.5)
  double ee;        // ee = sqrt(4*Pi*alpha)
  double g_weak;    // g = ee/sw
  double gp_hyper;  // g' = ee/cw
  double Yt;        // mt/vev;
};

// template <unsigned N5, unsigned N6, bool GEN = false>
class CSVMODEL : public SM {
 public:
  CSVMODEL(int N5, int N6, bool GEN = false);
  ~CSVMODEL(){};
  double GetLikelihood(double *cube);
  int Get_NRandoms() { return NRandoms; }
  int Get_NPARAMETERS() { return 2 * NVEVS + 7; }
  void Set_PARAMETERS(double *cube);
  void Print_PARAMETERS();

 private:
  void Random_to_vev_ui(double *cube);
  double Calc_Rho();
  void SetUpKappas();
  double SetUpChi2();
  int N5;
  int N6;
  int GEN;
  int NVEVS;
  int NRandoms;
  std::vector<double> vevs;
  std::vector<double> ui;
  std::vector<double> kfs;
  std::vector<double> kWs;
  std::vector<double> kZs;
  double kf, kW, kZ, kh;
  double lamfZ, kfZ, lamWZ;
};

class CSVGM : public SM {
 public:
  CSVGM(){};
  ~CSVGM(){};

  double GetLikelihood(double thetaH, double v1, double v2, double v3);
  void SetUpKappas();
  double vphi, vchi, vxi;
  double v1, v2, v3;
  double kfphi, kfchi, kfxi;
  double kWphi, kWchi, kWxi;
  double kZphi, kZchi, kZxi;
  double kf, kW, kZ, kh;
  double lamfZ, kfZ, lamWZ;
};

#endif  // End of SM_GAC_H
