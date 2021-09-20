#include "SM_GAC_SEXPEN.h"

#include <iostream>

using namespace std;
#define EPS pow(10, -6)

SM::SM() {
  alpha = 1.0 / alpha1;
  ee = sqrt(4 * Pi * alpha);
  vev = pow(sqrt(2) * GF, -0.5);
  double A = sqrt(Pi * alpha) * vev;
  thetaW = asin(2 * A / MZ) / 2.0;
  MW = A / sin(thetaW);
  g_weak = ee / sin(thetaW);
  gp_hyper = ee / cos(thetaW);
  // Yt = Mt/vev;
#ifdef DEBUG
  cout << "A:  " << A << endl;
  cout << "vev: " << vev << endl;
  cout << "thetaW: " << thetaW << endl;
  cout << "MW: " << MW << endl;
  cout << "Yt: " << Yt << endl;
  cout << "g_weak: " << g_weak << endl;
  cout << "gp_hyper: " << gp_hyper << endl;
#endif
}

double CSVGM::GetLikelihood(double thetaH, double _v1, double _v2, double _v3) {
  double cH = cos(thetaH);
  double sH = sin(thetaH);
  vphi = vev * cH;
  vchi = vev * sH / 2.0 / sqrt(2);
  vxi = vchi;
  v1 = _v1;
  v2 = _v2;
  v3 = _v3;
  SetUpKappas();
  if (lamWZ > -0.69 || lamWZ < -1.44) return pow(-10, 99);
  double chi2 = 0;
  double dlamfZ = fabs(lamfZ) - 0.98;
  double dkfZ = fabs(kfZ) - 0.99;
  double dx[2] = {dlamfZ, dkfZ};
  double inverse_cov[2][2] = {{109.239, 29.4947}, {29.4947, 507.964}};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      chi2 += dx[i] * inverse_cov[i][j] * dx[j];
    }
  }
  chi2 = -chi2;
  return chi2;
}
void CSVGM::SetUpKappas() {
  kfphi = vphi / vev;
  kWphi = vphi / vev;
  kZphi = vphi / vev;
  kfchi = 0;
  kWchi = 2 * sqrt(2) * vchi / vev;
  kZchi = 4 * sqrt(2) * vchi / vev;
  kfxi = 0;
  kWxi = 4 * vxi / vev;
  kZxi = 0;
  kZ = v1 * kZphi + v2 * kZchi + v3 * kZxi;
  kW = v1 * kWphi + v2 * kWchi + v3 * kWxi;
  kf = v1 * kfphi + v2 * kfchi + v3 * kfxi;
  kh = sqrt(0.75 * kf * kf + 0.22 * kW * kW + 0.03 * kZ * kZ);
  lamfZ = kf / kZ;
  kfZ = kf * kZ / kh;
  lamWZ = kW / kZ;
}

const double C_W_5[3] = {4 * sqrt(2), 10 * sqrt(2), 12};
const double C_Z_5[3] = {16 * sqrt(2), 4 * sqrt(2), 0};
const double C_V_5[3] = {sqrt(2), sqrt(2), 1};
const double C_V_T_5 = 40;
const double C_W_6[3] = {5 * sqrt(2), 13 * sqrt(2), 17 * sqrt(2)};
const double C_Z_6[3] = {25 * sqrt(2), 9 * sqrt(2), sqrt(2)};
const double C_V_6[3] = {sqrt(2), sqrt(2), sqrt(2)};
const double C_V_T_6 = 70;

// template <unsigned N5, unsigned N6, bool GEN>
CSVMODEL::CSVMODEL(int _N5, int _N6, bool _GEN) : N5(_N5), N6(_N6), GEN(_GEN) {
  NVEVS = 1 + 3 * N5 + 3 * N6;
  // * First calculate the # of randoms for vevs
  NRandoms = N5 + N6;
  if (GEN) NRandoms = 3 * N5 + 3 * N6;
  // * Then the # of randoms for the unity vector
  NRandoms += 3 * N5 + 3 * N6;
  vevs = vector<double>(NVEVS, 0);
  ui = vector<double>(NVEVS, 0);
}

double ScaleRange(double input, double min, double max) {
  return min + (max - min) * input;
}

void SetSphere(const double total, double *cubics, const int dim,
               vector<double> &res) {
  // * dim is the dimension of res;
  res.clear();
  double tmp_cur = 0;
  double tmp_next = total;
  double angle_tmp;
  if (dim == 1) {
    res.push_back(total);
    return;
  }
  for (int i = 0; i < dim - 2; i++) {
    angle_tmp = ScaleRange(cubics[i], 0, M_PI);
    cubics[i] = angle_tmp;
    res.push_back(tmp_next * cos(angle_tmp));
    tmp_next = tmp_next * sin(angle_tmp);
  }
  angle_tmp = ScaleRange(cubics[dim - 2], -M_PI, M_PI);
  cubics[dim - 2] = angle_tmp;
  res.push_back(tmp_next * cos(angle_tmp));
  res.push_back(tmp_next * sin(angle_tmp));
}

void CSVMODEL::Random_to_vev_ui(double *cube) {
  // the number in cube[] are all in [0,1]
  double tmp_cur;
  double tmp_next = vev;
  double angle_tmp;
  double *cube_u;
  if (GEN) {
    vector<double> vgen;
    SetSphere(vev, cube, NVEVS, vgen);
    cube_u = cube + NVEVS - 1;
    vevs[0] = vgen[0];
    for (int i = 0; i < 3; i++) {
      // ! Setting Pentet vevs;
      for (int j = 0; j < N5; j++) {
        int id5 = 1 + i + 3 * j;
        vevs[id5] = vgen[id5] / sqrt(C_W_5[i] * C_V_5[i]);
      }
      // ! Setting Sextet vevs;
      for (int j = 0; j < N6; j++) {
        int id6 = 1 + 3 * N5 + i + 3 * j;
        vevs[id6] = vgen[id6] / sqrt(C_W_6[i] * C_V_6[i]);
      }
    }
  } else {
    vector<double> vac;
    SetSphere(vev, cube, N5 + N6 + 1, vac);
    cube_u = cube + N5 + N6;
    vevs[0] = vac[0];
    for (int i = 0; i < 3; i++) {
      // ! Setting Pentet vevs;
      for (int j = 0; j < N5; j++) {
        int id5 = 1 + i + 3 * j;
        vevs[id5] = vac[1 + j] / sqrt(C_V_T_5);
      }
      // ! Setting Sextet vevs;
      for (int j = 0; j < N6; j++) {
        int id6 = 1 + 3 * N5 + i + 3 * j;
        vevs[id6] = vac[1 + N5 + j] / sqrt(C_V_T_6);
      }
    }
  }

  SetSphere(1, cube_u, NVEVS, ui);
}

void CSVMODEL::SetUpKappas() {
  kfs.clear();
  kWs.clear();
  kZs.clear();
  kfs.push_back(vevs[0] / vev);
  kWs.push_back(vevs[0] / vev);
  kZs.push_back(vevs[0] / vev);
  int idv = 1;
  for (int i = 0; i < N5; i++) {
    for (int j = 0; j < 3; j++) {
      kfs.push_back(0);
      kWs.push_back(C_W_5[j] * vevs[idv] / vev);
      kZs.push_back(C_Z_5[j] * vevs[idv] / vev);
      ++idv;
    }
  }
  for (int i = 0; i < N6; i++) {
    for (int j = 0; j < 3; j++) {
      kfs.push_back(0);
      kWs.push_back(C_W_6[j] * vevs[idv] / vev);
      kZs.push_back(C_Z_6[j] * vevs[idv] / vev);
      ++idv;
    }
  }

  kf = 0;
  kW = 0;
  kZ = 0;
  for (int i = 0; i < NVEVS; i++) {
    kf += ui[i] * kfs[i];
    kW += ui[i] * kWs[i];
    kZ += ui[i] * kZs[i];
  }

  kh = sqrt(0.75 * kf * kf + 0.22 * kW * kW + 0.03 * kZ * kZ);
  lamfZ = kf / kZ;
  kfZ = kf * kZ / kh;
  lamWZ = kW / kZ;
}

double CSVMODEL::SetUpChi2() {
  double chi2 = 0;
  double dlamfZ = fabs(lamfZ) - 0.98;
  double dkfZ = fabs(kfZ) - 0.99;
  double dx[2] = {dlamfZ, dkfZ};
  double inverse_cov[2][2] = {{109.239, 29.4947}, {29.4947, 507.964}};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      chi2 += dx[i] * inverse_cov[i][j] * dx[j];
    }
  }
  chi2 = -chi2;
  return chi2;
}

double CSVMODEL::Calc_Rho() {
  double v2_w = 0;
  double v2_z = 0;
  v2_w += vevs[0] * vevs[0];
  v2_z += vevs[0] * vevs[0];
  int idv = 1;
  for (int i = 0; i < N5; i++) {
    for (int j = 0; j < 3; j++) {
      v2_w += vevs[idv] * vevs[idv] * C_W_5[j] * C_V_5[j];
      v2_z += vevs[idv] * vevs[idv] * C_Z_5[j] * C_Z_5[j];
      ++idv;
    }
  }
  for (int i = 0; i < N6; i++) {
    for (int j = 0; j < 3; j++) {
      v2_w += vevs[idv] * vevs[idv] * C_W_6[j] * C_V_6[j];
      v2_z += vevs[idv] * vevs[idv] * C_Z_6[j] * C_V_6[j];
      ++idv;
    }
  }
  return v2_w / v2_z;
}

double CSVMODEL::GetLikelihood(double *cube) {
  Random_to_vev_ui(cube);
  double rho = Calc_Rho();
  if (rho > 1.00134 || rho < 0.99944) return -pow(10, 99);
  SetUpKappas();
  if (lamWZ > -0.69 || lamWZ < -1.44) return -pow(10, 99);
  return SetUpChi2();
}

void CSVMODEL::Set_PARAMETERS(double *cube) {
  int id = NRandoms;
  for (int i = 0; i < NVEVS; i++) {
    cube[id] = vevs[i];
    ++id;
  }
  for (int i = 0; i < NVEVS; i++) {
    cube[id] = ui[i];
    ++id;
  }
  cube[id] = kf;
  ++id;
  cube[id] = kW;
  ++id;
  cube[id] = kZ;
  ++id;
  cube[id] = kh;
  ++id;
  cube[id] = lamfZ;
  ++id;
  cube[id] = kfZ;
  ++id;
  cube[id] = lamWZ;
  ++id;
}

void CSVMODEL::Print_PARAMETERS() {
  cout << "postp lognn ";
  for (int i = 0; i < NRandoms; i++) {
    cout << "th" << i << " ";
  }
  for (int i = 0; i < NVEVS; i++) {
    cout << "vev" << i << " ";
  }
  for (int i = 0; i < NVEVS; i++) {
    cout << "u" << i << " ";
  }
  cout << "kf kW kZ kh lamfZ kfZ lamWZ" << endl;
}
