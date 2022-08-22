// const int pbin = 29;
// const double p_mid[pbin] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 20, 23, 26, 30, 35, 40, 45, 50};
// static double p_lo[pbin] = {0};
// static double p_hi[pbin] = {0};

// const int etabin = 15;
// const double eta_mid[etabin] = {-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5};
// static double eta_lo[etabin] = {0};
// static double eta_hi[etabin] = {0};

// const int pbin = 1;
// const double p_mid[pbin] = {10};
// static double p_lo[pbin] = {9.5};
// static double p_hi[pbin] = {10.5};

// const int etabin = 1;
// const double eta_mid[etabin] = {0};
// static double eta_lo[etabin] = {-0.5};
// static double eta_hi[etabin] = {0.5};

const int pbin = 27;
const double p_lo[pbin] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 20, 23, 26};
static double p_hi[pbin] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 20, 23, 26, 30};

const int etabin = 14;
const double eta_lo[etabin] = {-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3};
static double eta_hi[etabin] = {-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5};

const int sysbins = 5;
const char* system_name[sysbins] = {"DIRC","pfRICH aerogel","dRICH aerogel","dRICH gas mid-point","barrel TOF"};

