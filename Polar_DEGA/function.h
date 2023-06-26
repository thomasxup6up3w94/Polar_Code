void Precomputation(int *, int *, int *);
void frozenlookup(int *, int *, double *, float , int , int *);

void PolarEncode(double *, double *, int *);
void SystematicPolarEncode(double *, double *, int *);

void BPSK_AWGN(double *, double *, int , float );
void QPSK_AWGN_channel(double *, double *, int , float );
void QPSK_BernoulliGauss(double *, double *, int , float , double *);
void QPSK_AWAN(double *, double *, int , float , double *, double *, double *, int *);

void PolarDecode(double *, double *, float , int *, int *, int *, int *, double *);
void SystematicPolarDecode(double *, double *, float , int *, int *, int *, int *, double *, double *, double *);
void SystematicPolarDecode_BG(double *, double *, float , int *, int *, int *, int *, double *, double );
void SystematicPolarDecode_AWAN_QPSK(double *, double *, float , int *, int *, int *, int *, double *, double *, int *, double );
void SystematicBPPolarDecode_BG(double *, double *, float , int *, int *, int *, int *, double *, double );
void SystematicBPPolarDecode_AWAN(double *, double *, float , int *, int *, int *, int *, double *, double *, int *, double *, double );


int factorial(int );