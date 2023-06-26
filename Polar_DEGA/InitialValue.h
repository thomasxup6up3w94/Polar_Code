#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//#define N 8
//#define K 4
//#define n 3

#define N 1024
#define K 512
#define n 10

//#define N 2048
//#define K 1024
//#define n 11

//#define N 256
//#define K 128
//#define n 8



#define iteration     50				//BP iteration
#define Rm            2.0				//Es=Rm*Rate*Eb
#define Rate          (double)K/N
//#define loop          781250
//#define loop          24415
//#define loop          97657
//#define loop          58039
#define loop          1000000
#define Pi            3.14159265359
#define OneBySqrt2    1/sqrt(2)
//#define ChannelState  2.5

// Bernoulli-Gauss 把计砞﹚
#define P			0.01
#define Ratio_BG	100.0				//Ratio between average noise power in bad channel and in good channel
#define symbol_len  N / Rm

// AWAN 把计
#define  A_real     0.1            //A
#define  GIR_real   0.01           //tao
#define  m_length   4              // m!
#define  AWAN_LLR_max_count 4


// Markov-Gaussian 把计砞﹚
#define P_B			0.01
#define P_G			(1.0-P_B)
#define D_b			2.0                 //侥硈尿キА
#define Rho			(1.0-(1.0/(D_b*P_G)))	
#define P_GG		(double)(P_G+(Rho*P_B))
#define P_BB		(double)(P_B+(Rho*P_G))
#define P_GB		(1.0-P_GG)				
#define P_BG		(1.0-P_BB)	