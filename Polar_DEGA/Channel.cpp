#include "InitialValue.h"

void BPSK_AWGN(double *Cw, double *RCV, int len, float db)
{

	double u1, u2, r_max, rndm, s;// radom number variable
	float *noise;
	float scalar;
	int i = 0;

	noise = (float *)malloc(len * sizeof(float));
	r_max = RAND_MAX;          //for pc 

	for (i = 0; i < N; i++)
	{
		//scalar = (double)sqrt(2.0 * Rm * Rate * pow(10.0, db / 10.0));                   // = 1 / (sigma)
		//u1 = rand() / r_max;    //for pc
		//u2 = rand() / r_max;
		//noise[i] = (float)(sqrt(-2 * log(u1)) * cos(2 * Pi * u2)) / scalar;
		//i++;
		//if (i >= len)
		//	break;
		//noise[i] = (float)(sqrt(-2 * log(u1)) * sin(2 * Pi * u2)) / scalar;      
		////printf("(%f,%f)\n", noise[ i-1], noise[(i)]);

		do
		{
			rndm = rand() / (double)(RAND_MAX);
			u1 = rndm * 2 - 1.0;
			rndm = rand() / (double)(RAND_MAX);
			u2 = rndm * 2 - 1.0;
			s = u1 * u1 + u2 * u2;
		} while (s == 0.0 || s >1.0);

		scalar = (float)(1 / (2.0 * Rm * Rate * pow(10.0, db / 10.0)));                   // = (sigma) ^ 2
																						  //
																						  // (sigma) ^ 2 = N0/2;  SNR = Eb/N0;  N0 = Eb/SNR;  Eb = Es/(Rm*Rate)
																						  // (sigma) ^ 2 = Eb/2SNR
		 
																						// N0/2, the noise-variance
																						// R = log2M;     Es = R*Eb;    Es = m*Ec;    Es = m*R*Eb;                     
																						// SNR = Eb/No = Es/R*m*2*(sigma)^2 
		noise[i] = sqrt(scalar) * u1 * sqrt(-2.0*log(s) / s);
	}

	for (int j = 0; j < len; j++)
	{
		if (Cw[j] == 0)
			RCV[j] = noise[j] + (-1.0);			//bit 0 = -sqrt(Es)
		else
			RCV[j] = noise[j] + (1.0);		    //bit 1 = sqrt(Es)
	}
	free(noise);

}

void QPSK_AWGN_channel(double *Cw, double *RCV, int len, float db)
{

	double u1, u2, r_max, rndm, s;// radom number variable
	float *noise;
	float scalar;
	int i = 0;

	noise = (float *)malloc(len * sizeof(float));
	r_max = RAND_MAX;          //for pc 

	for (i = 0; i < len / 2; i++)
	{

		do
		{
			rndm = rand() / (double)(RAND_MAX);
			u1 = rndm * 2 - 1.0;
			rndm = rand() / (double)(RAND_MAX);
			u2 = rndm * 2 - 1.0;
			s = u1 * u1 + u2 * u2;
		} while (s == 0.0 || s >1.0);

		scalar = (float)(1 / (2.0 * Rm * Rate * pow(10.0, db / 10.0)));                   // = (sigma) ^ 2
																						  //
																						  // (sigma) ^ 2 = N0/2;  SNR = Eb/N0;  N0 = Eb/SNR;  Eb = Es/(Rm*Rate)
																						  // (sigma) ^ 2 = Eb/2SNR

		noise[2 * i] = sqrt(scalar) * u1 * sqrt(-2.0*log(s) / s);
		noise[(2 * i) + 1] = sqrt(scalar) * u2 * sqrt(-2.0*log(s) / s);

	}

	for (int j = 0; j < len; j++)
	{
		//printf("(noise[%d],noise[%d])    (%f,%f)\n",j,j+1, noise[j], noise[(j) + 1]);

		if (Cw[j] == 0)
			RCV[j] = noise[j] + (-OneBySqrt2);			//bit 0
		else
			RCV[j] = noise[j] + (OneBySqrt2);		    //bit 1

														//printf("(%f,%f)\n", RCV[2*i], RCV[2 * i+1]);
	}
	free(noise);

}

void QPSK_BernoulliGauss(double *Cw, double *RCV, int len, float db, double *BG)
{
	//FILE *error1;
	double u1, u2, r_max, rndm, s;// radom number variable
	float *noise;
	float scalar;
	float noisesum = 0;
	int i = 0;
	int T;

	noise = (float *)malloc(len * sizeof(float));
	r_max = RAND_MAX;          //for pc 
							   //error1 = fopen("BG.txt", "w");

	for (i = 0; i < len / 2; i++)
	{
		do
		{
			rndm = rand() / (double)(RAND_MAX);
			u1 = rndm * 2 - 1.0;
			rndm = rand() / (double)(RAND_MAX);
			u2 = rndm * 2 - 1.0;
			s = u1 * u1 + u2 * u2;
		} while (s == 0.0 || s >1.0);

		scalar = (float)(1 / (2.0 * Rm * Rate * pow(10.0, db / 10.0)));                     // = (sigma) ^ 2

		if (BG[2 * i] >= P)                                                                   																							// SNR=pow(10.0, snr(db) / 10.0);  db=R*Eb/No;   db=logSNR;  Ec=R*Eb;   Es=m*Ec;   Es=m*R*Eb;   Es:Energy per symbol;   Eb:Energy per bit;
																							// SNR = Eb/N0;  N0 = Eb/SNR;   No/2 = variance;																								
		{																					// (sigma) ^ 2 = N0/2;  SNR = Eb/N0;  N0 = Eb/SNR;  Eb = Es/(Rm*Rate)																						
			noise[2 * i] = sqrt(scalar) * u1 * sqrt(-2.0*log(s) / s);						// (sigma) ^ 2 = Eb/2SNR																							
			T = 0;																			//noise
		}																					//10*log(snr) = channelstate = db																			
		else if (BG[2 * i] < P)
		{
			noise[2 * i] = sqrt((Ratio_BG + 1) * scalar) * u1 * sqrt(-2.0*log(s) / s);
			T = 1;
		}


		if (BG[(2 * i) + 1] >= P)
		{
			noise[(2 * i) + 1] = sqrt(scalar) * u2 * sqrt(-2.0*log(s) / s);
			T = 0;
		}
		else if (BG[(2 * i) + 1] < P)
		{
			noise[(2 * i) + 1] = sqrt((Ratio_BG + 1) * scalar) * u2 * sqrt(-2.0*log(s) / s);
			T = 1;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (int j = 0; j < len; j++)
	//{
	//	if (noise[j] > 0) noisesum = noisesum + noise[j];
	//	else noisesum = noisesum - noise[j];
	//}
	//printf("%f", noisesum / len);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (int j = 0; j < len; j++)
	{
		//printf("(noise[%d],noise[%d])    (%f,%f)\n",j,j+1, noise[j], noise[(j) + 1]);
		//printf("noise[%d] %f\n", j, noise[j]);

		if (Cw[j] == 0)
			RCV[j] = noise[j] + (-1 * OneBySqrt2);			//bit 0
		else
			RCV[j] = noise[j] + ( 1 * OneBySqrt2);		    //bit 1

		//printf("(%f,%f)\n", RCV[2*i], RCV[2 * i+1]);
	}
	free(noise);
	//fclose(error1);
}

void QPSK_AWAN(double *Cw, double *RCV, int len, float db, double *c_prob, double *p_prob, double *sigmam, int *AWANstate)
{
	double u1, u2, r_max, rndm, s, R_V;// radom number variable
	float *noise;
	float scalar;
	float noisesum = 0;
	int i = 0;
	int T;

	noise = (float *)malloc(len * sizeof(float));
	r_max = RAND_MAX;          //for pc 
							   //error1 = fopen("BG.txt", "w");

	for (i = 0; i < len / 2; i++)
	{
		do
		{
			rndm = rand() / (double)(RAND_MAX);
			u1 = rndm * 2 - 1.0;
			rndm = rand() / (double)(RAND_MAX);
			u2 = rndm * 2 - 1.0;
			s = u1 * u1 + u2 * u2;
		} while (s == 0.0 || s >1.0);


		noise[2 * i] = sigmam[AWANstate[i * 2]] * u1 * sqrt(-2.0*log(s) / s);           //AWANstate[i * 2]=sigmam[m_pa]
		noise[(2 * i) + 1] = sigmam[AWANstate[i * 2 + 1]] * u2 * sqrt(-2.0*log(s) / s);
	
	}

	for (int j = 0; j < len; j++)
	{
		//printf("(noise[%d],noise[%d])    (%f,%f)\n",j,j+1, noise[j], noise[(j) + 1]);
		//printf("noise[%d] %f\n", j, noise[j]);

		if (Cw[j] == 0)
			RCV[j] = noise[j] + (-1 * OneBySqrt2);			//bit 0
		else
			RCV[j] = noise[j] + (1 * OneBySqrt2);		    //bit 1

															//printf("(%f,%f)\n", RCV[2*i], RCV[2 * i+1]);
	}
	free(noise);
}