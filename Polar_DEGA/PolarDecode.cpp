#include "InitialValue.h"
#include <iostream>
#include <math.h>
double fabsmin(double x, double y)
{
	double z;
	if (fabs(x) > fabs(y))
		z = fabs(y);
	else
		z = fabs(x);

	return z;
}
double sgn(double x)
{

	if (x < 0)
		return -1;
	else if (x == 0)
		return 0;
	else
		return 1;
}
double atanh(double x)
{
	double y;
	y = (1.0 / 2.0)*log((1.0 + x) / (1.0 - x));
	return(y);
}
double logdomain_sum(double x, double y);
void PolarEncode(double *Message, double *Code, int *frozenbits);
void PolarDecode(double *RCV, double *RBit, float SNR, int *frozenbits, int *bitreversedindices, int *index_of_first1_from_MSB, int *index_of_first0_from_MSB, double *BG)
{
	double EbN0, nextlevel, lastlevel;
	int indx, st, ed, lev, j, i;
	int R = 0;
	double initialLLRs[N] = { 0 };
	double LLR[2 * N - 1] = { 0 };
	double d_hat[N] = { 0 };
	double info;
	int BITS[2][N - 1] = { 0 };

	int m_count, con_num;
	double sum, sqr_0, sqr_1;
	double tp[4];

	EbN0 = pow(10, (SNR / 10));
	double variance = (float)(1 / (2.0 * Rm * Rate * EbN0));       // = (sigma) ^ 2


	// BPSK AWGN receiver////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (i = 0; i < N; i++)
	//{
	//	initialLLRs[i] = -2 * (2 * Rm * Rate * EbN0) * RCV[i]; // QPSK要多乘以0.707(OneBySqrt2)
	//														   //
	//														   //                 Pr{y(i) | x(i) = 0} 
	//														   //  LLR(i) = log -----------------------
	//														   //                 Pr{y(i) | x(i) = 1}
	//}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// QPSK AWGN receiver ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },              //0
						  { OneBySqrt2,-OneBySqrt2 },              //1
						  { -OneBySqrt2, OneBySqrt2 },             //2
						  { -OneBySqrt2,-OneBySqrt2 } };	       //3

	for (i = 0; i < N / 2; i++)
	{
		for (con_num = 0; con_num<Rm*Rm; con_num++)
		{
			sum = 0, sqr_0 = 0, sqr_1 = 0;

			sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
			sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);

			tp[con_num] = 1/(2*Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance)));
		}

		initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
		initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	}

	//pdecode_LLRs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i = N; i <= 2 * N - 1; i++)
	{
		LLR[i - 1] = initialLLRs[i - N];
	}

	for (j = 1; j <= N; j++)
	{
		i = bitreversedindices[j - 1] + 1;

		//updateLLR///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == 1)
			nextlevel = n;
		else
		{
			lastlevel = index_of_first1_from_MSB[i - 1];

			st = pow(2, lastlevel - 1);
			ed = pow(2, lastlevel) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//lowerconv  g ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (BITS[0][indx - 1] == 0)
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] + LLR[ed + 2 * (indx - st) ];
				else
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] - LLR[ed + 2 * (indx - st) ];
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
			nextlevel = lastlevel - 1;
		}

		for (lev = nextlevel; lev >= 1; lev--)
		{
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//upperconv  f ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//LLR[indx - 1] = logdomain_sum(LLR[ed + 2 * (indx - st) ] + LLR[ed + 2 * (indx - st) + 1], 0)
				//	- logdomain_sum(LLR[ed + 2 * (indx - st) ], LLR[ed + 2 * (indx - st) + 1]);
				//upperconv  f ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				info = tanh(LLR[ed + 2 * (indx - st)] / 2) * tanh(LLR[ed + 2 * (indx - st) + 1] / 2);
				if (info == 1) info = 0.9999999999;
				if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
				LLR[indx - 1] = 2 * atanh(info);
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (frozenbits[i - 1] == -1)
		{
			if (LLR[0] > 0) d_hat[i - 1] = 0;
			else d_hat[i - 1] = 1;
		}
		else
			d_hat[i - 1] = frozenbits[i - 1];

		//updateBITS///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == N) continue;
		else if (i <= N / 2) BITS[0][0] = d_hat[i - 1];
		else
		{
			lastlevel = index_of_first0_from_MSB[i - 1];
			BITS[1][0] = d_hat[i - 1];

			for (lev = 1; lev <= lastlevel - 2; lev++)
			{
				st = pow(2, lev - 1);
				ed = pow(2, lev) - 1;
				for (indx = st; indx <= ed; indx++)
				{
					BITS[1][ed + 2 * (indx - st) ] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
					BITS[1][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
				}
			}

			lev = lastlevel - 1;
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				BITS[0][ed + 2 * (indx - st) ] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
				BITS[0][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = d_hat[i];
			R++;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

void SystematicPolarDecode(double *RCV, double *RBit, float SNR, int *frozenbits, int *bitreversedindices, int *index_of_first1_from_MSB, int *index_of_first0_from_MSB, double *BG, double *p_prob, double *sigmam)
{
	double EbN0, nextlevel, lastlevel;
	int indx, st, ed, lev, j, i;
	int R = 0;
	double initialLLRs[N] = { 0 };
	double LLR[2 * N - 1] = { 0 };
	double d_hat[N] = { 0 };
	double RBitCode[N] = { 0 };
	double info;
	int BITS[2][N - 1] = { 0 };

	int m_count, con_num;
	double sum, sqr_0, sqr_1;
	double tp[4];
	EbN0 = pow(10, (SNR / 10));
	double variance = (float)(1 / (2.0 * Rm * Rate * EbN0));       // = (sigma) ^ 2


	// BPSK AWGN receiver////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (i = 0; i < N; i++)
	//{
	//	initialLLRs[i] = -2 * (2 * Rm * Rate * EbN0) * RCV[i]; // QPSK要多乘以0.707(OneBySqrt2)
	//														   //
	//														   //                 Pr{y(i) | x(i) = 0} 
	//														   //  LLR(i) = log -----------------------
	//														   //                 Pr{y(i) | x(i) = 1}
	//}

	// QPSK AWGN receiver///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },              //0
						  { OneBySqrt2,-OneBySqrt2 },              //1
						  { -OneBySqrt2, OneBySqrt2 },             //2
						  { -OneBySqrt2,-OneBySqrt2 } };	       //3

	for (i = 0; i < N / 2; i++)
	{
		for (con_num = 0; con_num<Rm*Rm; con_num++)
		{
			sum = 0, sqr_0 = 0, sqr_1 = 0;

			sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
			sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);

			tp[con_num] = 1 / (2 * Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance)));
		}

		initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
		initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	}

	//pdecode_LLRs/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i = N; i <= 2 * N - 1; i++)
	{
		LLR[i - 1] = initialLLRs[i - N];
	}

	/*for (int P = 0; P < 2*N; P++)
	{
	printf("%d, ", bitreversedindices[P]);
	}*/

	for (j = 1; j <= N; j++)
	{
		i = bitreversedindices[j - 1] + 1;

		//updateLLR///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == 1)
			nextlevel = n;
		else
		{
			lastlevel = index_of_first1_from_MSB[i - 1];

			st = pow(2, lastlevel - 1);
			ed = pow(2, lastlevel) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//lowerconv  g ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (BITS[0][indx - 1] == 0)
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] + LLR[ed + 2 * (indx - st) ];
				else
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] - LLR[ed + 2 * (indx - st) ];
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
			}
			nextlevel = lastlevel - 1;
		}

		for (lev = nextlevel; lev >= 1; lev--)
		{
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//upperconv  f ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//LLR[indx-1] = logdomain_sum(LLR[ed + 2 * (indx - st) ] + LLR[ed + 2 * (indx - st) + 1], 0)
				//- logdomain_sum(LLR[ed + 2 * (indx - st) ], LLR[ed + 2 * (indx - st) + 1]);
				//soft upperconv  f ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				info = tanh(LLR[ed + 2 * (indx - st)] / 2) * tanh(LLR[ed + 2 * (indx - st) + 1] / 2);
				if (info == 1) info = 0.9999999999;
				if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
				LLR[indx - 1] = 2 * atanh(info);
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (frozenbits[i - 1] == -1)
		{
			if (LLR[0] > 0) d_hat[i - 1] = 0;
			else d_hat[i - 1] = 1;
		}
		else
			d_hat[i - 1] = frozenbits[i - 1];

		//updateBITS///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == N) continue;
		else if (i <= N / 2) BITS[0][0] = d_hat[i - 1];
		else
		{
			lastlevel = index_of_first0_from_MSB[i - 1];
			BITS[1][0] = d_hat[i - 1];

			for (lev = 1; lev <= lastlevel - 2; lev++)
			{
				st = pow(2, lev - 1);
				ed = pow(2, lev) - 1;
				for (indx = st; indx <= ed; indx++)
				{
					BITS[1][ed + 2 * (indx - st)] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
					BITS[1][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
				}
			}

			lev = lastlevel - 1;
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				BITS[0][ed + 2 * (indx - st) ] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
				BITS[0][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = d_hat[i];
			R++;
		}
	}

	PolarEncode(RBit, RBitCode, frozenbits);

	R = 0;
	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = RBitCode[i];
			R++;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

void SystematicPolarDecode_BG(double *RCV, double *RBit, float SNR, int *frozenbits, int *bitreversedindices, int *index_of_first1_from_MSB, int *index_of_first0_from_MSB, double *BG, double p_scan)
{
	double EbN0, nextlevel, lastlevel;
	int indx, st, ed, lev, j, i;
	int R = 0;
	double initialLLRs[N] = { 0 };
	double LLR[2 * N - 1] = { 0 };
	double d_hat[N] = { 0 };
	double RBitCode[N] = { 0 };
	double info;
	int BITS[2][N - 1] = { 0 };

	double pos_sum = 0, neg_sum = 0, sum, sqr_0, sqr_1;
	int m_count, con_num;
	double tp[4];

	EbN0 = pow(10, SNR / 10);
	double variance = (float)(1 / (2.0 * Rm * Rate * EbN0));       // = (sigma) ^ 2

	//For QPSK BG Rx known statistical features ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },              //0
	//					  { OneBySqrt2,-OneBySqrt2 },              //1
	//					  { -OneBySqrt2, OneBySqrt2 },             //2
	//					  { -OneBySqrt2,-OneBySqrt2 } };	       //3     //星座圖

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num<Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;

	//		sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
	//		sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);             //奇(偶)數位元距離平方

	//		tp[con_num] = (1-P)*( 1/(2*Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance))) ) + P*( 1/(2*Pi*(Ratio_BG+1)*variance) * (exp(-(sqr_0 + sqr_1) / (2*(Ratio_BG+1)*variance))) );
	//	}                                //白努力定理

	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));        //(2+3)/(0+1)機率
	//	//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//
	//
	//}

	//For QPSK BG Rx known impulse positions /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },              //0
	//					  { OneBySqrt2,-OneBySqrt2 },              //1
	//					  { -OneBySqrt2, OneBySqrt2 },             //2
	//					  { -OneBySqrt2,-OneBySqrt2 } };	       //3

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num<Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;

	//		sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
	//		sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);

	//		if (BG[2 * i] < P)
	//			tp[con_num] = 1 / (2 * Pi*(Ratio_BG + 1)*variance) * (exp(-(sqr_0 + sqr_1) / (2 * (Ratio_BG + 1)*variance))); //bg狀態   
	//		else
	//			tp[con_num] = 1 / (2 * Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance)));
	//	}

	//	initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));       //(2+3)/(0+1)機率
	//	initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//}

	//For QPSK BG  Rx unknown anything and do clip ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { {   OneBySqrt2, OneBySqrt2 },                //0
	//					  {   OneBySqrt2,-OneBySqrt2 },                //1
	//					  {  -OneBySqrt2, OneBySqrt2 },                //2
	//					  {  -OneBySqrt2,-OneBySqrt2 } };	           //3

	//double threshold_BG = sqrt(-2.0 * (variance)*log((p_scan / (1 - p_scan))*(sqrt(2 * Pi*variance))));         //BG裁剪臨界值掃描
	
	//double p_guess = 0.005;
	//double threshold_BG = sqrt(-2.0 * (variance)*log((p_guess / (1 - p_guess))*(sqrt(2 * Pi*variance))));   //BG裁剪臨界值
	
	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;

	//		double root_0 = RCV[2 * i] - sgnl[con_num][0];
	//		double root_1 = RCV[(2 * i) + 1] - sgnl[con_num][1];

	//		if (fabs(root_0) > threshold_BG) root_0 = sgn(root_0)*threshold_BG;    //奇數位元裁剪臨界值判斷
	//		if (fabs(root_1) > threshold_BG) root_1 = sgn(root_1)*threshold_BG;    //偶數位元裁剪臨界值判斷

	//		double root_0_root_1 = pow(root_0, 2) + pow(root_1, 2);

	//		tp[con_num] = 1 / (2 * Pi*variance) *(exp(-(root_0_root_1) / (2 * (variance))));
	//	}

	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//}



	//pdecode_LLRs/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i = N; i <= 2 * N - 1; i++)
	{
		LLR[i - 1] = initialLLRs[i - N];
	}

	/*for (int P = 0; P < 2*N; P++)
	{
	printf("%d, ", bitreversedindices[P]);
	}*/

	for (j = 1; j <= N; j++)
	{
		i = bitreversedindices[j - 1] + 1;

		//updateLLR///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == 1)
			nextlevel = n;
		else
		{
			lastlevel = index_of_first1_from_MSB[i - 1];

			st = pow(2, lastlevel - 1);
			ed = pow(2, lastlevel) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//lowerconv  g ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (BITS[0][indx - 1] == 0)
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] + LLR[ed + 2 * (indx - st) ];
				else
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] - LLR[ed + 2 * (indx - st) ];
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
			nextlevel = lastlevel - 1;
		}

		for (lev = nextlevel; lev >= 1; lev--)
		{
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//upperconv  f ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//LLR[indx-1] = logdomain_sum(LLR[ed + 2 * (indx - st) ] + LLR[ed + 2 * (indx - st) + 1], 0)
				//- logdomain_sum(LLR[ed + 2 * (indx - st) ], LLR[ed + 2 * (indx - st) + 1]);

				//soft upperconv  f ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				info = tanh(LLR[ed + 2 * (indx - st)] / 2) * tanh(LLR[ed + 2 * (indx - st) + 1] / 2);
				if (info == 1) info = 0.9999999999;
				if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
				LLR[indx - 1] = 2 * atanh(info);
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (frozenbits[i - 1] == -1)
		{
			if (LLR[0] > 0) d_hat[i - 1] = 0;
			else d_hat[i - 1] = 1;
		}
		else
			d_hat[i - 1] = frozenbits[i - 1];

		//updateBITS///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == N) continue;
		else if (i <= N / 2) BITS[0][0] = d_hat[i - 1];
		else
		{
			lastlevel = index_of_first0_from_MSB[i - 1];
			BITS[1][0] = d_hat[i - 1];

			for (lev = 1; lev <= lastlevel - 2; lev++)
			{
				st = pow(2, lev - 1);
				ed = pow(2, lev) - 1;
				for (indx = st; indx <= ed; indx++)
				{
					BITS[1][ed + 2 * (indx - st) ] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
					BITS[1][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
				}
			}

			lev = lastlevel - 1;
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				BITS[0][ed + 2 * (indx - st) ] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
				BITS[0][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = d_hat[i];
			R++;
		}
	}

	PolarEncode(RBit, RBitCode, frozenbits);

	R = 0;
	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = RBitCode[i];
			R++;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

void SystematicPolarDecode_AWAN_QPSK(double *RCV, double *RBit, float SNR, int *frozenbits, int *bitreversedindices, int *index_of_first1_from_MSB, int *index_of_first0_from_MSB, double *am_prob, double *sigmam, int *AWANstate, double A_real_scan)
{
	double EbN0, nextlevel, lastlevel;
	int indx, st, ed, lev, j, i;
	int R = 0;
	double initialLLRs[N] = { 0 };
	double LLR[2 * N - 1] = { 0 };
	double d_hat[N] = { 0 };
	double info;
	double RBitCode[N] = { 0 };
	int BITS[2][N - 1] = { 0 };

	double pos_sum = 0, neg_sum = 0, sum, sqr_0, sqr_1;
	int m_count, con_num;
	double tp[4];

	EbN0 = pow(10, (SNR / 10));
	double variance = (float)(1 / (2.0 * Rm * Rate * EbN0));       // = (sigma) ^ 2


	//AWAN Rx known statistical features////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2  },              //0
	//					  { OneBySqrt2,-OneBySqrt2  },              //1
	//					  { -OneBySqrt2, OneBySqrt2 },              //2
	//					  { -OneBySqrt2,-OneBySqrt2 } };	        //3

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num<Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;

	//		sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
	//		sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);

	//		for (m_count = 0; m_count <= AWAN_LLR_max_count; m_count++)     // m_coun = AWAN_LLR_max_count sigma 累加 最大值
	//		{
	//			sum +=  exp(-A_real)*(am_prob[m_count]) *(1/(2*Pi*(sigmam[AWANstate[i * 2]]*sigmam[AWANstate[i * 2]]))) * (exp(-(sqr_0 + sqr_1) / (2 * pow(sigmam[AWANstate[i * 2]], 2))));
	//		}

	//		tp[con_num] = sum;
	//	}

	//	initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//}

	//AWAN Rx known impulse positions//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },               //0
						  { OneBySqrt2,-OneBySqrt2 },				//1
						  {-OneBySqrt2, OneBySqrt2 },				//2
						  {-OneBySqrt2,-OneBySqrt2 } };				//3

	for (i = 0; i < N / 2; i++)
	{
		for (con_num = 0; con_num<Rm*Rm; con_num++)
		{
			sum = 0, sqr_0 = 0, sqr_1 = 0;

			sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
			sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);

			tp[con_num] = 1 / ((2 * Pi)*sigmam[AWANstate[i * 2]] * sigmam[AWANstate[i * 2]]) * (exp(-(sqr_0 + sqr_1) / (2 * pow(sigmam[AWANstate[i * 2]], 2))));
		}

		initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
		initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	}

	//AWAN Rx unknown any information and do clip ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2   },              //0
	//					  {   OneBySqrt2,-OneBySqrt2 },              //1
	//					  {  -OneBySqrt2, OneBySqrt2 },              //2
	//					  {  -OneBySqrt2,-OneBySqrt2 } };	         //3

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;
	//		double root_0 = RCV[2 * i] - sgnl[con_num][0];
	//		double root_1 = RCV[(2 * i) + 1] - sgnl[con_num][1];

	//		double threshold_0 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_scan)) / (exp(-A_real_scan))))));      //奇數位元裁剪臨界值掃描
	//		double threshold_1 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_scan)) / (exp(-A_real_scan))))));      //偶數位元裁剪臨界值掃描

	//		//double A_real_guess = 0.005;
	//		//double threshold_0 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_guess)) / (exp(-A_real_guess))))));        //裁剪臨界值判斷
	//		//double threshold_1 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_guess)) / (exp(-A_real_guess))))));


	//		if (fabs(root_0) > threshold_0)  root_0 = sgn(root_0)*threshold_0;
	//		if (fabs(root_1) > threshold_1)	 root_1 = sgn(root_1)*threshold_1;

	//		double root_0_root_1 = pow(root_0, 2) + pow(root_1, 2);

	//		tp[con_num] = 1 / (2 * Pi*(variance))*(exp(-(root_0_root_1) / (2 * (variance))));
	//	}


	//	initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));

	//}


	//pdecode_LLRs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (i = N; i <= 2 * N - 1; i++)
	{
		LLR[i - 1] = initialLLRs[i - N];
	}

	for (j = 1; j <= N; j++)
	{
		i = bitreversedindices[j - 1] + 1;

		//updateLLR///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == 1)
			nextlevel = n;
		else
		{
			lastlevel = index_of_first1_from_MSB[i - 1];

			st = pow(2, lastlevel - 1);
			ed = pow(2, lastlevel) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//lowerconv  g ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (BITS[0][indx - 1] == 0)
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] + LLR[ed + 2 * (indx - st) ];
				else
					LLR[indx - 1] = LLR[ed + 2 * (indx - st) + 1] - LLR[ed + 2 * (indx - st) ];
				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
			nextlevel = lastlevel - 1;
		}

		for (lev = nextlevel; lev >= 1; lev--)
		{
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				//upperconv  f ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//LLR[indx - 1] = logdomain_sum(LLR[ed + 2 * (indx - st) ] + LLR[ed + 2 * (indx - st) + 1], 0)
				//	- logdomain_sum(LLR[ed + 2 * (indx - st) ], LLR[ed + 2 * (indx - st) + 1]);
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//soft upperconv  f ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
				info = tanh(LLR[ed + 2 * (indx - st)] / 2) * tanh(LLR[ed + 2 * (indx - st) + 1] / 2);
				if (info == 1) info = 0.9999999999;
				if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
				LLR[indx - 1] = 2 * atanh(info);
				///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (frozenbits[i - 1] == -1)
		{
			if (LLR[0] > 0) d_hat[i - 1] = 0;
			else d_hat[i - 1] = 1;
		}
		else
			d_hat[i - 1] = frozenbits[i - 1];

		//updateBITS///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (i == N) continue;
		else if (i <= N / 2) BITS[0][0] = d_hat[i - 1];
		else
		{
			lastlevel = index_of_first0_from_MSB[i - 1];
			BITS[1][0] = d_hat[i - 1];

			for (lev = 1; lev <= lastlevel - 2; lev++)
			{
				st = pow(2, lev - 1);
				ed = pow(2, lev) - 1;
				for (indx = st; indx <= ed; indx++)
				{
					BITS[1][ed + 2 * (indx - st) ] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
					BITS[1][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
				}
			}

			lev = lastlevel - 1;
			st = pow(2, lev - 1);
			ed = pow(2, lev) - 1;
			for (indx = st; indx <= ed; indx++)
			{
				BITS[0][ed + 2 * (indx - st) ] = (BITS[0][indx - 1] + BITS[1][indx - 1]) % 2;
				BITS[0][ed + 2 * (indx - st) + 1] = BITS[1][indx - 1];
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	for (i = 0; i < N; i++)
	{
		if (frozenbits[i] == -1)
		{
			RBit[R] = d_hat[i];
			R++;
		}
	}

	PolarEncode(RBit, RBitCode, frozenbits);

	R = 0;
	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = RBitCode[i];
			R++;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}


void SystematicBPPolarDecode_BG(double *RCV, double *RBit, float SNR, int *frozenbits, int *bitreversedindices, int *index_of_first1_from_MSB, int *index_of_first0_from_MSB, double *BG, double p_scan)
{
	double EbN0, nextlevel, lastlevel;
	int j, i, base, nB, B, s;
	int R = 0;
	double RBitCode[N] = { 0 };
	double initialLLRs[N] = { 0 };
	double d_hat[N] = { 0 };
	double info;

	double L0[N][n + 1] = { 1 };     //LLR(LR)值
	double R0[N][n + 1] = { 1 };     //R值

	int R1[N][n + 1] = { 1 };        //G-matrix檢查矩陣
	int R2[N] = { 1 };               //做x估測
	int L1[N] = { 1 };               //做u估測
	int error = { 0 };               //判斷x=Gu


	double pos_sum = 0, neg_sum = 0, sum, sqr_0, sqr_1;
	int m_count, con_num;
	double tp[4];
	
	EbN0 = pow(10.0, SNR / 10.0);
	double variance = (float)(1 / (2.0 * Rm * Rate * EbN0));       // = (sigma) ^ 2



	//BPSK AWGN receiver //////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (i = 0; i < N; i++)
	//{
	//	
	//	initialLLRs[i] = -2.0 * 2.0 * Rm * Rate * EbN0 * RCV[i] ;
	//	//initialLLRs[i] = exp(-2.0 * 2.0 * Rm * Rate * EbN0 * RCV[i]);
	//															////////		 Pr{y(i) | x(i) = 0} 
	//													////////LLR(i) = log -----------------------
	//															////////		 Pr{y(i) | x(i) = 1}
	//}
	
	// QPSK AWGN receiver ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2   },              //0
	//					  {   OneBySqrt2,-OneBySqrt2 },              //1
	//					  {  -OneBySqrt2, OneBySqrt2 },              //2
	//					  {  -OneBySqrt2,-OneBySqrt2 } };	         //3
	//
	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;
	//
	//		sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
	//		sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);
	//
	//		tp[con_num] = 1 / (2 * Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance)));
	//	}
	//
	//	initialLLRs[2 * i] = (tp[2] + tp[3]) / (tp[0] + tp[1]);
	//	initialLLRs[2 * i + 1] = (tp[1] + tp[3]) / (tp[0] + tp[2]);

	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//}

	//For QPSK BG Rx known statistical features////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },              //0
						  { OneBySqrt2,-OneBySqrt2 },              //1
						  { -OneBySqrt2, OneBySqrt2 },             //2
						  { -OneBySqrt2,-OneBySqrt2 } };	       //3     //星座圖

	for (i = 0; i < N / 2; i++)
	{
		for (con_num = 0; con_num<Rm*Rm; con_num++)
		{
			sum = 0, sqr_0 = 0, sqr_1 = 0;

			sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
			sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);             //奇(偶)位元距離

			tp[con_num] = (1-P)*( 1/(2*Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance))) ) + P*( 1/(2*Pi*(Ratio_BG+1)*variance) * (exp(-(sqr_0 + sqr_1) / (2*(Ratio_BG+1)*variance))) );
		}                                //白努力定理

		//initialLLRs[2 * i] = (tp[2] + tp[3]) / (tp[0] + tp[1]);        //(2+3)/(0+1)機率
		//initialLLRs[2 * i + 1] = (tp[1] + tp[3]) / (tp[0] + tp[2]);
	
		initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));        //(2+3)/(0+1)機率
		initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	
	}

	//For QPSK BG Rx known impulse positions/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },              //0
	//					  { OneBySqrt2,-OneBySqrt2 },              //1
	//					  { -OneBySqrt2, OneBySqrt2 },             //2
	//					  { -OneBySqrt2,-OneBySqrt2 } };	       //3

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;

	//		sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
	//		sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);

	//		if (BG[2 * i] < P)
	//			tp[con_num] = 1 / (2 * Pi*(Ratio_BG + 1)*variance) * (exp(-(sqr_0 + sqr_1) / (2 * (Ratio_BG + 1)*variance))); //bg狀態   
	//		else
	//			tp[con_num] = 1 / (2 * Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance)));
	//	}

	//	initialLLRs[2 * i] = (tp[2] + tp[3]) / (tp[0] + tp[1]);
	//	initialLLRs[2 * i + 1] = (tp[1] + tp[3]) / (tp[0] + tp[2]);
	//
	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));

	//}

	//For QPSK BG  Rx unknown anything and do clip ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { {   OneBySqrt2, OneBySqrt2 },                //0
	//					  {   OneBySqrt2,-OneBySqrt2 },                //1
	//					  {  -OneBySqrt2, OneBySqrt2 },                //2
	//					  {  -OneBySqrt2,-OneBySqrt2 } };	           //3

	////double threshold_BG = sqrt(-2.0 * (variance)*log((p_scan / (1 - p_scan))*(sqrt(2 * Pi*variance))));   //BG裁剪臨界值掃描
	//
	//double p_guess = 0.005;
	//double threshold_BG = sqrt(-2.0 * (variance)*log((p_guess / (1 - p_guess))*(sqrt(2 * Pi*variance))));   //BG裁剪臨界值
	//
	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;

	//		double root_0 = RCV[2 * i] - sgnl[con_num][0];
	//		double root_1 = RCV[(2 * i) + 1] - sgnl[con_num][1];

	//		if (fabs(root_0) > threshold_BG) root_0 = sgn(root_0)*threshold_BG;
	//		if (fabs(root_1) > threshold_BG) root_1 = sgn(root_1)*threshold_BG;

	//		double root_0_root_1 = pow(root_0, 2) + pow(root_1, 2);

	//		tp[con_num] = 1 / (2 * Pi*variance) *(exp(-(root_0_root_1) / (2 * (variance))));
	//	}

	//	initialLLRs[2 * i] = (tp[2] + tp[3]) / (tp[0] + tp[1]);
	//	initialLLRs[2 * i + 1] = (tp[1] + tp[3]) / (tp[0] + tp[2]);

	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//}



	//pdecode_LRs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (j = 1; j < n; j++)                  //L和R初值(無log)
	//{
	//	for (i = 0; i < N; i++)
	//	{
	//		L0[i][j] = 1.0;
	//		R0[i][j] = 1.0;
	//
	//	}
	//}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	for (int i = 0; i < N; i++)
	{
		L0[i][n] = initialLLRs[i];      //LLR
		if (frozenbits[i] == 0)

			//R0[i][0] = INT_MAX;
			R0[i][0] = 1000000000;
		else if (frozenbits[i] == -1)

			R0[i][0] = 0;               //message 機率比值取log 
		    //R0[i][0] = 1.0;
	}



	for (s = 0; s < iteration; s++)                //BP疊代次數
	{
		for (int i = 1; i <= n; i++)                //編碼端 R 訊息做疊代
		{
			B = pow(2, (n - i + 1));
			nB = pow(2, (i - 1));
			for (int j = 1; j <= nB; j++)
			{
				base = (j - 1) * B;
				for (int l = 1; l <= B / 2; l++)
				{
					
					//R0[base + l - 1][i] = (1 + R0[base + l - 1][i - 1] * L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]) / (R0[base + l - 1][i - 1] + L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]);
					//R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] * ((1 + R0[base + l - 1][i - 1] * L0[base + l - 1][i]) / (R0[base + l - 1][i - 1] + L0[base + l - 1][i]));


					//取 log /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					R0[base + l - 1][i] = 0.9375*sgn(R0[base + l - 1][i - 1])*sgn(L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1])*fabsmin(R0[base + l - 1][i - 1], L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] + 0.9375*sgn(R0[base + l - 1][i - 1])*sgn(L0[base + l - 1][i])*fabsmin(R0[base + l - 1][i - 1], L0[base + l - 1][i]);
					
					

					//info = tanh(R0[base + l - 1][i - 1] / 2) * tanh((L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]) / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//R0[base + l - 1][i] = 2 * atanh(info);
					//info = tanh(R0[base + l - 1][i - 1] / 2) * tanh(L0[base + l - 1][i] / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] + 2 * atanh(info);


					//R0[base + l - 1][i] = logdomain_sum(R0[base + l - 1][i - 1] + L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1], 0)
					//                       - logdomain_sum(R0[base + l - 1][i - 1], L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					//R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] + logdomain_sum(R0[base + l - 1][i - 1] + L0[base + l - 1][i], 0)
					//                       - logdomain_sum(R0[base + l - 1][i - 1] , L0[base + l - 1][i]);

					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


				}
			}

		}

		for (int i = n; i >= 1; i--)                   //解碼端 L 訊息做疊代
		{

			B = pow(2, (n - i + 1));
			nB = pow(2, (i - 1));
			for (int j = 1; j <= nB; j++)
			{
				base = (j - 1) * B;
				for (int l = 1; l <= B / 2; l++)
				{

					//L0[base + l - 1][i - 1] = (1 + L0[base + l - 1][i] * L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]) / (L0[base + l - 1][i] + L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]);         
					//L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] * ((1 + L0[base + l - 1][i] * R0[base + l - 1][i - 1]) / (L0[base + l - 1][i] + R0[base + l - 1][i - 1]));


					//取 log /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


					L0[base + l - 1][i - 1] = 0.9375*sgn(L0[base + l - 1][i])*sgn(L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1])*fabsmin(L0[base + l - 1][i], L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] + 0.9375*sgn(L0[base + l - 1][i])*sgn(R0[base + l - 1][i - 1])*fabsmin(L0[base + l - 1][i], R0[base + l - 1][i - 1]);


					//info = tanh(L0[base + l - 1][i] / 2) * tanh((L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]) / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//L0[base + l - 1][i - 1] = 2 * atanh(info);
					//info = tanh(L0[base + l - 1][i] / 2) * tanh(R0[base + l - 1][i - 1] / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] + 2 * atanh(info);


					//L0[base + l - 1][i - 1] = logdomain_sum(L0[base + l - 1][i] + L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1], 0)
					//							- logdomain_sum(L0[base + l - 1][i] , L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					//L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] + logdomain_sum(L0[base + l - 1][i] + R0[base + l - 1][i - 1], 0)
					//							- logdomain_sum(L0[base + l - 1][i] , R0[base + l - 1][i - 1]);


				}
			}
		}
		

		//G-matrix檢查矩陣////////////////////////////////////////////////////////////////////////////////////////////////////
		for (i = 0; i < N; i++)             
		{

			//if (R0[i][n]>= 1.0) L1[i] = 0;
			if (R0[i][n]>= 0) L1[i] = 0;                      //R0取log
			else L1[i] = 1;

		}
		for (i = 0; i < N; i++)
		{
			if (frozenbits[i] == -1)
			{
				//if (L0[i][0] >= 1.0) R1[i][0] = 0;         
				if (L0[i][0]>= 0) R1[i][0] = 0;              //L0取log
				else R1[i][0] = 1;
			}
			else
				R1[i][0] = frozenbits[i];
		}

		for (int i = 1; i <= n; i++)                
		{
			B = pow(2, (n - i + 1));
			nB = pow(2, (i - 1));
			for (int j = 1; j <= nB; j++)
			{
				base = (j - 1) * B;
				for (int l = 1; l <= B / 2; l++)
				{
					R1[base + l - 1][i] = (R1[base + l - 1][i - 1] + R1[base + B / 2 + l - 1][i - 1]) % 2;
					R1[base + B / 2 + l - 1][i] = R1[base + B / 2 + l - 1][i - 1] % 2;
				}
			}
		}
		for (int i = 0; i < N; i++)
		{
			R2[i] = R1[i][n];
		}



		for (i = 0; i < N; i++)
		{

			if (L1[i] != R2[i])
				error++;


		}

		if (error == 0 || s >= iteration - 1)
			break;


	}
	
	
	
	for (i = 0; i < N; i++)                                    //軟式決策
	{
		if (frozenbits[i] == -1)
		{
			//if (L0[i][0] >= 1.0) RBit[R] = 0;
			if (L0[i][0]  >= 0) RBit[R] = 0;                   //L0取log
			else RBit[R] = 1;
			R++;
		}

	}

	//BP 另一種寫法/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (int BPloop = 1; BPloop <= iteraion; BPloop++)
	//{
	//	for (j = n - 1; j >= 0; j--)
	//	{
	//		for (i = 0; i < N; i++)
	//		{
	//			/*
	//			I
	//			0--+--.--+-----.
	//			|     |
	//			1--.--.-----+--.
	//			|  |
	//			2--+--.--.-----.
	//			|        |
	//			3--.--.-----.--.
	//			L 0     1        2        第L個Level是每2^(L-1)個Bits就換一種運算(%2是為了因應當L=0時的狀況)
	//			*/	
	//			int r = pow(2, j);
	//			if ((i / r) % 2 == 0) R1[i][j] = (R1[i][j + 1] * R1[i + r][j + 1] * LR1[i + r][j] + 1) / (R1[i][j + 1] + R1[i + r][j + 1] * LR1[i + r][j]);
	//			else				  R1[i][j] = R1[i][j + 1] * (R1[i - r][j + 1] * LR1[i - r][j] + 1) / (R1[i - r][j + 1] + LR1[i - r][j]);
	//		}
	//	}
	//	for (j = 0; j <= n + 1; j++)
	//	{
	//		for (i = 0; i < N; i++)
	//		{
	//			/*
	//			I
	//			0--+--.--+-----.
	//			|     |
	//			1--.--.-----+--.
	//			|  |
	//			2--+--.--.-----.
	//			|        |
	//			3--.--.-----.--.
	//			L 0     1        2        第L個Level是每2^(L-1)個Bits就換一種運算(%2是為了因應當L=0時的狀況)
	//			*/
	//			int r = pow(2, j);
	//			if ((i / r) % 2 == 0) LR1[i][j] = (LR1[i][j - 1] * LR1[i + r][j - 1] * R1[i + r][j] + 1) / (LR1[i][j - 1] + LR1[i + r][j - 1] * R1[i + r][j]);
	//			else				  LR1[i][j] = LR1[i][j - 1] * (LR1[i - r][j - 1] * R1[i - r][j] + 1) / (LR1[i - r][j - 1] + R1[i - r][j]);
	//		}
	//	}
	//}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (i = 0; i < N; i++)
	//{
	//	if (frozenbits[i] == -1)
	//	{
	//      //if (L0[i][0]>= 0) d_hat[i] = 0;                 //取log
	//		if (L0[i][0]>= 1.0) d_hat[i] = 0;
	//		else d_hat[i] = 1;
	//	}
	//	else
	//	d_hat[i] = frozenbits[i];
	//}
	//for (i = 0; i < N; i++)
	//{
	//	//printf("%f\t", d_hat[i]);
	//	if (frozenbits[i] == -1)
	//	{
	//		RBit[R] = d_hat[i];
	//		R++;
	//	}
	//}
	////////////////////////////////////////////////////////////////////////////////////

	//再encode一次///////////////////////////////////////////////////////////////////////////////////
	PolarEncode(RBit, RBitCode, frozenbits);

	R = 0;
	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = RBitCode[i];
			R++;
		}
	}


}///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SystematicBPPolarDecode_AWAN(double *RCV, double *RBit, float SNR, int *frozenbits, int *bitreversedindices, int *index_of_first1_from_MSB, int *index_of_first0_from_MSB, double *am_prob, double *sigmam, int *AWANstate, double *c_prob, double A_real_scan)

{
	double EbN0, nextlevel, lastlevel;
	int j, i, base, nB, B, s;
	int R = 0;
	double RBitCode[N] = { 0 };
	double initialLLRs[N] = { 0 };
	double d_hat[N] = { 0 };
	double info;

	double L0[N][n + 1] = { 1 };     //LLR(LR)值
	double R0[N][n + 1] = { 1 };     //R值

	int R1[N][n + 1] = { 1 };        //G-matrix檢查矩陣
	int R2[N] = { 1 };               //做x估測
	int L1[N] = { 1 };               //做u估測
	int error = { 0 };               //判斷x=Gu

	double pos_sum = 0, neg_sum = 0, sum, sqr_0, sqr_1;
	int m_count, con_num;
	double tp[4];

	EbN0 = pow(10.0, SNR / 10.0);
	double variance = (float)(1 / (2.0 * Rm * Rate * EbN0));       // = (sigma) ^ 2


	//AWAN Rx unknown any information and do clip ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2   },              //0
	//					  {   OneBySqrt2,-OneBySqrt2 },              //1
	//					  {  -OneBySqrt2, OneBySqrt2 },              //2
	//					  {  -OneBySqrt2,-OneBySqrt2 } };	         //3

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;
	//		double root_0 = RCV[2 * i] - sgnl[con_num][0];
	//		double root_1 = RCV[(2 * i) + 1] - sgnl[con_num][1];

	//		//double threshold_0 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_scan)) / (exp(-A_real_scan))))));        //裁剪臨界值掃描
	//		//double threshold_1 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_scan)) / (exp(-A_real_scan))))));

	//		double A_real_guess = 0.005;
	//		double threshold_0 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_guess)) / (exp(-A_real_guess))))));        //裁剪臨界值判斷
	//		double threshold_1 = sqrt(-2.0 * (sigmam[0] * sigmam[0])*log((sqrt(2 * Pi * sigmam[0] * sigmam[0])*((1 - exp(-A_real_guess)) / (exp(-A_real_guess))))));



	//		if (fabs(root_0) > threshold_0)  root_0 = sgn(root_0)*threshold_0;			
	//		if (fabs(root_1) > threshold_1)	 root_1 = sgn(root_1)*threshold_1;
	//			
	//		double root_0_root_1 = pow(root_0,2) + pow(root_1,2);

	//		tp[con_num] = 1 / (2 * Pi*(variance))*(exp(-(root_0_root_1) / (2 * (variance))));
	//	}

	//	initialLLRs[2 * i] = (((tp[2] + tp[3]) / (tp[0] + tp[1])));
	//	initialLLRs[2 * i + 1] = (((tp[1] + tp[3]) / (tp[0] + tp[2])));
	//
	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
    //  //initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//}

	//BPSK AWGN receiver ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//for (i = 0; i < N; i++)
	//{
	//	initialLLRs[i] = exp(-2.0 * 2.0 * Rm * Rate * EbN0 * RCV[i]);
	//	//initialLLRs[i] = -2.0 * 2.0 * Rm * Rate * EbN0 * RCV[i];
	//																				//
	//																				//                 Pr{y(i) | x(i) = 0} 
	//																				//  LLR(i) = log -----------------------
	//																				//                 Pr{y(i) | x(i) = 1}
	//}

	//AWAN Rx known statistical features /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { {  OneBySqrt2, OneBySqrt2 },             //0
	//					  {  OneBySqrt2,-OneBySqrt2 },             //1
	//					  { -OneBySqrt2, OneBySqrt2 },             //2
	//					  { -OneBySqrt2,-OneBySqrt2 } };		   //3
	//						

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;
	//		sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
	//		sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);
	//		for (m_count = 0; m_count <= AWAN_LLR_max_count; m_count++)     // m_coun = AWAN_LLR_max_count sigma 累加 最大值
	//		{
	//			
	//			sum += exp(-A_real)*(am_prob[m_count]) *(1 / (2 * Pi*(sigmam[AWANstate[i * 2]] * sigmam[AWANstate[i * 2]]))) * (exp(-(sqr_0 + sqr_1) / (2 * pow(sigmam[AWANstate[i * 2]], 2))));

	//		}
	//		tp[con_num] = sum;
	//	}
	//	
	//	initialLLRs[2 * i] = ((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	initialLLRs[2 * i + 1] = ((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//
	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//
	//}

	//AWAN Rx known impulse positions //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },               //0
						{   OneBySqrt2,-OneBySqrt2 },               //1
						{  -OneBySqrt2, OneBySqrt2 },               //2
						{  -OneBySqrt2,-OneBySqrt2 } };	            //3
						

	for (i = 0; i < N / 2; i++)
	{
		for (con_num = 0; con_num < Rm*Rm; con_num++)
		{
			sum = 0, sqr_0 = 0, sqr_1 = 0;
			sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
			sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);
			tp[con_num] = 1 / ((2 * Pi)*sigmam[AWANstate[i * 2]] * sigmam[AWANstate[i * 2]]) * (exp(-(sqr_0 + sqr_1) / (2 * pow(sigmam[AWANstate[i * 2]], 2))));
		
		}
		initialLLRs[2 * i] = ((tp[2] + tp[3]) / (tp[0] + tp[1]));
		initialLLRs[2 * i + 1] = ((tp[1] + tp[3]) / (tp[0] + tp[2]));

		//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
		//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	}

	// QPSK AWGN receiver ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//double sgnl[4][2] = { { OneBySqrt2, OneBySqrt2 },               //0
	//					  {   OneBySqrt2,-OneBySqrt2 },             //1
	//					  {  -OneBySqrt2, OneBySqrt2 },             //2
	//					  {  -OneBySqrt2,-OneBySqrt2 } };	        //3

	//for (i = 0; i < N / 2; i++)
	//{
	//	for (con_num = 0; con_num < Rm*Rm; con_num++)
	//	{
	//		sum = 0, sqr_0 = 0, sqr_1 = 0;
	//		sqr_0 = pow(RCV[2 * i] - sgnl[con_num][0], 2);
	//		sqr_1 = pow(RCV[(2 * i) + 1] - sgnl[con_num][1], 2);
	//		tp[con_num] = 1 / (2 * Pi*variance) * (exp(-(sqr_0 + sqr_1) / (2 * variance)));
	//	}
	//	initialLLRs[2 * i] = ((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	initialLLRs[2 * i + 1] = ((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//
	//	//initialLLRs[2 * i] = log((tp[2] + tp[3]) / (tp[0] + tp[1]));
	//	//initialLLRs[2 * i + 1] = log((tp[1] + tp[3]) / (tp[0] + tp[2]));
	//
	//}



	//pdecode_LRs////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//for (j = 1; j < n; j++)
	//{
	//	for (i = 0; i < N; i++)
	//	{
	//		L0[i][j] = 1.0;
	//		R0[i][j] = 1.0;
	//
	//	}
	//}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	for (int i = 0; i < N; i++)
	{
		L0[i][n] = initialLLRs[i];          //LLR
		if (frozenbits[i] == 0)

			//R0[i][0] = INT_MAX;
			R0[i][0] = 1000000000;
		else if (frozenbits[i] == -1)

			R0[i][0] = 0;                     //R0取log
			//R0[i][0] = 1.0;
	}



	for (s = 0; s < iteration; s++)
	{
		for (int i = 1; i <= n; i++)
		{
			B = pow(2, (n - i + 1));
			nB = pow(2, (i - 1));
			for (int j = 1; j <= nB; j++)
			{
				base = (j - 1) * B;
				for (int l = 1; l <= B / 2; l++)
				{
					//R0[base + l - 1][i] = (1 + R0[base + l - 1][i - 1] * L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]) / (R0[base + l - 1][i - 1] + L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]);
					//R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] * ((1 + R0[base + l - 1][i - 1] * L0[base + l - 1][i]) / (R0[base + l - 1][i - 1] + L0[base + l - 1][i]));

					//R0 取 log ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					R0[base + l - 1][i] = 0.9375*sgn(R0[base + l - 1][i - 1])*sgn(L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1])*fabsmin(R0[base + l - 1][i - 1], L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] + 0.9375*sgn(R0[base + l - 1][i - 1])*sgn(L0[base + l - 1][i])*fabsmin(R0[base + l - 1][i - 1], L0[base + l - 1][i]);


					//info = tanh(R0[base + l - 1][i - 1] / 2) * tanh((L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]) / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//R0[base + l - 1][i] = 2 * atanh(info);
					//info = tanh(R0[base + l - 1][i - 1] / 2) * tanh(L0[base + l - 1][i] / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] + 2 * atanh(info);

					
					//R0[base + l - 1][i] = logdomain_sum(R0[base + l - 1][i - 1] + L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1], 0)
					//                       - logdomain_sum(R0[base + l - 1][i - 1], L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					//R0[base + B / 2 + l - 1][i] = R0[base + B / 2 + l - 1][i - 1] + logdomain_sum(R0[base + l - 1][i - 1] + L0[base + l - 1][i], 0)
					//                       - logdomain_sum(R0[base + l - 1][i - 1] , L0[base + l - 1][i]);



				}
			}

		}

		for (int i = n; i >= 1; i--)
		{
			B = pow(2, (n - i + 1));
			nB = pow(2, (i - 1));
			for (int j = 1; j <= nB; j++)
			{
				base = (j - 1) * B;
				for (int l = 1; l <= B / 2; l++)
				{

					//L0[base + l - 1][i - 1] = (1 + L0[base + l - 1][i] * L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]) / (L0[base + l - 1][i] + L0[base + B / 2 + l - 1][i] * R0[base + B / 2 + l - 1][i - 1]);
					//L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] * ((1 + L0[base + l - 1][i] * R0[base + l - 1][i - 1]) / (L0[base + l - 1][i] + R0[base + l - 1][i - 1]));

					//L0 取 log ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					L0[base + l - 1][i - 1] = 0.9375*sgn(L0[base + l - 1][i])*sgn(L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1])*fabsmin(L0[base + l - 1][i], L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] + 0.9375*sgn(L0[base + l - 1][i])*sgn(R0[base + l - 1][i - 1])*fabsmin(L0[base + l - 1][i], R0[base + l - 1][i - 1]);


					//info = tanh(L0[base + l - 1][i] / 2) * tanh((L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]) / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//L0[base + l - 1][i - 1] = 2 * atanh(info);


					//info = tanh(L0[base + l - 1][i] / 2) * tanh(R0[base + l - 1][i - 1] / 2);
					//if (info == 1) info = 0.9999999999;
					//if (info == -1) info = -0.9999999999; //atanh內不能為1或-1
					//L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] + 2 * atanh(info);


					//L0[base + l - 1][i - 1] = logdomain_sum(L0[base + l - 1][i] + L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1], 0)
					//							- logdomain_sum(L0[base + l - 1][i] , L0[base + B / 2 + l - 1][i] + R0[base + B / 2 + l - 1][i - 1]);
					//L0[base + B / 2 + l - 1][i - 1] = L0[base + B / 2 + l - 1][i] + logdomain_sum(L0[base + l - 1][i] + R0[base + l - 1][i - 1], 0)
					//							- logdomain_sum(L0[base + l - 1][i] , R0[base + l - 1][i - 1]);



				}
			}
		}
		//G-matrix///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for (i = 0; i < N; i++)
		{

			//if (R0[i][n] >= 1.0) L1[i] = 0;
			if (R0[i][n] >= 0) L1[i] = 0;                        //R0 取 log
			else L1[i] = 1;

		}
		for (i = 0; i < N; i++)
		{
			if (frozenbits[i] == -1)
			{
				//if (L0[i][0] >= 1.0) R1[i][0] = 0;
				if (L0[i][0] >= 0) R1[i][0] = 0;                //L0 取 log
				else R1[i][0] = 1;
			}
			else
				R1[i][0] = frozenbits[i];
		}

		for (int i = 1; i <= n; i++)
		{
			B = pow(2, (n - i + 1));
			nB = pow(2, (i - 1));
			for (int j = 1; j <= nB; j++)
			{
				base = (j - 1) * B;
				for (int l = 1; l <= B / 2; l++)
				{
					R1[base + l - 1][i] = (R1[base + l - 1][i - 1] + R1[base + B / 2 + l - 1][i - 1]) % 2;
					R1[base + B / 2 + l - 1][i] = R1[base + B / 2 + l - 1][i - 1] % 2;
				}
			}
		}
		for (int i = 0; i < N; i++)
		{
			R2[i] = R1[i][n];
		}



		for (i = 0; i < N; i++)
		{

			if (L1[i] != R2[i])
				error++;


		}

		if (error == 0 || s >= iteration - 1)
			break;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	}
	for (i = 0; i < N; i++)                                 //軟式決策
	{
		if (frozenbits[i] == -1)
		{
			//if (L0[i][0] >= 1.0) RBit[R] = 0;
			if (L0[i][0]>= 0) RBit[R] = 0;                 //L0 取 log
			else RBit[R] = 1;
			R++;
		}
	
	}

	////////////////////////////////////////////////////////////////////////////////
	//for (i = 0; i < N; i++)
	//{
	//	if (frozenbits[i] == -1)
	//	{
	//      //if (L0[i][0]>= 0) d_hat[i] = 0;
	//		if (L0[i][0]>= 1.0) d_hat[i] = 0;
	//		else d_hat[i] = 1;
	//	}
	//	else
	//	d_hat[i] = frozenbits[i];
	//}
	//for (i = 0; i < N; i++)
	//{
	//	//printf("%f\t", d_hat[i]);
	//	if (frozenbits[i] == -1)
	//	{
	//		RBit[R] = d_hat[i];
	//		R++;
	//	}
	//}
	////////////////////////////////////////////////////////////////////////////////////
	PolarEncode(RBit, RBitCode, frozenbits);

	R = 0;
	for (i = 0; i < N; i++)
	{
		//printf("%f\t", d_hat[i]);

		if (frozenbits[i] == -1)
		{
			RBit[R] = RBitCode[i];
			R++;
		}
	}

}///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////