#include "InitialValue.h"
#include "function.h"
#include <iostream>
#include <math.h>

double SNR[26] = { 0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0 };
double A_real_scan[9] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };               //AWAN�{�ɭȱ��y
double p_scan[9] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };                    //BG�{�ɭȱ��y


float sigma;
int frozenbits[N] = { 0 };

void main()
{
	FILE *error;
	int a, i, sum, sum1, fsum, l, D;
	double Message[K] = { 0 };
	double RBit[K] = { 0 };
	double Code[N] = { 0 };
	double RCV[N] = { 0 };
	int bitreversedindices[N] = { 0 };
	int index_of_first1_from_MSB[N] = { 0 };
	int index_of_first0_from_MSB[N] = { 0 };
	int indices[N] = { 0 };
	double BG[N] = { 0 };                        // Bernoulli-Gauss�߽Ħ�m
	int AWANstate[N] = { 0 };					 // AWAN�߽Ħ�m
	double p_prob[m_length + 1] = { 0 };
	double sigmam[m_length + 1] = { 0 };
	double c_prob[m_length + 1] = { 0 };
	double am_prob[m_length + 1] = { 0 };
	double variance, BER, rndm;

	/*error = fopen("05known.txt", "w");*/
	srand(time(NULL));

	Precomputation(bitreversedindices, index_of_first1_from_MSB, index_of_first0_from_MSB);
	//frozenlookup(indices, frozenbits, BG, 0, 0, AWANstate);


	for (a = 0; a < 10; a++)
	{
		//error = fopen("SPolar_1024_QPSK_BG0d1_BA0dB_RxKnownPosition.txt", "a");
		error = fopen("TEST_1.txt", "a");
		sum = 0;
		fsum = 0;

		variance = (float)(1 / (2.0 * Rm * Rate * pow(10, (SNR[a] / 10))));       // = (sigma) ^ 2   sigma_g

		// AWAN Impulse Model////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for (int m_count = 0; m_count <= m_length; m_count++)                               //poisson distribution
		{
			p_prob[m_count] = exp(-A_real)*pow(A_real, m_count) / factorial(m_count);		//array�����C��index�ӧO�x�sm=0,1,...10���� 
			//printf("p_prob[%d] = %.5f\n",m_count,p_prob[m_count]);
			am_prob[m_count] = pow(A_real, m_count) / factorial(m_count);                   //poisson distribution  �Lexp        
			//printf("am_prob[%d] = %.5f\n",m_count,am_prob[m_count]);
		}

		for (int m_count = 0; m_count <= m_length; m_count++)        //c_prob��, �C��index�s�J����m=0,1,2,...10��Poisson prob 
		{
			c_prob[m_count] = 0;                                      //poisson distribution m=0 awgn m=1~�L���j
			for (int j = 0; j <= m_count; j++)
				c_prob[m_count] += p_prob[j];                         //poisson distribution ���[
		}

		for (int m_count = 0; m_count <= m_length; m_count++)
		{
			sigmam[m_count] = 0;
			sigmam[m_count] = sqrt(variance + (m_count / A_real)*(variance / GIR_real));	//��²�L�᪺�Φ� �mm , ���O����         sigma_m���ڸ� 
																							//printf("sigma_%d = %.5f\n",m_count,sigmam[m_count]);
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for (l = 1; l <= loop; l++)
		{
			sum1 = 0;

			// �w���yBG QPSK /////////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int i = 0; i < N; i = i + 2)
			{
				rndm = (double)rand() / RAND_MAX;   //�ü�
				BG[i] = rndm;						//�q�D���A
				BG[i + 1] = BG[i];					 //qpsk 2BITS �_��
			}
			// �w���yAWAN QPSK ///////////////////////////////////////////////////////////////////////////////////////////////
			for (int A = 0; A < N; A = A + 2)
			{
				int m_pa = 0;
				double R_V = rand() / (double)(RAND_MAX);

				while ((R_V > c_prob[m_pa]) && (m_pa != (m_length))) //sigma_m�֥[���� = m_length
				{
					m_pa = m_pa + 1;                     //SIGMA_M ���v�֥[ =c_PROB 
				}
				AWANstate[A] = m_pa;                      
				AWANstate[A + 1] = m_pa;                     //�]��qpsk�ҥH2BITS�@��
			}
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			frozenlookup(indices, frozenbits, BG, SNR[a], a, AWANstate);

			for (i = 0; i < K; i++)
			{
				Message[i] = (rand() % 2);
			}

			//PolarEncode(Message, Code, frozenbits);
			SystematicPolarEncode(Message, Code, frozenbits);

			//BPSK_AWGN(Code, RCV, N, SNR[a]);
			//QPSK_AWGN_channel(Code, RCV, N, SNR[a]);
			QPSK_BernoulliGauss(Code, RCV, N, SNR[a], BG);
			//QPSK_AWAN(Code, RCV, N, SNR[a], c_prob, p_prob, sigmam, AWANstate);
			
			//PolarDecode(RCV, RBit, SNR[a], frozenbits, bitreversedindices, index_of_first1_from_MSB, index_of_first0_from_MSB, BG);
			//SystematicPolarDecode(RCV, RBit, SNR[a], frozenbits, bitreversedindices, index_of_first1_from_MSB, index_of_first0_from_MSB, BG, p_prob, sigmam);
			//SystematicPolarDecode_BG(RCV, RBit, SNR[a], frozenbits, bitreversedindices, index_of_first1_from_MSB, index_of_first0_from_MSB, BG, p_scan[a]);
			//SystematicPolarDecode_AWAN_QPSK(RCV, RBit, SNR[a], frozenbits, bitreversedindices, index_of_first1_from_MSB, index_of_first0_from_MSB, am_prob, sigmam, AWANstate, A_real_scan[a]);
			SystematicBPPolarDecode_BG(RCV, RBit, SNR[a], frozenbits, bitreversedindices, index_of_first1_from_MSB, index_of_first0_from_MSB, BG, p_scan[a]);
			//SystematicBPPolarDecode_AWAN(RCV, RBit, SNR[a], frozenbits, bitreversedindices, index_of_first1_from_MSB, index_of_first0_from_MSB, am_prob, sigmam, AWANstate, c_prob, A_real_scan[a]);

			
			for (i = 0; i < K; i++)
			{
				if (Message[i] != RBit[i])
				{
					sum++;
					//sum1++;               //FER
				}
			}
			//if (sum1 > 0) fsum++;         //FER

			BER = (double)sum / (K * l);
			printf("\r%1f\t %.2E\t %d\t %d", (double)SNR[a], BER, l, sum);                            //�H�۰j��@����s    //BER
			//printf("\r%1f\t %.2E\t %d\t %d", (double)A_real_scan[a], BER, l, sum);�@�@�@�@�@�@�@�@�@//AWAN�@A sacn
			//printf("\r%1f\t %.2E\t %d\t %d", (double)p_scan[a], BER, l, sum);�@�@�@�@�@�@�@�@�@�@�@ //BG p scan�@�@�@�@�@�@�@

			//fram error rate /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//printf("\r%1f\t %.2E\t %d\t %d", (double)SNR[a], (double)fsum / (l + 1), l + 1, fsum);   //�H�۰j��@����s    //frame error rate
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//printf("%1f\t %.2E\t %d\t %d\t %d\n", (double)SNR[a], (double)sum / (K * (l + 1)), l + 1, sum, D);

			if (sum > 8000 && l > 8000) break;
			//if (fsum > 10000 && l+1 > 10000) break;
		}
		
		fprintf(error, "%1f\t %.9lf\n", (double)SNR[a], BER);                                     //��X���ɮ�      //BER		
		//fprintf(error, "%1f\t %.9lf\n", (double)A_real_scan[a], BER);                           //AWAN�@A sacn
		//fprintf(error, "%1f\t %.9lf\n", (double)p_scan[a], BER);                                //BG p scan

		
		//frame error rate /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//fprintf(error, "%1f\t %.8lf\n", (double)SNR[a], (double)fsum / l);                      //��X���ɮ�      //frame error rate
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		printf("\r%1f\t %.2E\t %d\n", (double)SNR[a], BER, l);                                    //�̲׭�      //BER	
		//printf("\r%1f\t %.2E\t %d\n", (double)A_real_scan[a], BER, l);�@�@�@�@�@�@�@�@�@�@�@�@�@//AWAN�@A sacn
		//printf("\r%1f\t %.2E\t %d\n", (double)p_scan[a], BER, l);�@�@�@�@�@�@�@�@�@�@�@�@�@�@�@ //BG p scan

		//frame error rate /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printf("\r%1f\t %.2E\t %d\t %d\n", (double)SNR[a], (double)fsum / l, l, fsum);          //�̲׭�   //frame error rate              
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		fclose(error);

	}
	/*fclose(error);*/
	system("pause");

}