#include "InitialValue.h"
#include <iostream>

void PolarEncode(double *Message, double *Code, int *frozenbits)
{
	int B, nB, base;
	int m = 0;

	for (int i = 0; i < N; i++)
	{
		Code[i] = frozenbits[i];
		if (Code[i] == -1)
		{
			Code[i] = Message[m];
			m++;
		}
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
				int x = Code[base + l - 1] + Code[base + B / 2 + l - 1];
				Code[base + l - 1] = x % 2;
			}
		}
	}
}

double min(double x, double y)
{
	double z;
	if (x > y)
		z = y;
	else
		z = x;

	return z;
}
void SystematicPolarEncode(double *Message, double *Code, int *frozenbits)
{
	int m = 0;
	int s, delta;
	double y[N] = { 0 };
	int z[N][n + 1] = { 0 };
	int bit[n];

	for (int i = 0; i < N; i++)
	{
		Code[i] = frozenbits[i];

		if (Code[i] == -1)
		{
			Code[i] = Message[m];
			m++;
		}
		else
		{
			Code[i] = -1;
		}
	}

	for (int i = 0; i < N; i++)
	{
		y[i] = frozenbits[i];
	}

	for (int i = 0; i < N; i++)
	{
		z[i][0] = y[i];
		z[i][n] = Code[i];
	}

	for (int i = N; i >= 1; i--)
	{
		if (frozenbits[i - 1] == -1)
		{
			s = n + 1;
			delta = -1;
		}
		else
		{
			s = 1;
			delta = 1;
		}

		int a = i - 1;
		for (int x = n - 1; x >= 0; x--)
		{
			bit[x] = a % 2;
			a /= 2;
			if (z == 0)
			{
				for (int y = x + 1; y < n; y++) bit[y] = 0;
				break;
			}
		}

		for (int j = 1; j <= n; j++)
		{
			int t = s + j*delta;
			int l = min(t, t - delta);
			int k = pow(2, n - l);

			if (bit[l - 1] == 0)
				z[i - 1][t - 1] = (z[i - 1][t - delta - 1] + z[i + k - 1][t - delta - 1]) % 2;
			else
				z[i - 1][t - 1] = z[i - 1][t - delta - 1];
		}
	}

	for (int i = 0; i < N; i++)
	{
		Code[i] = z[i][n];
	}
}