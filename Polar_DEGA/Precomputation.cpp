#include "InitialValue.h"

void Precomputation(int *bitreversedindices, int *index_of_first1_from_MSB, int *index_of_first0_from_MSB)
{
	int bit[n] = { 0 };
	int a, index;
	
	for (int i = 0; i < N; i++)
	{
		int z = i;
		for (int x = 0; x < n; x++)
		{
			bit[x] = z % 2;
			z /= 2;
			if (z == 0)
			{
				for (int y = x + 1; y < n; y++) bit[y] = 0;
				break;
			}
		}
		
		a = 0;
		for (int x = n - 1; x >= 0; x--)
		{
			bitreversedindices[i] = bitreversedindices[i] + pow(2, a) * bit[x];
			a++;
		}
		//printf("%d, ", bitreversedindices[i]);

		index = 1;
		for (int x = n - 1; x >= 0; x--)
		{
			if (bit[x] == 1)
			{
				index_of_first1_from_MSB[i] = index;
				break;
			}
			index_of_first1_from_MSB[i] = index;
			index++;
		}

		index = 1;
		for (int x = n - 1; x >= 0; x--)
		{
			if (bit[x] == 0)
			{
				index_of_first0_from_MSB[i] = index;
				break;
			}
			index_of_first0_from_MSB[i] = index;
			index++;
		}
	}
}