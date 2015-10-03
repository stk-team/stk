#include <stk/stk.hpp>

int g_nPts;

extern "C" void module_sampler_init(int nPts)
{
	g_nPts = nPts;
}

extern "C" void module_sampler_exec(stk::PointSet2dd& pts)
{
	const int D = 2;
	
	// L = max number of bits needed 
	unsigned L = (unsigned)std::ceil(std::log((double)g_nPts)/log(2.0)); 

	// C[i] = index from the right of the first zero bit of i
	unsigned *C = new unsigned [g_nPts];
	C[0] = 1;
	for (unsigned i=1;i<=g_nPts-1;i++)
	{
		C[i] = 1;
		unsigned value = i;
		while (value & 1)
		{
			value >>= 1;
			C[i]++;
		}
	}

	// POINTS[i][j] = the jth component of the ith point
	//                with i indexed from 0 to N-1 and j indexed from 0 to D-1
	double *POINTS = new double[g_nPts*D];
	POINTS[0] = 0.0;
	POINTS[1] = 0.0;

	// ----- Compute the first dimension -----

	// Compute direction numbers V[1] to V[L], scaled by pow(2,32)
	unsigned *V = new unsigned [L+1]; 
	for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

	// Evalulate X[0] to X[N-1], scaled by pow(2,32)
	unsigned *X = new unsigned [g_nPts];
	X[0] = 0;
	for (unsigned i=1;i<=g_nPts-1;i++)
	{
		X[i] = X[i-1] ^ V[C[i-1]];
		POINTS[i*D+0] = (double)X[i]/pow(2.0,32); // *** the actual points
		//        ^ 0 for first dimension
	}

	// Clean up
	delete [] V;
	delete [] X;

	//Computing remaning dimensions (just the second one in fact)
	{
		const int j=1;
		const unsigned d=2;
		const unsigned s=1;
		const unsigned a=0;

		// Read in parameters from file
		unsigned *m = new unsigned [s+1];
		m[1] = 1;

		// Compute direction numbers V[1] to V[L], scaled by pow(2,32)
		unsigned *V = new unsigned [L+1];
		if (L <= s)
		{
			for (unsigned i=1;i<=L;i++)
			{
				V[i] = m[i] << (32-i); 
			}
		}
		else
		{
			for (unsigned i=1;i<=s;i++)
			{
				V[i] = m[i] << (32-i);
			}
			for (unsigned i=s+1;i<=L;i++)
			{
				V[i] = V[i-s] ^ (V[i-s] >> s); 
				
				for (unsigned k=1;k<=s-1;k++) 
				{
					V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
				}
			}
		}

		// Evalulate X[0] to X[N-1], scaled by pow(2,32)
		unsigned *X = new unsigned [g_nPts];
		X[0] = 0;
		for (unsigned i=1;i<=g_nPts-1;i++)
		{
			X[i] = X[i-1] ^ V[C[i-1]];
			POINTS[i*D+j] = (double)X[i]/pow(2.0,32); // *** the actual points
			//        ^ j for dimension (j+1)
		}

		// Clean up
		delete [] m;
		delete [] V;
		delete [] X;
	}
	
	delete [] C;
	
	for(int i=0; i<g_nPts; i++)
	{
		pts.push_back(stk::Point2dd(stk::Vector2d(POINTS[i*D+0], POINTS[i*D+1]), 1.0));
	}
	
	delete [] POINTS;
}

extern "C" void module_sampler_free(int nPts)
{
	
}
