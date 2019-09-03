//Modified Bessel function of order 0
//multipled by exponential as it appears in gbarvv
//from Numerical Recipes in C second edition, available online

#include <math.h>
#include <stdio.h>

float expBessi0(float x)
//Returns exp(-x)*I0(x) for any real x >= 0.
{
	float ans;
	double y;  //Accumulate polynomials in double precision.
	if(x < 0) printf("*** error in expBessi0, argment <0 ***\n");
	else if(x < 3.75){ //Polynomial fit.
		y = x/3.75;
		y*=y;
		ans = 1. + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
			+ y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
		ans *= exp(-x);
	}
	else{
		y = 3.75/x;
		ans = (1./sqrt(x))*(0.39894228 + y*(0.1328592e-1
			+y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}

