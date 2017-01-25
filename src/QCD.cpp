#include <iostream>
#include <math.h>
#include <R.h>
using namespace std;
//extern "C" __declspec(dllexport)
/*  
function beta_iteration for updating the beta. 
x is the predictor matrix, y is the reponse viable. 
beta1 is the coefficients for updating.
j is the subscript of the excluded column of x
tau is the quantile parameter.
return the vector prebeta for caculating the weighted median 
*/
//extern "C"{
//function for iteration to get the optimal solution. 
// y is the responses, x is the predictors, beta1 is the parameters, nyrow is the row length 
// of y, nxcol is the column length of x, tau is the quantile parameter.
extern "C"{
// void F77_NAME(xssort)(double*, double*, int*, int*);
void xssort_ (double*, double*, int*, int*);

// C WRAPPER FUNCTION FOR xssort SUBROUTINE
void quicksort(double *x, double *w, int *n){
	int kflag = 2;
//  CALLING FORTRAN	FUNCTION FROM C
	// F77_CALL(xssort)(x, w, n, &kflag);
	xssort_ (x, w, n, &kflag);

}

void QCD( double *x, double *beta, double *intval, double *penweight, double *residuals,
		  	  int *n, int *p, int *intercept, double *tau, double *eps, int *maxin)
// x is the design matrix (just covariates, no column of ones); n x p matrix (converted to vector)
// beta is the initial (and will be the final) value of coefficients; vector of length p
// intval is the initial (and will be the final) value of the intercept; can be anything if no intercept in model
// penweight is the penalty weight for each coefficient; vector of length p
// residuals is the vector of residuals from current value of betas; vector of length n ( = y - x%*%beta )
// n is the number of observations
// p is the number of coefficients (not including intercept)
// intercept is the intercept indicator, 1 is intercept, 0 is no intercept
// tau is the quantile of interest
// thresh is the convergence threshold
// maxin is maximum number of iterations
{	
	int iter = 0;                 // Iteration count
	int col;                  // Column index (p columns in x and then do intercept last)
	int row, rowcol, nonzerox;    // Row index, (row, col) entry of x, keep track of nonzero x entries for each column

	double *weight = new double[*n + 1]; // Weight vector
	int length_weight;                   // Length of weight vector

	double *pre_value_final=new double[*n+1]; // pre_value_final is the vector to save the weighted median vector

	// Used for computing weighted median
	int count3;
	double weight_sum, temp1;

	// Keep track whether beta vector has converged
	double betavecDiff;
	double betaDiff;
	double newBeta; // placeholder for new beta


	while( (iter < *maxin) ) //number of inside iterations
	{	
		col = 0;
		betavecDiff = 0;
		while( col < *p ) // This is the loop for the coefficients, do the intercept last
		{
			temp1=0;  // Used for calculation of weighted median
			weight_sum=0;  // sum of the weight vector
			row = 0;
			nonzerox = -1;
			rowcol = (*n)*col;
			while( row < *n )  // loop through the entries in x of column number col
			{	
				if( x[rowcol] != 0 )
				{
					nonzerox++;
					pre_value_final[nonzerox] = residuals[row];

					if( pre_value_final[nonzerox] > 0 )
						weight[nonzerox] = fabs(x[rowcol])*(*tau);
					else
						weight[nonzerox] = fabs(x[rowcol])*(1 - *tau);

					pre_value_final[nonzerox] = ( pre_value_final[nonzerox] + x[rowcol]*beta[col] )/( x[rowcol] );
					if ( *tau >= 0.5 )
						pre_value_final[nonzerox] =  -( pre_value_final[nonzerox] );

					weight_sum += weight[nonzerox];
				}

				row++;	
				rowcol++;
			}

			// Compute the extra pseudo-observation (penalty for the coefficient)
			nonzerox++;
			pre_value_final[nonzerox] = 0;

			weight[nonzerox] = penweight[col];
			weight_sum += weight[nonzerox];

			// Compute the weighted median
			length_weight = nonzerox+1;
			quicksort( pre_value_final, weight, &length_weight );
			count3 = -1;
			temp1 = 0;
			weight_sum = weight_sum/2; // Don't need sum anymore, just need sum/2
			while( temp1 < weight_sum )
			{
				count3++;
				temp1 = temp1 + weight[count3];
			}

			// Update the coefficient
			newBeta = pre_value_final[count3];
			if ( *tau >= 0.5)
				newBeta = -newBeta;

			betaDiff = beta[col] - newBeta;
			betavecDiff = betavecDiff + betaDiff*betaDiff;
			beta[col] = newBeta;

			// Update the residual vector if newBeta has changed
			if( betaDiff != 0 )
			{
				row = 0;
				rowcol = (*n)*col;
				while( row < *n )
				{
					residuals[row] = residuals[row] + x[rowcol]*betaDiff;
					row++;
					rowcol++;
				}
			}

			col++; // Next coefficient and repeat
		}

		if( *intercept == 1) // If intercept, update the intercept
		{
			temp1 = 0;  // Used for calculation of univariate tauth quantile
			weight_sum = 0;  // sum of the weight vector (for the intercept, everything has weight of 1)
			row = 0;
			nonzerox = -1;
			while( row < *n )  // loop through the entries in x of column number col
			{	
				nonzerox++;
				pre_value_final[nonzerox] = residuals[row] + intval[0];
				weight[nonzerox] = 1;
				row++;	
				weight_sum = weight_sum + weight[nonzerox];
			}

			// Compute the univariate tauth quantile
			length_weight = nonzerox+1;
			quicksort( pre_value_final, weight, &length_weight );
			count3 = -1;
			temp1 = 0;
			weight_sum = weight_sum*(*tau); // Don't need sum anymore, just need sum*tau
			while( temp1 < weight_sum )
			{
				count3++;
				temp1 = temp1 + weight[count3];
			}

			// Update the coefficient
			newBeta = pre_value_final[count3];

			betaDiff = intval[0] - newBeta;
			betavecDiff = betavecDiff + betaDiff*betaDiff;
			intval[0] = newBeta;

			// Update the residual vector if newBeta has changed
			if( betaDiff != 0 )
			{
				row = 0;
				while( row < *n )
				{
					residuals[row] = residuals[row] + betaDiff;
					row++;
				}
			}

		}


		if( sqrt(betavecDiff) < *eps ) // Check for convergence
		{
			break;
		}

		iter++; // Next iteration (col resets to 0 at the top)
	}

		// Final clean up
		delete [] pre_value_final;
		delete [] weight;

		
}



void penderiv( double *beta, int *p, double *a, double *lambda, int *pentype)
// beta is the vector of coefficients and will return the value of the derivative of the penalties
// n is the number of observations
// p is the number of coefficients
// a is a parameter for SCAD and MCP
// lambda is a parameter for SCAD, MCP, and LASSO
// pentype is the penalty (0 for SCAD, 1 for MCP, 2 for LASSO)
{
	int count = 0;
	double temp;
	while( count < *p )
	{
		if( *pentype == 0 ) // SCAD
		{
			if( fabs(beta[count]) < *lambda )
				temp = *lambda;
			else if ( fabs(beta[count]) < (*a)*(*lambda) )
				temp = ( (*a)*(*lambda) - fabs(beta[count]) )/( *a - 1.0 );
			else 
				temp = 0;
		}
		else if( *pentype == 1 ) // MCP
		{
			if ( fabs(beta[count]) < (*a)*(*lambda) )
				temp = ( *lambda - fabs(beta[count]) )/( *a );
			else 
				temp = 0;
		}
		else // LASSO
		{
			temp = *lambda;
		}

		beta[count] = temp;
		count++;
	}
}
	
}









