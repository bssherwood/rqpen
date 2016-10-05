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
void F77_NAME(xssort)(double*, double*, int*, int*);

// C WRAPPER FUNCTION FOR xssort SUBROUTINE
void quicksort(double *x, double *w, int *n){
	int kflag = 2;
//  CALLING FORTRAN	FUNCTION FROM C
	F77_CALL(xssort)(x, w, n, &kflag);
}

void QCD( double *y,  double *x,    double *beta0,  double *beta1, double *xbeta1,
		  int *nyrow, int *nxcol,   int *intercept, double *tau,  double *lambda,
		  double *a,  int *pentype, double *thresh, int *maxin)
// y is the response variable; vector of length n
// x is the design matrix (with intercept included at end if needed); n x p matrix (converted to vector)
// beta0 is the initial value of coefficients; vector of length p
// beta1 will be the final value of coefficients; vector of length p
// xbeta1 is the vector to calculate the weighted median; vector of length n ( = x%*%beta )
// nyrow is the length of y, which is n
// nxcol is the number of columns of x, which is p
// intercept is the intercept indicator, 1 for intercept, 0 for no intercept
// tau is the quantile
// lambda is the tuning parameter
// a is another tuning parameter for SCAD and MCP penalties
// pentype is the indicator for penalty types, 0 = SCAD, 1 = MCP.
// thresh is the convergence threshold
// maxin is maximal iteration numbers
{	
	int iter = 0;                 // Iteration count
	int col = 0;                  // Column index (Do intercept last)
	int row, rowcol, nonzerox;    // Row index, (row, col) entry of x, keep track of nonzero x entries for each column
	double nrow = double(*nyrow); // Turn *nyrow into double (from int)

	double *weight = new double[*nyrow + 1]; // Weight vector
	int length_weight;                       // Length of weight vector

	double *pre_value_final=new double[*nyrow+1]; // pre_value_final is the vector to save the weighted median vector

	// Used for computing weighted median
	int count3;
	double weight_sum, temp1;

	// Keep track whether beta vector has converged
	double betavecDiff=0;
	double betaDiff=0;
	double newBeta; // placeholder for new beta


	while( (iter < *maxin) ) //number of inside iterations
	{	

		temp1=0;
		// Used for calculation of weighted median
		
		weight_sum=0;
		// sum of the weight vector
		row = 0;
		nonzerox = -1;
		rowcol = (*nyrow)*col;
		while( row < *nyrow )
		// loop through the entries in x of column number col
		{	

			if( betaDiff != 0 ){
				if(col > 0)
					xbeta1[row] += betaDiff*x[rowcol - *nyrow];
				else
					xbeta1[row] += betaDiff*x[rowcol + (*nxcol-1)*(*nyrow)];
			}

			if( x[rowcol] != 0 )
			{
				nonzerox++;
				pre_value_final[nonzerox] = y[row] - xbeta1[row];

				if( pre_value_final[nonzerox] > 0 )
					weight[nonzerox] = fabs(x[rowcol])*(*tau)/nrow;
				else
					weight[nonzerox] = fabs(x[rowcol])*(1 - *tau)/nrow;
				// update the weight_i

				if ( *tau >= 0.5 )
					pre_value_final[nonzerox] = -( pre_value_final[nonzerox] + x[rowcol]*beta1[col] )/double( x[rowcol] );
				else
					pre_value_final[nonzerox] =  ( pre_value_final[nonzerox] + x[rowcol]*beta1[col] )/double( x[rowcol] );

				weight_sum += weight[nonzerox];
			}

			row++;	
			rowcol++;
		}
		if ( (col < (*nxcol -1)) || (*intercept == 0) ) // Not Intercept coefficient
		{
			nonzerox++;
			pre_value_final[nonzerox] = 0;

			/////////////////////////////////////////////////SCAD
			if ( *pentype == 0 )
			{
      			if ( fabs(beta0[col]) < *lambda )
					weight[nonzerox] = *lambda;
				else if ( fabs(beta0[col]) < (*a)*(*lambda) )
					weight[nonzerox] = ( (*a)*(*lambda) - fabs(beta0[col]) )/(double(*a-1));
				else
					weight[nonzerox] = 0;
			}
			/////////////////////////////////////////////////
			/////////////////////////////////////////////////MCP
			else if ( *pentype == 1 )
			{
				if ( fabs(beta0[col]) < (*a)*(*lambda) )
					weight[nonzerox] = *lambda - fabs(beta0[col])/ *a;
				else
					weight[nonzerox] = 0;
			/////////////////////////////////////////////////
			}
			weight_sum += weight[nonzerox];
			length_weight = nonzerox+1;
			quicksort( pre_value_final, weight, &length_weight );
			for (count3=0 ; count3 <= nonzerox; count3++)
			{
				temp1 += weight[count3];
				if ( temp1 > weight_sum/2 )
				{
					break;
				}
			}
		}
		else // Intercept coefficient
		{
			length_weight = nonzerox+1;
			quicksort( pre_value_final, weight, &length_weight );
			for ( count3=0; count3 <= nonzerox; count3++ )
			{
				temp1 += weight[count3];
				if (temp1> weight_sum/2)
				{
					break;
				}
			}
		}

		if ( *tau >= 0.5)
			newBeta = -pre_value_final[count3];
		else
			newBeta =  pre_value_final[count3];
		//change from current beta1 to new beta1

		if( (fabs(newBeta) < *thresh) && ( (col < (*nxcol -1)) || (*intercept == 0) ) )
		{
			newBeta = 0;
		}

		betaDiff = newBeta - beta1[col];
		beta1[col] = newBeta;


		betavecDiff += fabs(betaDiff);
		col++;

		if( col == *nxcol )
		{
			if( sqrt(betavecDiff) < *thresh )
			{
				break;
			}

			iter++;
			col = 0;
			betavecDiff = 0;
		}

		
	}
	
		delete [] pre_value_final;
		delete [] weight;
}







void QCDgroup( double *y,  double *x,    double *beta0,  double *beta1, double *xbeta1,
		  int *nyrow, int *nxcol,   int *intercept, double *tau, double *lambda,
		  double *a,  int *pentype, double *thresh, int *maxin, double *groupl1)
// y is the response variable; vector of length n
// x is the design matrix (with intercept included at end if needed); n x p matrix (converted to vector)
// beta0 is the initial value of coefficients; vector of length p
// beta1 will be the final value of coefficients; vector of length p
// xbeta1 is the vector to calculate the weighted median; vector of length n ( = x%*%beta )
// nyrow is the length of y, which is n
// nxcol is the number of columns of x, which is p
// intercept is the intercept indicator, 1 for intercept, 0 for no intercept
// tau is the quantile
// lambda is the tuning parameter
// a is another tuning parameter for SCAD and MCP penalties
// pentype is the indicator for penalty types, 0 = SCAD, 1 = MCP, 2 = LASSO.
// thresh is the convergence threshold
// maxin is maximal iteration numbers
// groupl1 is vector of length p of group l1 norm for each covariate (not unique as each member of a group has the same group l1 norm)
{	
	int iter = 0;                 // Iteration count
	int col = 0;                  // Column index (Do intercept last)
	int row, rowcol, nonzerox;    // Row index, (row, col) entry of x, keep track of nonzero x entries for each column
	double nrow = double(*nyrow); // Turn *nyrow into double (from int)

	double *weight = new double[*nyrow + 1]; // Weight vector
	int length_weight;                       // Length of weight vector

	double *pre_value_final=new double[*nyrow+1]; // pre_value_final is the vector to save the weighted median vector

	// Used for computing weighted median
	int count3;
	double weight_sum, temp1;

	// Keep track whether beta vector has converged
	double betavecDiff=0;
	double betaDiff=0;
	double newBeta; // placeholder for new beta


	while( (iter < *maxin) ) //number of inside iterations
	{	

		temp1=0;
		// Used for calculation of weighted median
		
		weight_sum=0;
		// sum of the weight vector
		row = 0;
		nonzerox = -1;
		rowcol = (*nyrow)*col;
		while( row < *nyrow )
		// loop through the entries in x of column number col
		{	

			if( betaDiff != 0 ){
				if(col > 0)
					xbeta1[row] += betaDiff*x[rowcol - *nyrow];
				else
					xbeta1[row] += betaDiff*x[rowcol + (*nxcol-1)*(*nyrow)];
			}

			if( x[rowcol] != 0 )
			{
				nonzerox++;
				pre_value_final[nonzerox] = y[row] - xbeta1[row];

				if( pre_value_final[nonzerox] > 0 )
					weight[nonzerox] = fabs(x[rowcol])*(*tau)/nrow;
				else
					weight[nonzerox] = fabs(x[rowcol])*(1 - *tau)/nrow;
				// update the weight_i

				if ( *tau >= 0.5 )
					pre_value_final[nonzerox] = -( pre_value_final[nonzerox] + x[rowcol]*beta1[col] )/double( x[rowcol] );
				else
					pre_value_final[nonzerox] =  ( pre_value_final[nonzerox] + x[rowcol]*beta1[col] )/double( x[rowcol] );

				weight_sum += weight[nonzerox];
			}

			row++;	
			rowcol++;
		}
		if ( (col < (*nxcol -1)) || (*intercept == 0) ) // Not Intercept coefficient
		{
			nonzerox++;
			pre_value_final[nonzerox] = 0;

			/////////////////////////////////////////////////SCAD
			if ( *pentype == 0 )
			{
      			if ( fabs(beta0[col]) < *lambda )
					weight[nonzerox] = *lambda;
				else if ( fabs(beta0[col]) < (*a)*(*lambda) )
					weight[nonzerox] = ( (*a)*(*lambda) - groupl1[col] )/(double(*a-1));
				else
					weight[nonzerox] = 0;
			}
			/////////////////////////////////////////////////
			/////////////////////////////////////////////////MCP
			else if ( *pentype == 1 )
			{
				if ( fabs(beta0[col]) < (*a)*(*lambda) )
					weight[nonzerox] = *lambda - groupl1[col]/ *a;
				else
					weight[nonzerox] = 0;
			/////////////////////////////////////////////////
			}
			/////////////////////////////////////////////////LASSO
			else if ( *pentype == 2 )
			{
				weight[nonzerox] = *lambda;
			}
			/////////////////////////////////////////////////

			weight_sum += weight[nonzerox];
			length_weight = nonzerox+1;
			quicksort( pre_value_final, weight, &length_weight );
			for (count3=0 ; count3 <= nonzerox; count3++)
			{
				temp1 += weight[count3];
				if ( temp1 > weight_sum/2 )
				{
					break;
				}
			}
		}
		else // Intercept coefficient
		{
			length_weight = nonzerox+1;
			quicksort( pre_value_final, weight, &length_weight );
			for ( count3=0; count3 <= nonzerox; count3++ )
			{
				temp1 += weight[count3];
				if (temp1> weight_sum/2)
				{
					break;
				}
			}
		}

		if ( *tau >= 0.5)
			newBeta = -pre_value_final[count3];
		else
			newBeta =  pre_value_final[count3];
		//change from current beta1 to new beta1

		if( (fabs(newBeta) < *thresh) && ( (col < (*nxcol -1)) || (*intercept == 0) ) )
		{
			newBeta = 0;
		}

		betaDiff = newBeta - beta1[col];
		beta1[col] = newBeta;


		betavecDiff += fabs(betaDiff);
		col++;

		if( col == *nxcol )
		{
			if( sqrt(betavecDiff) < *thresh )
			{
				break;
			}

			iter++;
			col = 0;
			betavecDiff = 0;
		}

		
	}
	
		delete [] pre_value_final;
		delete [] weight;
}




}


// void QCDgrouplasso( double *y,  double *x,    double *beta0,  double *beta1, double *xbeta1,
// 		  int *nyrow, int *nxcol,   int *intercept, double *tau, double *lambda,
// 		  double *thresh, int *maxin, double *groupl1)
// // y is the response variable; vector of length n
// // x is the design matrix (with intercept included at end if needed); n x p matrix (converted to vector)
// // beta0 is the initial value of coefficients; vector of length p
// // beta1 will be the final value of coefficients; vector of length p
// // xbeta1 is the vector to calculate the weighted median; vector of length n ( = x%*%beta )
// // nyrow is the length of y, which is n
// // nxcol is the number of columns of x, which is p
// // intercept is the intercept indicator, 1 for intercept, 0 for no intercept
// // tau is the quantile
// // lambda is the tuning parameter
// // pentype is the indicator for penalty types, 0 = SCAD, 1 = MCP.
// // thresh is the convergence threshold
// // maxin is maximal iteration numbers
// // groupl1 is vector of length p of group l1 norm for each covariate (not unique as each member of a group has the same group l1 norm)
// {	
// 	int iter = 0;                 // Iteration count
// 	int col = 0;                  // Column index (Do intercept last)
// 	int row, rowcol, nonzerox;    // Row index, (row, col) entry of x, keep track of nonzero x entries for each column
// 	double nrow = double(*nyrow); // Turn *nyrow into double (from int)

// 	double *weight = new double[*nyrow + 1]; // Weight vector
// 	int length_weight;                       // Length of weight vector

// 	double *pre_value_final=new double[*nyrow+1]; // pre_value_final is the vector to save the weighted median vector

// 	// Used for computing weighted median
// 	int count3;
// 	double weight_sum, temp1;

// 	// Keep track whether beta vector has converged
// 	double betavecDiff=0;
// 	double betaDiff=0;
// 	double newBeta; // placeholder for new beta


// 	while( (iter < *maxin) ) //number of inside iterations
// 	{	

// 		temp1=0;
// 		// Used for calculation of weighted median
		
// 		weight_sum=0;
// 		// sum of the weight vector
// 		row = 0;
// 		nonzerox = -1;
// 		rowcol = (*nyrow)*col;
// 		while( row < *nyrow )
// 		// loop through the entries in x of column number col
// 		{	

// 			if( betaDiff != 0 ){
// 				if(col > 0)
// 					xbeta1[row] += betaDiff*x[rowcol - *nyrow];
// 				else
// 					xbeta1[row] += betaDiff*x[rowcol + (*nxcol-1)*(*nyrow)];
// 			}

// 			if( x[rowcol] != 0 )
// 			{
// 				nonzerox++;
// 				pre_value_final[nonzerox] = y[row] - xbeta1[row];

// 				if( pre_value_final[nonzerox] > 0 )
// 					weight[nonzerox] = fabs(x[rowcol])*(*tau)/nrow;
// 				else
// 					weight[nonzerox] = fabs(x[rowcol])*(1 - *tau)/nrow;
// 				// update the weight_i

// 				if ( *tau >= 0.5 )
// 					pre_value_final[nonzerox] = -( pre_value_final[nonzerox] + x[rowcol]*beta1[col] )/double( x[rowcol] );
// 				else
// 					pre_value_final[nonzerox] =  ( pre_value_final[nonzerox] + x[rowcol]*beta1[col] )/double( x[rowcol] );

// 				weight_sum += weight[nonzerox];
// 			}

// 			row++;	
// 			rowcol++;
// 		}
// 		if ( (col < (*nxcol -1)) || (*intercept == 0) ) // Not Intercept coefficient
// 		{
// 			nonzerox++;
// 			pre_value_final[nonzerox] = 0;

// 			/////////////////////////////////////////////////LASSO
// 			weight[nonzerox] = *lambda;
// 			/////////////////////////////////////////////////
// 			weight_sum += weight[nonzerox];
// 			length_weight = nonzerox+1;
// 			quicksort( pre_value_final, weight, &length_weight );
// 			for (count3=0 ; count3 <= nonzerox; count3++)
// 			{
// 				temp1 += weight[count3];
// 				if ( temp1 > weight_sum/2 )
// 				{
// 					break;
// 				}
// 			}
// 		}
// 		else // Intercept coefficient
// 		{
// 			length_weight = nonzerox+1;
// 			quicksort( pre_value_final, weight, &length_weight );
// 			for ( count3=0; count3 <= nonzerox; count3++ )
// 			{
// 				temp1 += weight[count3];
// 				if (temp1> weight_sum/2)
// 				{
// 					break;
// 				}
// 			}
// 		}

// 		if ( *tau >= 0.5)
// 			newBeta = -pre_value_final[count3];
// 		else
// 			newBeta =  pre_value_final[count3];
// 		//change from current beta1 to new beta1

// 		if( (fabs(newBeta) < *thresh) && ( (col < (*nxcol -1)) || (*intercept == 0) ) )
// 		{
// 			newBeta = 0;
// 		}

// 		betaDiff = newBeta - beta1[col];
// 		beta1[col] = newBeta;


// 		betavecDiff += fabs(betaDiff);
// 		col++;

// 		if( col == *nxcol )
// 		{
// 			if( sqrt(betavecDiff) < *thresh )
// 			{
// 				break;
// 			}

// 			iter++;
// 			col = 0;
// 			betavecDiff = 0;
// 		}

		
// 	}
	
// 		delete [] pre_value_final;
// 		delete [] weight;
// }
