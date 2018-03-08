#include "calcMatrix.h"
#include <math.h>

bool getMatrixInverse(double** src, int n, double** dst)
{
	int i, j;
	double flag = getMatrixDet(src, n);
	double** t = new double*[n];
	for(i = 0; i < n; i++)
		t[i] = new double[n];
	if(flag == 0)
		return false;
	else
	{
		getMatrixStar(src, n, t);
		for(i = 0; i < n; i++)
			for(j = 0; j < n; j++)
				dst[i][j] = t[i][j] / flag;
		return true;
	}
	
	for(i = 0; i < n; i++)
		delete[] t[i];
	delete[] t;
}

double getMatrixDet(double** src, int n)
{
	int i, j, k;
	if(n == 1)
		return src[0][0];
	double ans = 0;
	double** temp = new double*[n];
	for(i = 0; i < n; i++)
		temp[i] = new double[n];

	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			temp[i][j] = 0.0;

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n - 1; j++)
			for(k = 0; k < n - 1; k++)
				temp[j][k] = src[j + 1][(k >= i)? k + 1 : k];
		double t = getMatrixDet(temp, n - 1);
		if(i % 2 == 0)
			ans += src[0][i] * t;
		else
			ans -= src[0][i] * t;
	}

	for(i = 0; i < n; i++)
		delete[] temp[i];
	delete[] temp;

	return ans;
}

void getMatrixStar(double** src, int n, double** dst)
{
	if(n == 1)
	{
		dst[0][0] =  1;
		return;
	}

	int i, j, k, t;
	double** temp = new double*[n];
	for(i = 0; i < n; i++)
		temp[i] = new double[n];
	
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
		{
			for(k = 0; k < n - 1; k++)
				for(t = 0; t < n - 1; t++)
					temp[k][t] = src[(k >= i)? k + 1 : k][(t >= j)? t + 1 : t];
			
			dst[j][i] = getMatrixDet(temp, n - 1);
			if((i + j) % 2 == 1)
				dst[j][i] = -dst[j][i];
		}

	for(i = 0; i < n; i++)
		delete[] temp[i];
	delete[] temp;
}

void getMatrixSqrt(double** src, int n, double** dst1, double** dst2)
{
	int i, j;
	double** Y = new double*[n];
	double** Z = new double*[n];
	double** invY = new double*[n];
	double** invZ = new double*[n];
	double** temp = new double*[n];
	double** temp1 = new double*[n];
	for(i = 0; i < n; i++)
	{
		Y[i] = new double[n];
		Z[i] = new double[n];
		invY[i] = new double[n];
		invZ[i] = new double[n];
		temp[i] = new double[n];
		temp1[i] = new double[n];
	}

	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
		{
			Y[i][j] = src[i][j];
			if(i == j)
				Z[i][j] = 1.0;
			else
				Z[i][j] = 0.0;
		}
	
	double eps = 1e-5;
	double ans = 1.0;
	while(ans > eps)
	{
		getMatrixEqual(Z, n, temp);
		getMatrixInverse(Z, n, invZ);
		getMatrixInverse(Y, n, invY);
		for(i = 0; i < n; i++)
			for(j = 0; j < n; j++)
			{
				Y[i][j] = (Y[i][j] + invZ[i][j]) / 2;
				Z[i][j] = (Z[i][j] + invY[i][j]) / 2;
			}
	    getMatrixSubstract(Z, temp, n, temp1);
		ans = fabs(getMatrixDet(temp1, n));
	}

	getMatrixEqual(Y, n, dst1);
	getMatrixEqual(Z, n, dst2);

	for(i = 0; i < n; i++)
	{
		delete[] Y[i], Z[i], invY[i], invZ[i], temp[i], temp1[i];
	}
	delete[] Y, Z, invY, invZ, temp, temp1;
}

void getMatrixSubstract(double** src1, double** src2, int n, double** dst)
{
	for(int i = 0; i < n; i++)
		for(int j =0; j < n; j++)
			dst[i][j] = src1[i][j] - src2[i][j];
}

void getMatrixAdd(double** src1, double** src2, int n, double** dst)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			dst[i][j] = src1[i][j] + src2[i][j];
}

void getMatrixEqual(double **src, int n, double** dst)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			dst[i][j] = src[i][j];
}

void getMatrixMultiply(double** src1, double** src2, int m, int n, int p, double** dst)
{
	double* row = new double[n];
	double* column = new double[n];
	for(int i = 0; i < m; i++)
		for(int j = 0; j < p; j++)
		{
			getMatrixRow(src1, m, n, i, row);
			getMatrixColumn(src2, n, p, j, column);
			dst[i][j] = getVectorMultiply(row, column, n);
		}
	delete[] row, column;
}

void getMatrixScale(double** src, int m, int n, double scale, double** dst)
{
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			dst[i][j] = scale * src[i][j];
}

double getVectorMultiply(double* vec1, double* vec2, int n)
{
	double ans = 0;
	for(int i = 0; i < n; i++)
		ans += vec1[i] * vec2[i];
	return ans;
}

void getVectorScale(double* vec, int n, double scale, double* dst)
{
	for(int i = 0; i < n; i++)
		dst[i] = scale * vec[i];
}

void getMatrixRow(double** src, int m, int n, int i, double* irow)
{
	for(int j = 0; j < n; j++)
		irow[j] = src[i][j];
}

void getMatrixColumn(double** src, int m, int n, int i, double* icolumn)
{
	for(int j = 0; j < m; j++)
		icolumn[j] = src[j][i];
}

void getVectorCrossMultiply(double* vec1, double* vec2, double* dst)
{
	dst[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	dst[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	dst[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

double getVectorNorm(double* vec, int n)
{
	double elemSquare = 0;
	for(int i = 0; i < n; i++)
		elemSquare += vec[i] * vec[i];
	return sqrt(elemSquare);
}


















