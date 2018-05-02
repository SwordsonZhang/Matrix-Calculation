#ifndef CALCMATRIX_GENERAL_H_
#define CALCMATRIX_GENERAL_H_

#include <math.h>

/****************************************************************************************
* Struct: 
*	typePoint2D: 2-Dimensional point in coordinate(x, y)
* 
* Variables:
*	T x: x-coordinate  
*	T y: y-coordinate
****************************************************************************************/
template<typename T>
struct typePoint2D
{
	T x;
	T y;
};


/****************************************************************************************
* Struct:
*	typePoint3D: 3-Dimensional point in coordinate(x, y, z)
* 
* Variables:
*	T x: x-coordinate  
*	T y: y-coordinate
*	T z: z-coordinate
****************************************************************************************/
template<typename T>
struct typePoint3D
{
	T x;
	T y;
	T z;
};


/***************************************************************************************
* Function
*	getMatrixInverse(): to calculate the inverse matrix
*
* Parameters
*	T** src: input matrix, square matrix
*	int n: dimension of input matrix
*	float** dst: output matrix, inverse square matrix
*
* Return
*	bool: inversed if true, or false if failed
***************************************************************************************/
template<typename T>
bool getMatrixInverse(T** src, int n, float** dst)
{
	int i, j;
	float flag = getMatrixDet(src, n);
	float** t = new T*[n];
	for(i = 0; i < n; i++)
		t[i] = new T[n];
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



/**************************************************************************************
* Function
*	getMatrixDet(): to calculate determinant of src
*
* Parameters
*	T** src: input matrix
*	int n: dimension of src
*
* Return
* 	T: value of the determinant
**************************************************************************************/
template<typename T>
float getMatrixDet(T** src, int n)
{
	int i, j, k;
	if(n == 1)
		return src[0][0];
	float ans = 0;
	T** temp = new T*[n];
	for(i = 0; i < n; i++)
		temp[i] = new T[n];

	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			temp[i][j] = 0.0;

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n - 1; j++)
			for(k = 0; k < n - 1; k++)
				temp[j][k] = src[j + 1][(k >= i)? k + 1 : k];
		float t = getMatrixDet(temp, n - 1);
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


/**************************************************************************************
* Function
*	getMatrixStar(): to calculate adjoint matrix of src
*
* Parameters
*	T** src: input square matrix
*	int n: dimension of src
*	float** dst: adjoint matrix of src
*
* Return
*	void
**************************************************************************************/
template<typename T>
void getMatrixStar(T** src, int n, float** dst)
{
	if(n == 1)
	{
		dst[0][0] =  1;
		return;
	}

	int i, j, k, t;
	T** temp = new T*[n];
	for(i = 0; i < n; i++)
		temp[i] = new T[n];
	
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


/*************************************************************************************
* Function
*	getMatrixSqrt(): to calculate square matrix of src
*
* Parameters
*	T** src: input matrix
*	int n: dimension of src
*	T** dst1: 1/2 of input matrix
*	T** dst2: -1/2 of input matrix
*
* Return
*	void
*************************************************************************************/
template<typename T>
void getMatrixSqrt(T** src, int n, T** dst1, T** dst2)
{
	int i, j;
	T** Y = new T*[n];
	T** Z = new T*[n];
	T** invY = new T*[n];
	T** invZ = new T*[n];
	T** temp = new T*[n];
	T** temp1 = new T*[n];
	for(i = 0; i < n; i++)
	{
		Y[i] = new T[n];
		Z[i] = new T[n];
		invY[i] = new T[n];
		invZ[i] = new T[n];
		temp[i] = new T[n];
		temp1[i] = new T[n];
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
	
	T eps = 1e-5;
	T ans = 1.0;
	while(ans > eps)
	{
		getMatrixEqual(Z, n, n, temp);
		getMatrixInverse(Z, n, invZ);
		getMatrixInverse(Y, n, invY);
		for(i = 0; i < n; i++)
			for(j = 0; j < n; j++)
			{
				Y[i][j] = (Y[i][j] + invZ[i][j]) / 2;
				Z[i][j] = (Z[i][j] + invY[i][j]) / 2;
			}
	    getMatrixSubstract(Z, temp, n, n, temp1);
		ans = fabs(getMatrixDet(temp1, n));
	}

	getMatrixEqual(Y, n, n, dst1);
	getMatrixEqual(Z, n, n, dst2);

	for(i = 0; i < n; i++)
	{
		delete[] Y[i];
		delete[] Z[i];
		delete[] invY[i];
		delete[] invZ[i];
		delete[] temp[i];
		delete[] temp1[i];
	}
	delete[] Y;
	delete[] Z;
	delete[] invY;
	delete[] invZ;
	delete[] temp;
	delete[] temp1;
}



/*************************************************************************************
* Function
*	getMatrixSubstract(): to calculate difference matrix between src1 and src2
*
* Parameters
*	T** src1: input matrix, m*n
*	T** src2: input matrix, m*n
*	int m: rows of matrices
* 	int n: columns of matrices
*	T** dst: output matrix, difference matrix between src1 and src2
*
* Return
*	void
*************************************************************************************/
template<typename T>
void getMatrixSubstract(T** src1, T** src2, int m, int n, T** dst)
{
	for(int i = 0; i < m; i++)
		for(int j =0; j < n; j++)
			dst[i][j] = src1[i][j] - src2[i][j];
}



/*************************************************************************************
* Function
*	getMatrixAdd(): to calculate summation matrix between src1 and src2
*
* Parameters
*	T** src1: input matrix
*	T** src2: input matrix
*	int m: rows of matrices
*	int n: columns of matrices
*	T** dst: output matrix, summation matrix between src1 and src2
* 
* Return
*	void
*************************************************************************************/
template<typename T>
void getMatrixAdd(T** src1, T** src2, int m, int n, T** dst)
{
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			dst[i][j] = src1[i][j] + src2[i][j];
}



/*************************************************************************************
* Function
*	getMatrixEqual(): to assign values from src to dst
*
* Parameters
*	T** src: input matrix
*	int m: rows of matrices
*	int n: columns of matrices
*	T** dst: output matrix, equal to src
*
* Return
*	void
*************************************************************************************/
template<typename T>
void getMatrixEqual(T** src, int m, int n, T** dst)
{
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			dst[i][j] = src[i][j];
}


/*************************************************************************************
* Function
*	getMatrixMultiply(): to calculate product matrix between src1 and src2
*
* Parameters
*	T** src1: input matrix
*	T** src2: input matrix
*	int m: rows of src1
*	int n: columns of src1 as well as rows of src2
*	int p: columns of src2
* 	T** dst: output matrix, product matrix between src1 and src2
*
* Return
*	void  
*************************************************************************************/
template<typename T>
void getMatrixMultiply(T** src1, T** src2, int m, int n, int p, T** dst)
{
	T* row = new T[n];
	T* column = new T[n];
	for(int i = 0; i < m; i++)
		for(int j = 0; j < p; j++)
		{
			getMatrixRow(src1, m, n, i, row);
			getMatrixColumn(src2, n, p, j, column);
			dst[i][j] = getVectorMultiply(row, column, n);
		}
	delete[] row;
	delete[] column;
}



/*************************************************************************************
* Function
*	getMatrixScale(): to calculate scaled matrix with factor scale 
*
* Parameters
*	T** src: input matrix
*	int m: rows of src
*	int n: columns of src
*	double scale: scale factor
* 	T** dst: output matrix, scaled matrix of src
*
* Return
*	void 
*************************************************************************************/
template<typename T>
void getMatrixScale(T** src, int m, int n, double scale, T** dst)
{
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			dst[i][j] = static_cast<T>(scale * src[i][j]);
}



/*************************************************************************************
* Function
*	getVectorMultiply(): to calculate vector-product between vec1 and vec2  
*
* Parameters
*	T* vec1: input vector
*	T* vec2: input vector
*	int n: dimension of vectors
*
* Return
*	T: product value of two vectors
*************************************************************************************/
template<typename T>
T getVectorMultiply(T* vec1, T* vec2, int n)
{
	T ans = 0;
	for(int i = 0; i < n; i++)
		ans += vec1[i] * vec2[i];
	return ans;
}


/*************************************************************************************
* Function
*	getVectorScale(): to calculate scaled vector with factor scale 
*
* Parameters
*	T* vec: input vector
*	int n: dimension of input vector
*	double scale: scale factor
* 	T* dst: output vector, scaled vector
*
* Return
*	void 
*************************************************************************************/
template<typename T>
void getVectorScale(T* vec, int n, double scale, T* dst)
{
	for(int i = 0; i < n; i++)
		dst[i] = static_cast<T>(scale * vec[i]);
}



/*************************************************************************************
* Function
*	getMatrixRow(): to get a row of src
*
* Parameters
*	T** src: input matrix
*	int m: rows of src
*	int n: columns of src
*	int i: index of the row extracted
*	T* irow: the ith row extracted from src
*
* Return
* 	void
*************************************************************************************/
template<typename T>
void getMatrixRow(T** src, int m, int n, int i, T* irow)
{
	for(int j = 0; j < n; j++)
		irow[j] = src[i][j];
}


/************************************************************************************
* Function
*	getMatrixColumn(): to get a column of src
*
* Parameters
*	T** src: input matrix
*	int m: rows of src
*	int n: columns of src
*	int i: index of the row extracted
*	T* icolumn: the ith column extracted from src
************************************************************************************/
template<typename T>
void getMatrixColumn(T** src, int m, int n, int i, T* icolumn)
{
	for(int j = 0; j < m; j++)
		icolumn[j] = src[j][i];
}


/***********************************************************************************
* Function
*	getVectorCrossMultiply(): to calculate cross product of two vectors
*
* Parameters
*	T* vec1: input vector, 3*1
*	T* vec2: input vector, 3*1
*	T* dst: output vector, 3*1
*
* Return
*	void
************************************************************************************/
template<typename T>
void getVectorCrossMultiply(T* vec1, T* vec2, T* dst)
{
	dst[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	dst[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	dst[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}


/***********************************************************************************
* Function
*	getVectorNorm(): to calculate norm of a vector
*
* Parameters
*	T* vec: input vector
*	int n: dimension of vector
*
* Return
*	double: norm of input vector
***********************************************************************************/
template<typename T>
double getVectorNorm(T* vec, int n)
{
	double elemSquare = 0;
	for(int i = 0; i < n; i++)
		elemSquare += vec[i] * vec[i];
	return sqrt(elemSquare);
}



/**********************************************************************************
* Function
*	getMatrixInterp2(): to interpolate points in src1
*
* Parameters
*	T** src1: input matrix
*	typePoint** src2: points coordinates interpolated
*	int m: rows of src1
*	int n: columns of src1
*	T** dst: output matrix in interpolated points
*
* Return
*	void
**********************************************************************************/
template<typename T>
void getMatrixInterp2(T** src1, typePoint2D<T>** src2, int m, int n, T** dst)
{
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < n; j++)
		{
			typePoint2D<T> uv = src2[i][j];
			T u = uv.x;
			T v = uv.y;
			int u1 = floor(u);
			int u2 = ceil(u);
			int v1 = floor(v);
			int v2 = ceil(v);
			if(u1 >= m)
			{
				u1 = m - 1; 
				u2 = m - 1;
			}
			if(u2 >= m)
				u2 = m - 1;
			if(v1 >= n)
			{
				v1 = n - 1;
				v2 = n - 1;
			}
			if(v2 >= n)
				v2 = n - 1;
			if(u2 < 0)
			{
				u1 = 0;
				u2 = 0;
			}
			if(u1 < 0)
				u1 = 0;
			if(v2 < 0)
			{
				v1 = 0; 
				v2 = 0;
			}
			if(v1 < 0)
				v1 = 0;
			if(u1 == u2 || v1 == v2)
				dst[i][j] = (src1[u1][v1] + src1[u1][v2] 
						+ src1[u2][v1] + src1[u2][v2]) / 4.0;
			else
				dst[i][j] = src1[u1][v1]*(u2 - u)*(v2 - v)
					+ src1[u2][v1]*(u - u1)*(v2 - v)
					+ src1[u1][v2]*(u2 - u)*(v - v1)
					+ src1[u2][v2]*(u - u1)*(v-v1);
		}
	}
}


/**********************************************************************************
* Function
*	quickSort(): sort an array from minimum to maximum
* 
* Parameters
*	T* array: an array with any data type
*	int nLow: beginning position while sorting
*	int nHigh: ending position while sorting
*
* Return
*	void
**********************************************************************************/
template<typename T>
void quickSort(T* array, int nLow, int nHigh)
{
	if(nLow < nHigh)
	{
		int i = nLow, j = nHigh;
		T temp = array[nLow];
		while(i < j)
		{
			while(i < j && array[j] >= temp)
				j--;
			if(i < j)
				array[i++] = array[j];

			while(i < j && array[i] < temp)
				i++;
			if(i < j)
				array[j--] = array[i];
		}
		array[i]= temp;
		quickSort(array, nLow, i - 1);
		quickSort(array, i + 1, nHigh);
	}
}


/**********************************************************************************
* Function
*	arrayMax(): find the maximum of an array
*
* Parameters
*	T* array: an array with any data type
*	int n: length of the array
*
* Return
*	T: maximum of the array
**********************************************************************************/
template<typename T>
T arrayMax(T* array, int n)
{
	T max = array[0];
	for(int i = 1; i < n; i++)
		if(array[i] > max)
			max = array[i];
	
	return max;
}


/**********************************************************************************
* Function
*	arrayMin(): find the minimum of an array
*
* Parameters
*	T* array: an array with any data type
*	int n: length of the array
*
* Return
*	T: minimum of the array
**********************************************************************************/
template<typename T>
T arrayMin(T* array, int n)
{
	T min = array[0];
	for(int i = 1; i < n; i++)
		if(array[i] < min)
			min = array[i];

	return min;
}


#endif





















