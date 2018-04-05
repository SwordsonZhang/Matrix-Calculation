#ifndef CALCMATRIX_H_
#define CALCMATRIX_H_

/****************************************************************************************
* Struct: 
*	typePoint: point in coordinate(u, v)
* 
* Variables:
*	double u: x-coordinate  
*	double v: y-coordinate
****************************************************************************************/
struct typePoint
{
	double u;
	double v;
};


/***************************************************************************************
* Function
*	getMatrixInverse(): to calculate the inverse matrix
*
* Parameters
*	double** src: input matrix, square matrix
*	int n: dimension of input matrix
*	double** dst: output matrix, inverse square matrix
*
* Return
*	bool: inversed if true, or false if failed
***************************************************************************************/
bool getMatrixInverse(double** src, int n, double** dst);


/**************************************************************************************
* Function
*	getMatrixDet(): to calculate determinant of src
*
* Parameters
*	double** src: input matrix
*	int n: dimension of src
*
* Return
* 	double: value of the determinant
**************************************************************************************/
double getMatrixDet(double** src, int n);


/**************************************************************************************
* Function
*	getMatrixStar(): to calculate adjoint matrix of src
*
* Parameters
*	double** src: input square matrix
*	int n: dimension of src
*	double** dst: adjoint matrix of src
*
* Return
*	void
**************************************************************************************/
void getMatrixStar(double** src, int n, double** dst);


/*************************************************************************************
* Function
*	getMatrixSqrt(): to calculate square matrix of src
*
* Parameters
*	double** src: input matrix
*	int n: dimension of src
*	double** dst1: 1/2 of input matrix
*	double** dst2: -1/2 of input matrix
*
* Return
*	void
*************************************************************************************/
void getMatrixSqrt(double** src, int n, double** dst1, double** dst2);


/*************************************************************************************
* Function
*	getMatrixSubstract(): to calculate difference matrix between src1 and src2
*
* Parameters
*	double** src1: input matrix, m*n
*	double** src2: input matrix, m*n
*	int m: rows of matrices
* 	int n: columns of matrices
*	double** dst: output matrix, difference matrix between src1 and src2
*
* Return
*	void
*************************************************************************************/
void getMatrixSubstract(double** src1, double** src2, int m, int n, double** dst);


/*************************************************************************************
* Function
*	getMatrixAdd(): to calculate summation matrix between src1 and src2
*
* Parameters
*	double** src1: input matrix
*	double** src2: input matrix
*	int m: rows of matrices
*	int n: columns of matrices
*	double** dst: output matrix, summation matrix between src1 and src2
* 
* Return
*	void
*************************************************************************************/
void getMatrixAdd(double** src1, double** src2, int m, int n, double** dst);


/*************************************************************************************
* Function
*	getMatrixEqual(): to assign values from src to dst
*
* Parameters
*	double** src: input matrix
*	int m: rows of matrices
*	int n: columns of matrices
*	double** dst: output matrix, equal to src
*
* Return
*	void
*************************************************************************************/
void getMatrixEqual(double** src, int m, int n, double** dst);


/*************************************************************************************
* Function
*	getMatrixMultiply(): to calculate product matrix between src1 and src2
*
* Parameters
*	double** src1: input matrix
*	double** src2: input matrix
*	int m: rows of src1
*	int n: columns of src1 as well as rows of src2
*	int p: columns of src2
* 	double** dst: output matrix, product matrix between src1 and src2
*
* Return
*	void  
*************************************************************************************/
void getMatrixMultiply(double** src1, double** src2, int m, int n, int p, double** dst);


/*************************************************************************************
* Function
*	getMatrixScale(): to calculate scaled matrix with factor scale 
*
* Parameters
*	double** src: input matrix
*	int m: rows of src
*	int n: columns of src
*	double scale: scale factor
* 	double** dst: output matrix, scaled matrix of src
*
* Return
*	void 
*************************************************************************************/
void getMatrixScale(double** src, int m, int n, double scale, double** dst);


/*************************************************************************************
* Function
*	getVectorMultiply(): to calculate vector-product between vec1 and vec2  
*
* Parameters
*	double* vec1: input vector
*	double* vec2: input vector
*	int n: dimension of vectors
*
* Return
*	double: product value of two vectors
*************************************************************************************/
double getVectorMultiply(double* vec1, double* vec2, int n);


/*************************************************************************************
* Function
*	getVectorScale(): to calculate scaled vector with factor scale 
*
* Parameters
*	double* vec: input vector
*	int n: dimension of input vector
*	double scale: scale factor
* 	double* dst: output vector, scaled vector
*
* Return
*	void 
*************************************************************************************/
void getVectorScale(double* vec, int n, double scale, double* dst);


/*************************************************************************************
* Function
*	getMatrixRow(): to get a row of src
*
* Parameters
*	double** src: input matrix
*	int m: rows of src
*	int n: columns of src
*	int i: index of the row extracted
*	double* irow: the ith row extracted from src
*
* Return
* 	void
*************************************************************************************/
void getMatrixRow(double** src, int m, int n, int i, double* irow);


/************************************************************************************
* Function
*	getMatrixColumn(): to get a column of src
*
* Parameters
*	double** src: input matrix
*	int m: rows of src
*	int n: columns of src
*	int i: index of the row extracted
*	double* icolumn: the ith column extracted from src
************************************************************************************/
void getMatrixColumn(double** src, int m, int n, int i, double* icolumn);


/***********************************************************************************
* Function
*	getVectorCrossMultiply(): to calculate cross product of two vectors
*
* Parameters
*	double* vec1: input vector, 3*1
*	double* vec2: input vector, 3*1
*	double* dst: output vector, 3*1
*
* Return
*	void
************************************************************************************/
void getVectorCrossMultiply(double* vec1, double* vec2, double* dst);


/***********************************************************************************
* Function
*	getVectorNorm(): to calculate norm of a vector
*
* Parameters
*	double* vec: input vector
*	int n: dimension of vector
*
* Return
*	double: norm of input vector
***********************************************************************************/
double getVectorNorm(double* vec, int n);


/**********************************************************************************
* Function
*	getMatrixInterp2(): to interpolate points in src1
*
* Parameters
*	double** src1: input matrix
*	typePoint** src2: points coordinates interpolated
*	int m: rows of src1
*	int n: columns of src1
*	double** dst: output matrix in interpolated points
*
* Return
*	void
**********************************************************************************/
void getMatrixInterp2(double** src1, typePoint** src2, int m, int n, double** dst);

#endif
