#ifndef CALCMATRIX_H_
#define CALCMATRIX_H_

bool getMatrixInverse(double** src, int n, double** dst);
double getMatrixDet(double** src, int n);
void getMatrixStar(double** src, int n, double** dst);
void getMatrixSqrt(double** src, int n, double** dst1, double** dst2);
void getMatrixSubstract(double** src1, double** src2, int n, double** dst);
void getMatrixAdd(double** src1, double** src2, int n, double** dst);
void getMatrixEqual(double** src, int n, double** dst);
void getMatrixMultiply(double** src1, double** src2, int m, int n, int p, double** dst);
void getMatrixScale(double** src, int m, int n, double scale, double** dst);
double getVectorMultiply(double* vec1, double* vec2, int n);
void getVectorScale(double* vec, int n, double scale, double* dst);
void getMatrixRow(double** src, int m, int n, int i, double* irow);
void getMatrixColumn(double** src, int m, int n, int i, double* icolumn);
void getVectorCrossMultiply(double* vec1, double* vec2, double* dst);
double getVectorNorm(double* vec, int n);

#endif
