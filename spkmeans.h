//
// Created by Administrator on 17/02/2022.
//

#ifndef SOFTWAREPROJECT_SPKMEANS_H
#define SOFTWAREPROJECT_SPKMEANS_H
double** weightedAdjMatrix(double** x, int n, int vec_len);
double** diagonalDegreeMatrix(double** w, int n);
double** NormalizedGraphLaplacian(double** d, double** w, int n);
struct eigen* jacobiAlgorithm(double** a, int n);
int determineK(double *eigenVals, int n);
double** kMeansAlgo(int k, int max_iter, int n, int vec_len ,double** x, double** clusters, double e);
double** formT(int k, int n, double **eigenVecs);
struct eigen{
    double value;
    double* vector;
};
#endif //SOFTWAREPROJECT_SPKMEANS_H

