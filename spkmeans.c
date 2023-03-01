#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


double distanceVectors(double *v1, double *v2, int arr_size){
    /*
   * computes the auclidian distance between two vectors
   */
    double sum = 0;
    int i;
    for(i=0; i < arr_size; i++){
        sum += pow(v1[i] - v2[i], 2);
    }
    return pow(sum, 0.5);
}

double clusterMean(int num_of_vecs, double **vectors, int coor){
    /*
   * computes vectors of a given cluster mean
   */
    double res;
    int i;
    res=0;
    for(i=0; i < num_of_vecs; i++){
        res += vectors[i][coor];
    }
    return res / num_of_vecs;
}

double** weightedAdjMatrix(double** x, int n, int vec_len){
    /*
   * gets n datapoints and returns their weighted adjancy matrix
   */
    int i, j;
    double** w = (double**)calloc(n, sizeof(double*));
    if(w==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        w[i] = (double*)calloc(n, sizeof(double));
        if(w[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i==j){
                w[i][j] = 0;
            }
            else{
                w[i][j] = exp((-1 * distanceVectors(x[i], x[j], vec_len))/2);
            }

        }
    }
    return w;
}

double** diagonalDegreeMatrix(double** w, int n){
  /*
   * gets weighted adjancy matrix and returns diagonal degree matrix
   */
    int i, j;
    double rowsum;
    double** d = (double**)calloc(n, sizeof(double*));
    if(d==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        d[i] = (double*)calloc(n, sizeof(double));
        if(d[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        rowsum = 0;
        for(j=0; j<n; j++){
            d[i][j] = 0;
            rowsum += w[i][j];
        }
        d[i][i] = rowsum;
    }
    return d;
}

void transformD(double** d, int n){
/*
* transform D to pow(D, -0.5)
*/
    int i;
    for(i=0; i<n; i++){
        d[i][i] = 1/(pow(d[i][i], 0.5));
    }
}

double** multiplyMatrices(double** m1 , double** m2, int r1, int c1, int r2, int c2){
  /*
   * gets two matrices A,B from size r1xc1 and r2xc2 and returns A*B from size r1xc2 (if c1=r2)
   */
    int i, j, k;
    double** res;
    if(r2 != c1){
        printf("r2 != c1");
        exit(1);
    }
    res = (double**)calloc(r1, sizeof(double*));
    if(res==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<r1; i++){
        res[i] = (double*)calloc(c2, sizeof(double));
        if(res[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }

    for (i = 0; i < r1; ++i) {
        for (j = 0; j < c2; ++j) {
            res[i][j] = 0;
        }
    }

    /*
    *  Multiplying first and second matrices and storing it in result
    */   
    for (i = 0; i < r1; ++i) {
        for (j = 0; j < c2; ++j) {
            for (k = 0; k < c1; ++k) {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return res;
}

double** NormalizedGraphLaplacian(double** d, double** w, int n){
  /*
   * gets diagonal degree matrix and weighted adjancy matrix, and returns their normalized graph laplacian
   */
    int i, j;
    double **multdw, **l;
    transformD(d, n);
    multdw = multiplyMatrices(d, w, n, n, n, n);
    l = multiplyMatrices(multdw, d, n, n, n, n);
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            l[i][j] = -1 * l[i][j];
            if(i==j){
                l[i][j] += 1;
            }
        }
    }
    free(multdw);
    return l;
}

struct eigen{
  /*
   * struct to represent eigenvalue and his eigenvector
   */
    double value;
    double* vector;
};

int compareEigenValues(const void *p1, const void *p2){
  /*
   * comparator function to compare between two instances of struct eigen
   */
    const struct eigen *s1 = p1;
    const struct eigen *s2 = p2;
    if(s1->value > s2->value){
        return 1;
    }
    if(s1->value < s2->value){
        return -1;
    }
    return 0;
}

double** transpose(double** m, int n){
  /*
   * gets matrix m from size nxn and returns m transpose
   */
    int i, j;
    double** mt = (double**)calloc(n, sizeof(double*));
    if(mt==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        mt[i] = (double*)calloc(n, sizeof(double));
        if(mt[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            mt[j][i] = m[i][j];
        }
    }
    return mt;
}

double** buildP(double** a, int n){
  /*
   * gets symetric matrix a and returns its rotation matrix p, for each iteration of jacobi algorithm
   */
    int i, j, imax=0, jmax=0, tmp, first=1;
    double teta, sign=1, t, c, s, maxval;
    double **p = (double**)calloc(n, sizeof(double*));
    if(p==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        p[i] = (double*)calloc(n, sizeof(double));
        if(p[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i!=j && (first == 1 || fabs(a[i][j]) > maxval)){
                jmax = j;
                imax = i;
                maxval = fabs(a[i][j]);
                first = 0;
            }
        }
    }

    if(jmax < imax){
        tmp = jmax;
        jmax = imax;
        imax = tmp;
    }

    teta = (a[jmax][jmax] - a[imax][imax])/(2 * a[imax][jmax]);
    if(teta < 0){
        sign = -1;
    }
    t = sign / (fabs(teta) + pow(pow(teta, 2) + 1, 0.5));
    c = 1 / pow(pow(t, 2) + 1, 0.5);
    s = t * c;

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i==j){
                p[i][j] = 1;
            }
            else{
                p[i][j] = 0;
            }
        }
    }
    p[imax][imax] = c;
    p[jmax][jmax] = c;
    p[imax][jmax] = s;
    p[jmax][imax] = -1*s;
    return p;
}

int isDiagonal(double** a, double** aTag, int n){
  /*
   * gets symetric matrix a and checks if it's diagonal (by our definition of diagonal)
   */
    int i, j;
    double aSum=0, aTagSum=0, e = exp(-15);
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i !=j ){
                aSum += pow(a[i][j], 2);
                aTagSum += pow(aTag[i][j], 2);
            }
        }
    }
    if (aSum - aTagSum <= e){
        return 1;
    }
    return 0;
}

struct eigen* jacobiAlgorithm(double** a, int n){
  /*
   * gets symetric matrix a and returns an array of struct eigen with a's n eigenvalues and their correspond eigenvectors
   */
    double **p, **v, **ptMultA, **tmpMat, **aTag, **pt;
    struct eigen *result;
    int first=1, maxIter = 100, countIter=0, i, j;
    double** copyA = (double**)calloc(n, sizeof(double*));
    if(copyA==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    /*
   * copy a to not destroy the original matrix
   */
    for(i=0; i<n; i++){
        copyA[i] = (double*)calloc(n, sizeof(double));
        if(copyA[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            copyA[i][j] = a[i][j];
        }
    }
    aTag = copyA;
    do{
        countIter++;
        p = buildP(aTag, n);
        pt = transpose(p, n);
        ptMultA = multiplyMatrices(pt, aTag, n, n, n, n);
        /*
       * not first iteration
       */
        if(first==0){
            free(copyA);
            copyA = aTag;
        }
        aTag = multiplyMatrices(ptMultA, p, n, n, n, n);
        free(ptMultA);
        free(pt);
        if(first==1){
            v = p;
            first = 0;
        }
        else{
            tmpMat = multiplyMatrices(v, p, n, n, n, n);
            for(i=0; i<n; i++){
                for(j=0; j<n; j++){
                    v[i][j] = tmpMat[i][j];
                }
            }
            free(tmpMat);
            free(p);
        }
    }while(isDiagonal(copyA, aTag, n) == 0 || countIter < maxIter);
    free(copyA);
    result = (struct eigen*)calloc(n, sizeof(struct eigen));
    if(result==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++) {
        result[i].vector = (double*)calloc(n, sizeof(double));
        if (result[i].vector == NULL) {
            printf("An Error Has Occurred");
            exit(1);
        }
    }

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            result[i].vector[j] = v[j][i];
        }
        result[i].value = aTag[i][i];
    }
    free(v);
    /*
   * sort the array of struct eigen by the eigenvalues
   */
    qsort(result, n, sizeof(struct eigen), compareEigenValues);
    if(aTag != NULL){
      free(aTag);
    }
    return result;
}


int determineK(double *eigenVals, int n){
  /*
   * gets n eigenvalues and determins k for the spk algorithm (in case k=0 in input)
   */
    double maxdiff = -1;
    double delta;
    int k=0, i;
    for(i=0; i<(int)n/2; i++){
        delta = fabs(eigenVals[i] - eigenVals[i+1]);
        if(delta > maxdiff){
            maxdiff = delta;
            k = i + 1;
        }
    }
    return k;
}

double** formT(int k, int n, double **eigenVecs){
  /*
   * gets matrix V of n eigenvectors, and forms matrix T according to step 5 in the algorithm
   */
    int i,j;
    double **t;
    double *rowSums = (double*)calloc(n,sizeof(double));
    if(rowSums==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    t = (double**)calloc(n, sizeof(double*));
    if(t==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        t[i] = (double*)calloc(k, sizeof(double));
        if(t[i]==NULL) {
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        rowSums[i] = 0;
    }
    for(i=0; i<n; i++){
        for(j=0; j<k;j++){
            t[i][j] = eigenVecs[j][i];
            rowSums[i] += pow(t[i][j], 2);
        }
    }
    for(i=0; i<n; i++){
        for(j=0; j<k;j++){
            t[i][j] = t[i][j] / pow(rowSums[i], 0.5);
        }
    }
    free(rowSums);
    return t;
}

double** kMeansAlgo(int k, int max_iter, int n, int vec_len ,double** x, double** clusters, double e){
  /*
   * gets n datapoints and clusters them into k clusters
   */
    int i,j,l, clus, counter, min_index;
    int *num_vecs_in_cluster;
    double maxdiff, **copyClusters, ***clusterObjects, *distanceList;
    double **clusterObjects_slice, dis;
    counter = 0;
    maxdiff = e + 1;
    while(maxdiff > e && counter <= max_iter){
        counter += 1;
        num_vecs_in_cluster = calloc(k, sizeof(int));
        if(num_vecs_in_cluster==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
        for(clus=0; clus < k; clus++){
            num_vecs_in_cluster[clus] = 0;
        }

        clusterObjects = calloc(k, sizeof(double**));
        if(clusterObjects==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
        for(i=0; i<k;i++){
            clusterObjects[i] = calloc(n,sizeof(double*));
            if(clusterObjects[i]==NULL){
                printf("An Error Has Occurred");
                exit(1);
            }
            for(j=0; j<n;j++){
                clusterObjects[i][j] = calloc(vec_len,sizeof(double));
                if(clusterObjects[i][j]==NULL){
                    printf("An Error Has Occurred");
                    exit(1);
                }
            }
        }
        for(i=0; i < n; i++){
            distanceList = calloc(k, sizeof(double));
            if(distanceList==NULL){
                printf("An Error Has Occurred");
                exit(1);
            }
            /*
            * saves all diffs from clusters to check which one is closest
            */
            for(j=0; j < k; j++){
                distanceList[j] = distanceVectors(x[i], clusters[j], vec_len);
            }

            min_index = 0;
            for(j=0; j < k; j++){
                if(distanceList[j] < distanceList[min_index]){
                    min_index = j;
                }
            }
            for(j=0; j < vec_len; j++){
                clusterObjects[min_index][num_vecs_in_cluster[min_index]][j] = x[i][j];
            }
            num_vecs_in_cluster[min_index]++;
            free(distanceList);
        }

        copyClusters = calloc(k, sizeof(double*));
        if(copyClusters==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
        for(i=0; i<k;i++){
            copyClusters[i] = calloc(vec_len, sizeof(double));
            if(copyClusters[i]==NULL){
                printf("An Error Has Occurred");
                exit(1);
            }
        }
        /*
        * save former iteration to check if diff is smaller than e
        */
        for(i=0; i < k; i++){
            for(j=0; j < vec_len; j++){
                copyClusters[i][j] = clusters[i][j];
            }
        }

        for(i=0; i < k; i++){
            clusterObjects_slice = calloc(num_vecs_in_cluster[i], sizeof(double*));
            if(clusterObjects_slice==NULL){
                printf("An Error Has Occurred");
                exit(1);
            }
            for(j=0; j<num_vecs_in_cluster[i]; j++){
                clusterObjects_slice[j] = calloc(vec_len, sizeof(double));
                if(clusterObjects_slice[j]==NULL){
                    printf("An Error Has Occurred");
                    exit(1);
                }
            }
            /*
            * save cluster's vectors in a fitting size array
            */
            for(j=0; j < num_vecs_in_cluster[i]; j++){

                for(l=0; l < vec_len; l++){
                    clusterObjects_slice[j][l] = clusterObjects[i][j][l];
                }
            }
            for(j=0; j < vec_len; j++){
                clusters[i][j] = clusterMean(num_vecs_in_cluster[i], clusterObjects_slice, j);
            }
            free(clusterObjects_slice);
        }
        for(i=0; i < k; i++){
            dis = distanceVectors(clusters[i], copyClusters[i], vec_len);
            if(dis > maxdiff){
                maxdiff = dis;
            }
        }
        maxdiff = pow(maxdiff, 0.5);
        free(clusterObjects);
        free(copyClusters);
        free(num_vecs_in_cluster);
    }

    return clusters;
}

enum goal{wam, ddg, lnorm, jacobi};

void printMatrix(double** matrix, int n){
    int i, j;
    for(i=0; i<n; i++){
        printf("%.4f", matrix[i][0]);
        for(j=1; j<n; j++){
            printf(",%.4f", matrix[i][j]);
        }
        printf("\n");
    }
}

void printJacobi(struct eigen *eigens, int n){
  /*
   * prints n eigenvalues and then their n corresponding eigenvectors
   */
    int i, j;
    for(i=0; i<n; i++){
        if(eigens[i].value < 0 && eigens[i].value > -0.0001){
            eigens[i].value = 0.0000;
        }
    }
    printf("%.4f",eigens[0].value);
    for(i=1; i<n; i++){
        printf(",%.4f", eigens[i].value);
    }
    printf("\n");
    for(i=0; i<n; i++){
        printf("%.4f", eigens[0].vector[i]);
        for(j=1; j<n; j++){
            printf(",%.4f", eigens[j].vector[i]);
        }
        printf("\n");
    }
}

int main(int argc ,char* argv[]) {
    char *filename, c;
    double **x, coor, **w, **d, **l;
    struct eigen *eigens;
    enum goal g;
    int n, vec_len, i, j;
    FILE* fp;
    /*
   * input validation
   */
    if(argc != 3){
        printf("Invalid Input!");
        exit(1);
    }

    if(!strcmp(argv[1],"wam")){
        g = wam;
    }
    else if(!strcmp(argv[1],"ddg")){
        g = ddg;
    }
    else if(!strcmp(argv[1],"lnorm")){
        g = lnorm;
    }
    else if(!strcmp(argv[1],"jacobi")){
        g = jacobi;
    }
    else{
        printf("Invalid Input!");
        exit(1);
    }
  /*
   * parsing the datapoints file
   */
    filename = argv[2];
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        perror("An Error Has Occurred\n");
        exit(EXIT_FAILURE);
    }

    n = 0;

    while(!feof(fp)){
        if(fgetc(fp)=='\n'){
            n++;
        }
    }
    fclose(fp);
    fp = fopen(filename, "r");
    vec_len = 1;

    while((c = fgetc(fp)) != '\n'){
        if(c == ','){
            vec_len++;
        }
    }
    fclose(fp);
    x = calloc(n, sizeof(double*));
    if(x==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n;i++){
        x[i] = calloc(vec_len,sizeof(double));
        if(x[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    fp = fopen(filename, "r");

    for(i=0; i < n; i++){
        for(j=0; j < vec_len; j++){
            fscanf(fp, "%lf", &coor);

            if(!feof(fp)){
                fgetc(fp);
            }
            x[i][j] = coor;
        }
    }

    fclose(fp);
    /*
   * parsing goal and run the corresponding algorithm
   */
    if(g==jacobi){
        eigens = jacobiAlgorithm(x, n);
        printJacobi(eigens, n);
        free(eigens);
        return 0;
    }

    w = weightedAdjMatrix(x, n ,vec_len);
    if(g==wam){
        printMatrix(w, n);
        free(w);
        return 0;
    }
    d = diagonalDegreeMatrix(w, n);
    if(g==ddg){
        printMatrix(d, n);
        free(w);
        free(d);
        return 0;
    }
    l = NormalizedGraphLaplacian(d, w, n);
    if(g==lnorm){
        printMatrix(l, n);
        free(w);
        free(d);
        free(l);
        return 0;
    }

    return 0;
}