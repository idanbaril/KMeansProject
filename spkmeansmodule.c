#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "spkmeans.h"
  /*
   * for each spkmeans.c function, parsing the arguments from python to C, running the function and parsing back the results to python
   */
static PyObject* weightedAdjMatrix_capi(PyObject *self, PyObject *args){
    double **x ,**w;
    int n, vec_len, i, j;
    PyObject *pyx, *sublistX;
    Py_ssize_t t,t2, len;
    PyObject *result;
    PyObject *item;
    if (!PyArg_ParseTuple(args, "O!ii", &PyList_Type, &pyx, &n, &vec_len))
    {
        return NULL;
    }
    x = (double**)calloc(n, sizeof(double*));
    if(x==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        x[i] = (double*)calloc(vec_len, sizeof(double));
        if(x[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        sublistX = PyList_GetItem(pyx, i);
        for(j=0; j<vec_len;j++){
            x[i][j] = PyFloat_AsDouble(PyList_GetItem(sublistX, j));
        }
    }
    w = weightedAdjMatrix(x, n, vec_len);
    len = n;
    /*
    * convert back from C array to python list
    */
    result = PyList_New(len);
    for (t = 0; t < len; t++) {
        item = PyList_New(len);
        for (t2 = 0; t2 < len; t2++)
            PyList_SET_ITEM(item, t2, PyFloat_FromDouble(w[t][t2]));
        PyList_SET_ITEM(result, t, item);
    }
    free(x);
    free(w);
    return result;
}

static PyObject* diagonalDegreeMatrix_capi(PyObject *self, PyObject *args){
    double **d, **w;
    int n, i, j;
    PyObject *pyw, *sublistW;
    Py_ssize_t t,t2, len;
    PyObject *result;
    PyObject *item;
    if (!PyArg_ParseTuple(args, "O!i", &PyList_Type, &pyw, &n))
    {
        return NULL;
    }

    w = (double**)calloc(n, sizeof(double*));
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
        sublistW = PyList_GetItem(pyw, i);
        for(j=0; j<n;j++){
            w[i][j] = PyFloat_AsDouble(PyList_GetItem(sublistW, j));
        }
    }
    d = diagonalDegreeMatrix(w, n);
    len = n;
    /*
    * convert back from C array to python list
    */
    result = PyList_New(len);
    for (t = 0; t < len; t++) {
        item = PyList_New(len);
        for (t2 = 0; t2 < len; t2++)
            PyList_SET_ITEM(item, t2, PyFloat_FromDouble(d[t][t2]));
        PyList_SET_ITEM(result, t, item);
    }
    free(d);
    free(w);
    return result;
}

static PyObject* NormalizedGraphLaplacian_capi(PyObject *self, PyObject *args){
    double **d, **w, **l;
    int n, i, j;
    PyObject *pyd, *pyw, *sublistD, *sublistW;
    Py_ssize_t t,t2, len;
    PyObject *result;
    PyObject *item;
    if (!PyArg_ParseTuple(args, "O!O!i", &PyList_Type, &pyd, &PyList_Type, &pyw, &n))
    {
        return NULL;
    }
    d = (double**)calloc(n, sizeof(double*));
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
    w = (double**)calloc(n, sizeof(double*));
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
        sublistD = PyList_GetItem(pyd, i);
        sublistW = PyList_GetItem(pyw, i);
        for(j=0; j<n;j++){
            d[i][j] = PyFloat_AsDouble(PyList_GetItem(sublistD, j));
            w[i][j] = PyFloat_AsDouble(PyList_GetItem(sublistW, j));
        }
    }
    l = NormalizedGraphLaplacian(d, w, n);
    len = n;
    /*
    * convert back from C array to python list
    */
    result = PyList_New(len);
    for (t = 0; t < len; t++) {
        item = PyList_New(len);
        for (t2 = 0; t2 < len; t2++)
            PyList_SET_ITEM(item, t2, PyFloat_FromDouble(l[t][t2]));
        PyList_SET_ITEM(result, t, item);
    }
    free(d);
    free(w);
    free(l);
    return result;
}

static PyObject* jacobiAlgorithm_capi(PyObject *self, PyObject *args){
    double **a;
    int n, i, j;
    struct eigen *eigens;
    PyObject *pya, *sublistA;
    Py_ssize_t t,t2, len;
    PyObject *result, *resultVecs, *resultVals;
    PyObject *item;
    if (!PyArg_ParseTuple(args, "O!i", &PyList_Type, &pya, &n))
    {
        return NULL;
    }

    a = (double**)calloc(n, sizeof(double*));
    if(a==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        a[i] = (double*)calloc(n, sizeof(double));
        if(a[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<n; i++){
        sublistA = PyList_GetItem(pya, i);
        for(j=0; j<n;j++){
            a[i][j] = PyFloat_AsDouble(PyList_GetItem(sublistA, j));
        }
    }
    eigens = jacobiAlgorithm(a, n);
    len = n;
    /*
    * convert back from C array to python list
    */
    resultVals = PyList_New(len);
    resultVecs = PyList_New(len);
    result = PyList_New(2);
    for (t = 0; t < len; t++) {
        PyList_SET_ITEM(resultVals, t, PyFloat_FromDouble(eigens[t].value));
        item = PyList_New(len);
        for (t2 = 0; t2 < len; t2++){
            PyList_SET_ITEM(item, t2, PyFloat_FromDouble(eigens[t].vector[t2]));
        }
        PyList_SET_ITEM(resultVecs, t, item);
    }
    PyList_SET_ITEM(result, 0, resultVals);
    PyList_SET_ITEM(result, 1, resultVecs);
    free(a);
    free(eigens);
    return result;
    
}

static PyObject* determineK_capi(PyObject *self, PyObject *args){
    double *x;
    int n, i;
    PyObject *pyx, *k;
    if (!PyArg_ParseTuple(args, "O!i", &PyList_Type, &pyx, &n))
    {
        return NULL;
    }

    x = (double*)calloc(n, sizeof(double));
    if(x==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        x[i] = PyFloat_AsDouble(PyList_GetItem(pyx, i));
    }

    k = Py_BuildValue("i", determineK(x, n));
    free(x);
    return k;
}

static PyObject* formT_capi(PyObject *self, PyObject *args){
    double **vecs, **t;
    int n, i, j, k;
    PyObject *pyv, *sublistV;
    Py_ssize_t t1,t2, len, innerLen;
    PyObject *result;
    PyObject *item;
    if (!PyArg_ParseTuple(args, "iiO!", &k, &n, &PyList_Type, &pyv))
    {
        return NULL;
    }
/*    t = (double**)calloc(n, sizeof(double*));
    if(t==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        t[i] = (double*)calloc(k, sizeof(double));
        if(t[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    */
    vecs = (double**)calloc(k, sizeof(double*));
    if(vecs==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<k; i++){
        vecs[i] = (double*)calloc(n, sizeof(double));
        if(vecs[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<k; i++){
        sublistV = PyList_GetItem(pyv, i);
        for(j=0; j<n;j++){
            vecs[i][j] = PyFloat_AsDouble(PyList_GetItem(sublistV, j));
        }
    }
    t = formT(k, n, vecs);
    len = n;
    innerLen = k;
    /*
    * convert back from C array to python list
    */
    result = PyList_New(len);
    for (t1 = 0; t1 < len; t1++) {
        item = PyList_New(innerLen);
        for (t2 = 0; t2 < innerLen; t2++)
            PyList_SET_ITEM(item, t2, PyFloat_FromDouble(t[t1][t2]));
        PyList_SET_ITEM(result, t1, item);
    }
    free(t);
    free(vecs);
    return result;
}

static PyObject* kmeans_capi(PyObject *self, PyObject *args){
    int k, max_iter, n, vec_len,i,j;
    PyObject *pyx;
    PyObject *pyclusters;
    double e;
    double **x;
    double **clusters;
    double** res_clusters;
    PyObject* sublist;
    Py_ssize_t t,t2, len, inner_len;
    PyObject *result;
    PyObject *item;
    if (!PyArg_ParseTuple(args, "iiiiO!O!d",&k,&max_iter,&n,&vec_len,&PyList_Type,&pyx,&PyList_Type,&pyclusters,&e))
    {
        return NULL;
    }
    x = (double**)calloc(n, sizeof(double*));
    if(x==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(i=0; i<n; i++){
        x[i] = (double*)calloc(vec_len, sizeof(double));
        if(x[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    /*
    * convert objects from python list to C array
    */
    for(i=0; i<n; i++){
        sublist = PyList_GetItem(pyx, i);
        for(j=0; j<vec_len;j++){
            x[i][j] = PyFloat_AsDouble(PyList_GetItem(sublist, j));
        }
    }
    clusters = (double**)calloc(k, sizeof(double*));
    if(clusters==NULL){
            printf("An Error Has Occurred");
            exit(1);
    }
    for(i=0; i<k; i++){
        clusters[i] = (double*)calloc(vec_len, sizeof(double));
        if(clusters[i]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    for(i=0; i<k; i++){
        sublist = PyList_GetItem(pyclusters, i);
        for(j=0; j<vec_len;j++){
            clusters[i][j] = PyFloat_AsDouble(PyList_GetItem(sublist, j));
        }
    }
    res_clusters = kMeansAlgo(k,max_iter,n,vec_len,x,clusters,e);
    len = k;
    inner_len = vec_len;
    result = PyList_New(len);
    /*
    * convert back from C array to python list
    */
    for (t = 0; t < len; t++) {
        item = PyList_New(inner_len);
        for (t2 = 0; t2 < inner_len; t2++)
            PyList_SET_ITEM(item, t2, PyFloat_FromDouble(res_clusters[t][t2]));
        PyList_SET_ITEM(result, t, item);
    }
    free(x);
    free(clusters);
    return result;
}

static PyMethodDef capiMethods[] = {
        {"weighted_adj_matrix", (PyCFunction) weightedAdjMatrix_capi, METH_VARARGS},
        {"kmeans_algo", (PyCFunction) kmeans_capi, METH_VARARGS},
        {"diagonal_degree_matrix", (PyCFunction) diagonalDegreeMatrix_capi, METH_VARARGS},
        {"normalized_graph_laplacian", (PyCFunction) NormalizedGraphLaplacian_capi, METH_VARARGS},
        {"jacobi_algorithm", (PyCFunction) jacobiAlgorithm_capi, METH_VARARGS},
        {"form_T", (PyCFunction) formT_capi, METH_VARARGS},
        {"determine_K", (PyCFunction) determineK_capi, METH_VARARGS},
        {NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        "doc",
        -1,
        capiMethods

};

PyMODINIT_FUNC PyInit_spkmeansmodule(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if(!m){return NULL;}
    return m;
}

