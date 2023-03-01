import sys
import numpy as np
import pandas as pd
import spkmeansmodule as sp
import enum

MAXITER = 300
EPS = 0

def distanceVectors(v1, v2):
    #compute auclidian distance between two vectors
    sum = 0
    for i in range(len(v1)):
        sum += (v1[i] - v2[i]) ** 2
    return sum


class Goal(enum.Enum):
    spk = "spk"
    wam = "wam"
    ddg = "ddg"
    lnorm = "lnorm"
    jacobi = "jacobi"


def main():
    # cmd arguments validation
    if len(sys.argv) != 4:
        print("Invalid Input!")
        quit()

    if not sys.argv[1].isdigit():
        print("Invalid Input!")
        quit()
    k = int(sys.argv[1])
    if(k==1):
        print("Invalid Input!")
        quit()
    
    if sys.argv[2] not in ["spk", "wam", "ddg", "lnorm", "jacobi"]:
        print("Invalid Input!")
        quit()
    goal = Goal(sys.argv[2])
    file_name = sys.argv[3]

    df = pd.read_csv(file_name, header=None)
    x = df.values.tolist()
    n = len(x)
    vecLen = len(x[0])
    
    if k >= len(x):
        print("Invalid Input!")
        quit()
        
    # parse enum value and choose the correspond algorithm
    
    if goal == Goal.jacobi:
        eigenValsAndVecs = sp.jacobi_algorithm(x, n)
        transVecs = transpose(eigenValsAndVecs[1])
        print_result(transVecs, indices=eigenValsAndVecs[0])
        quit()
        
    w = sp.weighted_adj_matrix(x, n, vecLen)
    if goal == Goal.wam:
        print_result(w)
        quit()
    
    d = sp.diagonal_degree_matrix(w, n)
    if goal == Goal.ddg:
        print_result(d)
        quit()
    
    l = sp.normalized_graph_laplacian(d, w, n)
    if goal == Goal.lnorm:
        print_result(l)
        quit()
    
    eigenValsAndVecs = sp.jacobi_algorithm(l, n)
        
    if goal == Goal.spk:        
        if(k==0):
            k = sp.determine_K(eigenValsAndVecs[0], n)
        t = sp.form_T(k, n, eigenValsAndVecs[1])
        initCentroids, indices = kmeanspp(t, k)
        clusters = sp.kmeans_algo(k, MAXITER, n, k, t, initCentroids, EPS)
        print_result(clusters, indices=indices)

# get matrix m and return m transpose
def transpose(matrix):
    matT = []
    for i in range(len(matrix[0])):
        matT.append([])
        for j in range(len(matrix)):
            matT[i].append(matrix[j][i])
    return matT

# print the results, with indices line on top, and then k clusters
# if goal = jacobi, prints eigenvalues line on top, and then n eigenvectors
def print_result(matrix, indices=[]):
    for index in indices:
        if "{:.4f}".format(index) == "-0.0000":
            index = 0.0000
    if indices != []:
        print("{:.4f}".format(indices[0]), end="")
        for i in range(1, len(indices)):
            print(",{:.4f}".format(indices[i]), end="")
        print("")
    for row in matrix:
        print("{:.4f}".format(row[0]), end="")
        for j in range(1, len(row)):
            print(",{:.4f}".format(row[j]), end="")
        print("")


def kmeanspp(x, k):
#initialize the centroids for the kmeans algorithm

    n = len(x)
    indices = []
    np.random.seed(0)
    i = 1
    clusters = []
    indexVec = np.random.choice(n)
    clusters.append(x[indexVec])
    indices.append(indexVec)
    numbers = [num for num in range(n)]
    while i < k:
        D = []
        for l in range(n):
            disClosestCluster = float("inf")
            for j in range(i):
                currdis = distanceVectors(x[l], clusters[j])
                disClosestCluster = min(disClosestCluster, currdis)
            D.append(disClosestCluster)
        P = []
        for l in range(n):
            P.append(D[l] / sum(D))
        i += 1
        indexVec = np.random.choice(numbers, p=P)
        clusters.append(x[indexVec])
        indices.append(indexVec)
    return clusters, indices


main()
