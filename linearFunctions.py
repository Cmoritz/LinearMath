#Chris Moritz Python Procedures
#Linear Algebra

#Functions for project 1
#Functions for project 2 are further down
def vplus(v,w): return [v[i] + w[i] for i in range(len(v))]

def smult(a,v): return [v[i]*a for i in range(len(v))]

def vzero(n): return [0 for i in range(n)]

def vneg(v): return [v[i]*-1 for i in range(len(v))]

def dot(v,w): return sum([v[i] * w[i] for i in range(len(v))])

def sbasis(j,n): return [1 if i == j - 1 else 0 for i in range(n)]

def vsum(vlist):
    v = [0]*len(vlist[0])
    for i in range(len(v)):
        for j in vlist:
            v[i] = j[i] + v[i]
    return v

def lincomb(clist, vlist):
    v = vlist
    for x in range(len(clist)):
        v[x] = smult(clist[x], v[x])
    v = vsum(v)
    return v

def mplus(A,B): return [[A[j][i] + B[j][i] for i in range(len(A))] for j in range(len(A[0]))]

def cmmult(c, A): return [[A[j][i] * c for i in range(len(A))] for j in range(len(A[0]))]

def mzero(m, n): return [[0 for i in range(n)] for j in range(m)]

def mneg(A): return [[A[j][i] * -1 for i in range(len(A))] for j in range(len(A[0]))]

def ID(n): return [[1 if i == j else 0 for i in range(n)] for j in range(n)]

def shape(A): return [len(A), len(A[0])]

def transpose(A): return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]

def mvmult(A, v):
    tempMatrix = A
    for i in range(len(A)):
        for j in range(len(A[i])):
            tempMatrix[i][j] = tempMatrix[i][j] * v[j]
    return [sum(i) for i in tempMatrix]

def mmult(A, B): return [[sum(a*b for a,b in zip(x,y)) for y in zip(*B)] for x in A]

def acompatible(A, B): return len(A) == len(B) and len(A[0]) == len(B[0])

def mcompatible(A,B): return len(A) == len(B[0]) and len(B) == len(A[0])

def mtov(A): return [A[i][j] for i in range(len(A)) for j in range(len(A[0]))]

def swap(j,k,A):
    tempVal = A[j-1]
    A[j-1] = A[k-1]
    A[k-1] = tempVal
    return A

def addrow(c,j,k,A):
    newMatrix = A
    newRow = smult(c, A[j-1])
    for i in range(len(A[k])):
        newMatrix[k-1][i] = A[k-1][i] + newRow[i]
    return newMatrix

def augment(A,B):
    tempMatrix = A
    for i in range(len(B)):
        for j in range(len(B[0])):
            tempMatrix[i].append(B[i][j])
    return tempMatrix






##Functions for Project 2
def norm(v):
    return dot(v, v)**(1/2)

def normc(v):
    return sum([(v[i] * (v[i].real + (v[i].imag*-1j))**1/2) for i in range(len(v))])

#Normalizes one vector
def normalize1(v):
    vecLength = norm(v)
    return [v[i] / vecLength for i in range(len(v))]

#Normalizes the set of vectors
def normalize(S):
    normedMatrix = S
    for i in range(len(S)):
        S[i] = normalize1(S[i])
    return normedMatrix

#Put the matrix in reduced row echelon form
def RRE(A):
    RRE = cleanMatrix(A)
    rowLen = len(A) #How many columns
    colLen = len(A[0]) #How many rows
    i = 0
    j = 0
    while(i < rowLen and j < colLen):
        isNonZero = 0
        #First, search the current column for the pivot
        pivotRow = 0
        pivotCol = 0
        for k in range(rowLen + i):
            if RRE[i][j] != 0:
                pivotRow = i
                pivotCol = j
                isNonZero = 1
                break
        if isNonZero is 1:
            #After the pivot has been found, normalize the row by dividing by the pivot
            pivotNum = RRE[pivotRow][j]
            for k in range(colLen):
                RRE[i][k] = RRE[i][k] / pivotNum
            #Now we must add multipes of the current row to all other values in the pivot column to make them 0
            l = i + 1
            while(l < rowLen):
                normNum = -1 * RRE[l][pivotCol]
                RRE[l] = vplus(smult(normNum, RRE[pivotRow]),RRE[l])
                l = l + 1
            #Now we must do the same so that above the pivot are 0 as well
            l = i - 1
            while(l >= 0):
                normNum = -1 * RRE[l][pivotCol]
                RRE[l] = vplus(smult(normNum, RRE[pivotRow]),RRE[l])
                l = l - 1
            i = i + 1
            j = j + 1
        else:
            #i = i + 1
            j = j + 1
    return RRE

#Find the orthonomal basis to a basis
def GS(B):
    orthBase = []
    firstBasis = normalize1(B[0])
    orthBase.append(firstBasis)
    i = 1
    while i < len(B):
        nextVec = B[i]
        for j in range(i):
            nextVec = normalize1(vsub(nextVec, smult(dot(nextVec, orthBase[j]), orthBase[j])))
        orthBase.append(nextVec)
        i = i + 1
    return orthBase

#Find the inverse of the given matrix
def inverse(A):
    #First need to append the identity matrix to the passed matrix
    if len(A) is len(A[0]):
        inverseMat = []
        Mat = augment(A, ID(len(A[0])))
        Mat = RRE(Mat)
        for i in range(len(Mat)):
            j = int(len(Mat[0])/2)
            c = 0
            row = []
            while j < len(Mat[0]):
                print(j)
                row.append(Mat[i][j])
                c = c + 1
                j = j + 1
            inverseMat.append(row)
        return inverseMat
    else:
        print("Non square matrices are not invertible")
        return [[]]


#Helper Functions
#Used to get rid of rows of 0
def cleanMatrix(A):
    #First scan the matrix for 0 rows
    cleaned = A
    seenZero = 1
    while seenZero is 1:
        seenZero = 0
        i = 0
        while seenZero is 0 and i < len(cleaned):
            hasZeroRow = 1
            for j in range(len(cleaned[0])):
                if cleaned[i][j] != 0:
                    hasZeroRow = 0
            if hasZeroRow is 1:
                seenZero = 1
                cleaned.pop(i)
            i = i + 1
    return cleaned

def vsub(a, b):
    return [a[i] - b[i] for i in range(len(a))]






























