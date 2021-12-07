import numpy as np
from numpy.linalg import inv
from numpy import genfromtxt

def input_standardform():
    import numpy as np

    #Defining the A matrix
    m = input("Enter the number of Linearly Independent Constriants: ")
    n = input("Enter the number of Variables: ")
    z = input("Maximize(1) or Minimize(0): ")

    m = int(m)
    n = int(n)
    A = np.zeros((m,n))
    b = np.zeros((m,1))
    C = np.zeros((n,1))

    i = 0
    j = 0
    for i in range(m):
        for j in range(n):
            print("Enter the value of X",i+1,j+1)
            A[i,j] = input()
            j+=1
        i+=1

    #Defining the b matrix
    k=0
    for k in range(m):
        print("Enter the value of b",k+1)
        b[k,0] = input()
        k+=1

    #Defining the Cost Coefficients matrix
    k=0
    for k in range(n):
        print("Enter the cost coefficent of X",k+1)
        C[k,0] = input()
        k+=1
    
    return A, b, C, z

def showequations(A, b, C, z):
    m,n = np.shape(A)

    if z == 1:
        print("Maximize")
    else:
        print("Minimize")

    k=0
    for k in range(n):
        print(int(C[k]),"X", k+1, "+")
        k+=1

    print(" ")

    i = 0
    j = 0
    
    for i in range(m):
        for j in range(n):
            print(A[i,j], "X", j+1, "+")
            j+=1
        print("=",int(b[i]))
        i+=1

def findslack(A):
    m,n = np.shape(A)
    l = 0
    k = 0
    Xb = np.zeros((m))

    for l in range(m):
        for k in range(n):
            if A[l,n-k-1]==1:
                Xb[l] = n-k
                break
            k+=1
        l+=1
    return Xb

def basic_matrix(A,Xb):
    m,n = np.shape(A)
    B = np.zeros((m, m))
    i = 0
    j = 0
    for i in range(m):
        for j in range(m):
            B[i,j] = A[i, int(Xb[j])-1]
    return B

def BFS(B):
    Binv = inv(B)
    Binv_b = Binv @ b
    return Binv_b, Binv

def find_Cb(C,Xb,m):
    Cb = np.zeros((m,1))
    i = 0
    for i in range(m):
        Cb[i] = int(C[int(Xb[i]-1)])
    return Cb

def initiate_tableau(A, Binv_b):
    m,n = np.shape(A)    
    tableau = np.zeros((m+1,n+1))
    tableau[1:m+1, 0:n] = A
    tableau[1:m+1, n] = np.asarray(Binv_b[:,0])
    return tableau

def reduced_cost(Cb,A,C,tableau,m,n):
    tableau[0,0:n] = np.subtract((np.transpose(Cb) @ A), np.transpose(C))
    return tableau

def obj(tableau,Cb,m,n):
    tableau[0,n] = (np.transpose(Cb) @ tableau[1:m+1, n])
    return tableau

def add_artifical_variables(A,C,M):
    m,n = np.shape(A)
    for i in range(m):
        x_a = np.zeros((m,1))
        x_a[i,0] = 1
        A = np.concatenate((A, x_a), 1)
        c_a =  np.array([[M]])
        C = np.concatenate((C, c_a), 0)
    return A, C

A = genfromtxt('/Users/kaustubhrajimwale/Desktop/Fall 21/IE535/Computing Project/Standard Form Convertor/A.csv', delimiter=',')
A = np.delete(A, 0, 1)

b = genfromtxt('/Users/kaustubhrajimwale/Desktop/Fall 21/IE535/Computing Project/Standard Form Convertor/b.csv', delimiter=',')
b = np.delete(b, 0, 1)

C = genfromtxt('/Users/kaustubhrajimwale/Desktop/Fall 21/IE535/Computing Project/Standard Form Convertor/C.csv', delimiter=',')
C = np.delete(C, 0, 1)

z = input("\nMaximize(1) or Minimize(0): ")
z = int(z)

if z != 0:
    C = C * (-1)


m,n = np.shape(A)
A,C = add_artifical_variables(A,C,100)
m,n = np.shape(A)
Xb = findslack(A)                               # DEFINING ARTIFICIAL VARIABLES
B = basic_matrix(A, Xb)                         # DEFINING B MATRIX
Binv_b, Binv = BFS(B)                           # DEFINING BINV*B
Cb = find_Cb(C,Xb,m)                            # DEFINING COST COEFFICIENTS OF BASIC VARIABLES
tableau = initiate_tableau(A, Binv_b)           # CONSTRUCTING TABLEAU
tableau = reduced_cost(Cb,A,C,tableau,m,n)      # DEFINING REDUCED COST
tableau = obj(tableau,Cb,m,n)                   # DEFINING OBJECTIVE FUCTION



iter = 0
print("\n===============================")
print("STARTING TABLEAU")
print("===============================\n")
print(tableau)
print("\n")


flag = int(1)
unbdd = 0

while(flag!=0):
    print("\n===============================")
    print("ITERATION ",iter+1)
    print("===============================\n")
    iter += 1

    #CHECKING FOR OPTIMALITY
    count=0
    for i in range (n):
        if tableau[0,i]<=0:
            count += 1
    if(count==n):
        optimal = 1
        flag = 0
        break
    
    else:
        z_max = 0
        for i in range(n):
            if tableau[0,i] >= z_max:
                z_max = tableau[0,i]
                z_max_i = i
    pivot_col = int(z_max_i)


    #CHECKING FOR UNBOUNDEDNESS
    y = (tableau[1:m+1,pivot_col]<=0)
    if (sum(y)==m):
        flag = 0
        unbdd = 1
        break


    #RATIO TEST
    ratio = 1000
    i = 0
    for i in range(m):
        if(tableau[i+1,pivot_col]>0):
            if (tableau[i+1,n]/tableau[i+1,pivot_col])<=ratio:
                ratio = (tableau[i+1,n]/tableau[i+1,pivot_col])
                pivot_row = i+1
    print("\nPIVOT ELEMENT = ", tableau[pivot_row,pivot_col])
    print("\nENTERING VARIABLE = X",pivot_col+1)
    print("LEAVING VARIABLE = X",int(Xb[pivot_row-1]),"\n")
    Xb[pivot_row-1] = pivot_col+1

    #ROW TRANSFORMATIONS
    i=0
    tableau[pivot_row,:] = tableau[pivot_row,:]/tableau[pivot_row,pivot_col]
    for i in range(m+1):
        if (i) != pivot_row:
            tableau[i,:] = tableau[i,:] - tableau[pivot_row,:]*int(tableau[i,pivot_col])
        i+=1
    
    print(tableau, Xb)

logic = (Xb<=(n-m))

if sum(logic)>=m:
    if optimal == 1:
        print(" ")
        print("ITERATIONS STOPPED\n")
        print("The current BFS is optimal:\n")
        print("Objective Function =",tableau[0,n])
    
        j=0
        for j in range(m):
            print("X",int(Xb[j]),"=",(tableau[j+1,n]))
        print(" ")
        print("All other variables are non-basic and are equal to 0.")
        print(" ")
        print(" ")

    if unbdd == 1:
        print("THE LP IS UNBOUNDED ALONG THE DIRECTION:")
        dir = np.zeros((n,1))
        x_unbdd = np.zeros((n,1))
        i = 0
        for i in range(m):
            x_unbdd[int(Xb[i]-1),0] = tableau[i+1,n]
            dir[int(Xb[i]-1),0] = -1*tableau[i+1,pivot_col]
        dir[pivot_col,0] = 1
        
        i = 0
        print(" ")
        for i in range(n):
            print('x',i+1,' = ',x_unbdd[i,0],"+",dir[i,0],'T')
        print(" ")

else:
    print("The linear program is infeasible, as the basic variables of the final tableau consists of an artificial variable with non-zero value.\nBasic Variables:")
    for i in range(m):
        print("X",int(Xb[i]),"=",(tableau[i+1,n]))