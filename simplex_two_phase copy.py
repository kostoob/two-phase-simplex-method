import numpy as np
from numpy.linalg import inv
from numpy import genfromtxt
from tabulate import tabulate

def header_names(tableau):
    m,n = np.shape(tableau)
    n = n-1
    headers = {}
    for i in range(n):
        headers[ "x" + str( i+1 ) ] = i
        i+=1
    headers["BFS"] = n
    return headers

def standard_input():
    A = genfromtxt('/Users/kaustubhrajimwale/Desktop/Fall 21/IE535/Computing Project/Model12/A.csv', delimiter=',')
    A = np.delete(A, 0, 1)

    b = genfromtxt('/Users/kaustubhrajimwale/Desktop/Fall 21/IE535/Computing Project/Model12/b.csv', delimiter=',')
    b = np.delete(b, 0, 1)

    C = genfromtxt('/Users/kaustubhrajimwale/Desktop/Fall 21/IE535/Computing Project/Model12/C.csv', delimiter=',')
    C = np.delete(C, 0, 1)

    z = input("\nMaximize(1) or Minimize(0): ")
    z = int(z)

    if z != 0:
        C = C * (-1)
    
    return(A,b,C,z)

def add_artifical_variables(A):
    m,n = np.shape(A)
    C_phase1 = np.zeros(((n+m),1))
    for i in range(m):
        x_a = np.zeros((m,1))
        x_a[i,0] = 1
        A = np.concatenate((A, x_a), 1)
        C_phase1[n+i,0] = 1
    return A, C_phase1

def initial_basic_variables(A):
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

def BFS(B,b):
    Binv = inv(B)
    Binv_b = Binv @ b
    return Binv_b, Binv

def find_Cb(C, Xb,A):
    m,n = np.shape(A)
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

def reduced_cost(Cb,A,C,tableau):
    m,n = np.shape(A)
    tableau[0,0:n] = np.subtract((np.transpose(Cb) @ A), np.transpose(C))
    return tableau

def obj(tableau,Cb,A):
    m,n = np.shape(A)
    tableau[0,n] = (np.transpose(Cb) @ tableau[1:m+1, n])
    return tableau

def initiate():
    A,b, C, z = standard_input()
    m,n = np.shape(A)
    A, C_phase1 = add_artifical_variables(A)
    Xb = initial_basic_variables(A)
    B = basic_matrix(A,Xb)
    Binv_b, Binv = BFS(B, b)
    Cb = find_Cb(C_phase1,Xb,A)
    tableau = initiate_tableau(A, Binv_b)
    tableau = reduced_cost(Cb,A,C_phase1,tableau)
    tableau = obj(tableau,Cb,A)
    return (tableau,Xb,C,z)

def simplex_iterations(tableau,Xb):


    ### PRINTER
    print("\n===============================")
    print("STARTING TABLEAU")
    print("===============================\n")
    print(tabulate(np.round(tableau,10), tablefmt="github", headers=header_names(tableau), numalign="right"))
    #PRINTER


    m,n = np.shape(tableau)
    m = m-1
    n = n-1
    iter = 0
    flag = 1
    optimal = 0
    unbdd = 0
    pivot_col = int(0)


    while(flag != 0):
        
        
        iter += 1
        print("\n===============================")
        print("ITERATION ",iter)
        print("===============================\n")
        print(tabulate(np.round(tableau,10), tablefmt="github", headers=header_names(tableau),numalign="right"),"\n\nBasic Variables:",Xb)

        if (np.all(np.round(tableau[0,0:n],10)<=0)):
            flag = 0
            flag = int(flag)
            optimal = 1
            break

        else:
            maxRC = 0
            i=0
            for i in range(n):
                if tableau[0,i]>=maxRC:
                    maxRC = tableau[0,i]
                    pivot_col = i
                i+=1
            
            y = (tableau[1:m+1,pivot_col]<=0)
            if (sum(y)==m):
                flag = 0
                flag = int(flag)
                unbdd = 1
                break
            
            #ratio_test
            ratio = 10000000001
            i=0
            for i in range(m):
                if (tableau[i+1,pivot_col]>0) and ((tableau[i+1,n]/tableau[i+1,pivot_col])<=ratio):
                    if (tableau[i+1,n]/tableau[i+1,pivot_col]) < ratio:
                        ratio = (tableau[i+1,n]/tableau[i+1,pivot_col])
                        pivot_row = i+1

                    else:
                        lexicographic = tableau[i+1,:] - tableau[pivot_row,:]
                        j=0
                        while(i>=0):
                            if lexicographic[j]>0:
                                ratio = (tableau[i+1,n]/tableau[i+1,pivot_col])
                                pivot_row = i+1
                                break
                            else:
                                if lexicographic[j]<0:
                                    break
                                else:
                                    j+=1
                i+=1
            print("\nPivot Element:",tableau[pivot_row,pivot_col],"\n")
            print("Entering Variable: x",pivot_col+1)
            print("Leaving Variable: x",int(Xb[pivot_row-1]),"\n")
            Xb[pivot_row-1] = pivot_col+1



            #Row Transformations
            i=0
            tableau[pivot_row,:] = tableau[pivot_row,:]/tableau[pivot_row,pivot_col]
            for i in range(m+1):
                if i!=pivot_row:
                    tableau[i,:] = tableau[i,:] - (tableau[pivot_row,:]*(tableau[i,pivot_col]))
                i+=1         

                    
    return(tableau,unbdd,optimal,pivot_col, Xb)

def result_phase1(unbdd, optimal, tableau, Xb, C):
    m,n = np.shape(tableau)
    m = m-1
    n = n-1
    phase2 = 1

    logic = Xb > n-m

    if(sum(logic)!=0):
        count = 0
        i = 0
        artificial_variables_in_basis =  np.zeros((0))
        temp = np.zeros((1))
        temp_row = np.zeros((1))
        rows_to_remove = np.zeros((0))
        val = np.zeros((0))

        for i in range(m):
            if Xb[i] > n-m:
                temp[0] = Xb[i]
                artificial_variables_in_basis = np.append(artificial_variables_in_basis,temp,axis= 0)
                temp[0] = round((tableau[i+1,n]),10)
                val = np.append(val,temp,axis= 0)
                if temp == 0:
                    temp_row[0] = i+1
                    rows_to_remove = np.append(rows_to_remove,temp_row,axis= 0)
            i+=1

        if (sum(val == 0)!=0):
            print("\nRedundancy Exists! Removing the redundant constraint.\n")
            l, = np.shape(rows_to_remove)
            i = 0
            print("...")
            for i in range(l):
                tableau = np.delete(tableau, int(rows_to_remove[i]), 0)
                Xb = np.delete(Xb, int(rows_to_remove[i]-1), 0)
                print("...")
                i+=1
            print("\nRedundant Row Deleted\n")

        else:
            print("\nThe linear program is infeasible, as the basic variables of the final tableau consists of one or more artificial variables with non-zero value.\nBasic Variables:")
            for i in range(m):
                print("X",int(Xb[i]),"=",(tableau[i+1,n]))
            print(" ")
            phase2 = 0
            return 0,0,0

    
    else:
        print("\nThe given LP is feasible. Starting Phase II now.\n")
    
    if phase2 == 1:
        i = 0
        for i in range(m):
            tableau = np.delete(tableau, n-i-1, 1)
            i+=1
    
        m,n = np.shape(tableau)
        m = m-1
        n = n-1    
        A = np.zeros((m,n))
        Cb = find_Cb(C,Xb,A)
        zero_col = np.array([[0]])
        temp = np.append(np.transpose(C), zero_col, 1)
        tableau[0,0:n+1] = np.subtract((np.transpose(Cb) @ tableau[1:m+1,:]), temp)
        print("Updated Tableau:\n")
        print(tabulate(np.round(tableau,10), tablefmt="github", headers=header_names(tableau), numalign="right"),"\n\nBasic Variables:",Xb)
        print(" ")
        print("\n End of Phase I. \n Phase II:\n")
        return tableau, Xb, 1

def result_phase2(unbdd,optimal,tableau,Xb,z):
    m,n = np.shape(tableau)
    m = m-1
    n = n-1
    if optimal == 1:
        print(" ")
        print("ITERATIONS STOPPED\n")
        print("The current BFS is optimal:\n")
        if z != 0:
            print("Objective Function =",(tableau[0,n]*(-1)))
        else:
            print("Objective Function =",(tableau[0,n]*(1)))
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


tableau,Xb,C, z = initiate()
tableau,unbdd,optimal,pivot_col, Xb = simplex_iterations(tableau,Xb)
tableau, Xb, next_phase = result_phase1(unbdd, optimal, tableau, Xb,C)
if next_phase == 1:
    tableau,unbdd,optimal,pivot_col, Xb = simplex_iterations(tableau,Xb)
    result_phase2(unbdd, optimal, tableau, Xb, z)