import cvxpy as cp
import numpy as np

def skew_symmetric(V):
    skew = cp.vstack([cp.hstack([0, -V[2, 0], V[1, 0]]),
                      cp.hstack([V[2, 0], 0, -V[0, 0]]),
                      cp.hstack([-V[1, 0], V[0, 0], 0])])
    return skew

def skew_symmetric_neg(V):
    skew = cp.vstack([cp.hstack([0, V[2, 0], -V[1, 0]]),
                      cp.hstack([-V[2, 0], 0, V[0, 0]]),
                      cp.hstack([V[1, 0], -V[0, 0], 0])])
    return skew

def my_kronecker(A, B):
    expr = []
    dimA = A.shape
    dimB = B.shape
    for i in range(dimA[0]):
        for k in range(dimB[0]):
            temp_list = []
            for j in range(dimA[1]):
                for l in range(dimB[1]):
                    temp_list.append(-A[i, j] * B[k, l])
            expr.append(temp_list)
    return expr

def np_to_cp(expr):
    num_rows = len(expr)
    rows = [cp.hstack(expr[i]) for i in range(num_rows)]
    full_expr = cp.vstack(rows)
    return full_expr

def add_lin_combination(X, list):
    #input must be a numpy matrix
    #adds the linear combination constraint
    list_len= len(list)
    [list[i].append(cp.Constant(0)) for i in range(list_len)]
    temp_list = [cp.hstack(list[i]) for i in range(list_len)]
    temp_list.append(cp.hstack([cp.Constant(0),
                           cp.Constant(0),
                           cp.Constant(0),
                           cp.Constant(0),
                           cp.Constant(0),
                           cp.Constant(0),
                           cp.Constant(0),
                           cp.Constant(0),
                           cp.Constant(0),
                           X[0,0] + X[1,1] + X[2,2]]))
    PX = cp.vstack(temp_list)
    return PX

def construct_C(X):
    ''' Column Orthogonality Constraints'''
    P1_kron_np = my_kronecker(X, np.eye(3))
    P1 = add_lin_combination(X, P1_kron_np)
    return P1

def construct_R(X):
    ''' Row Orthogonality Constraints'''
    P2_kron_np = my_kronecker(np.eye(3), X)
    P2 = add_lin_combination(X, P2_kron_np)
    return P2

def construct_H(X_ijk, X_jki, X_kij):
    '''Right-Handedness Constraints'''
    zero_matrix = np.zeros((3,3))
    zero_vector = np.zeros((3,1))
    row1 = cp.hstack([zero_matrix,
                      skew_symmetric_neg(X_ijk),
                      skew_symmetric(X_kij),
                      -X_jki])
    row2 = cp.hstack([skew_symmetric(X_ijk),
                      zero_matrix,
                      skew_symmetric_neg(X_jki),
                      -X_kij])
    row3 = cp.hstack([cp.vstack([skew_symmetric_neg(X_kij), -X_jki.T]),
                      cp.vstack([skew_symmetric(X_jki), -X_kij.T]),
                      cp.vstack([zero_matrix, -X_ijk.T]),
                      cp.vstack([-X_ijk, cp.Constant([[0]])])])
    P3 = cp.vstack([row1, row2, row3])
    return P3

def construct_y(X):
    expr = np.zeros((10,10)).tolist()
    expr[9][9] = -X
    P4 = np_to_cp(expr)
    return P4

def construct_P(constraints='RCH'):
    ''' Create the constraints requested'''
    #Create the y constraint
    X_4 = cp.Variable(1)
    P = construct_y(X_4)

    #Create list of requested constraints
    cons = list(constraints)

    #Implement constraints
    if 'R' in cons:
        X_2 = cp.Variable((3, 3), symmetric=True)
        P = P + construct_R(X_2)
    if 'C' in cons:
        X_1 = cp.Variable((3, 3), symmetric=True)
        P = P + construct_C(X_1)
    if 'H' in cons:
        X_ijk = cp.Variable((3, 1))
        X_jki = cp.Variable((3, 1))
        X_kij = cp.Variable((3, 1))
        P = P + construct_H(X_ijk, X_jki, X_kij)
    return P, X_4

def construct_A_and_b(constraints='RCH', debug=False):
    A = [] # trace(A@X) == b
    b = [] # trace(A@X) == b
    p = 0  # size of arrays

    # Create list of requested constraints
    cons = list(constraints)

    #Implement constraints
    # Iterate over each neu and give it a value of 1 while all others are 0
    if 'C' in cons:
        for i in range(3):
            for j in range(i,3):
                array = np.zeros([3,3])
                array[i,j] = 1
                array[j,i] = 1

                X_1 = cp.Constant(array)
                P_1 = construct_C(X_1)
                P_1_value = P_1.value
                if P_1_value[-1,-1] != 0:
                    P_1_value[-1,-1]  = 1

                A.append(P_1_value)
                b.append(0)
                p += 1

    if 'R' in cons:
        for i in range(3):
            for j in range(i,3):
                array = np.zeros([3,3])
                array[i,j] = 1
                array[j,i] = 1

                X_2= cp.Constant(array)
                P_2 = construct_R(X_2)
                P_2_value = P_2.value
                if P_2_value[-1,-1] != 0:
                    P_2_value[-1,-1]  = 1

                A.append(P_2_value)
                b.append(0)
                p += 1

    if 'H' in cons:
        for i in range(3):
            for j in range(3):
                array = np.zeros([3,3])
                array[i,j] = 1

                X_ijk = cp.Constant(array[:,0:1])
                X_jki = cp.Constant(array[:,1:2])
                X_kij = cp.Constant(array[:,2:3])
                P_3 = construct_H(X_ijk, X_jki, X_kij)

                A.append(P_3.value)
                b.append(0)
                p += 1

    X_4 = cp.Constant(1)
    P_4 = construct_y(X_4)

    A.append(P_4.value)
    b.append(-1)
    p += 1

    if debug:
        for i in range(p):
            print('#### Constraint {}'.format(i+1))
            print(A[i])
            print("\t=\n\t=")
            print(b[i])

    return A, b, p


def main():
    P = construct_P(constraints='RCH')
    A, b, p = construct_A_and_b(constraints='RCH')
    return

if __name__ == "__main__":
    main()