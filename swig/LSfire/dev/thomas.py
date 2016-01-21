def thomas(diagonals, d):
    """ THOMAS    Solves a tridiagonal linear system
    
        x = thomas([c, a, b], d) solves a tridiagonal linear system using the 
            very efficient Thomas algorith. The vector x is the returned answer.
    
           A*x = d;    | a0  b0   0   0   0   ...   0  |   | x0 |   | d0 |
                       | c0  a1  b1   0   0   ...   0  |   | x1 |   | d1 |
                       |  0  c1  a2  b2   0   ...   0  | x | x2 | = | d2 |
                       |  :   :   :   :   :    :    :  |   | x3 |   | d3 |
                       |  0   0   0   0 cn-3 an-2 bn-2 |   | :  |   | :  |
                       |  0   0   0   0   0  cn-2 an-1 |   |xn-1|   |dn-1|
    
       - The matrix A must be strictly diagonally dominant for a stable
         solution.
       
       - This algorithm solves this system on (5*n-4) multiplications/divisions 
         and (3*n-3) subtractions.
       
       - a is the diagonal of the matrix A. If size(a) = N, A is assumed to be
         a NxN matrix.
         
       - b and c are, respectively, the upper and lower diagonals. They should
         have (N-1) components each, i.e. it should be size(b) = size(c) = N-1.
         
       - If size(b) = size(c) = N, the first component of b and the last
         component of c are discarded, i.e. the vectors a, b and c are assumed
         to be the vectors that must be provided to the matlab/scipy spdiags 
         function to build the matrix A:
             A = spdiags(vstack([c, a, b]), [-1, 0, 1], N, N)
             
       Implemented in python by Undy on May 5, 2013, after the file thomas.m
    
    """
    
    c, a, b = diagonals[0], diagonals[1], diagonals[2]

    N = size(a)
        
    if size(b) == N and size(c) == N:
        b = b[1:N]
        c = c[:N-1]
    
    m, x, y = zeros(N), zeros(N), zeros(N)
    #l = zeros(N-1) # it is not necessary to store the whole l array 
    
    # step 1: LU decomposition & forward substitution (L*y=d, for y)
    #
    # L = | 1                |     U =  | m0  r0                |
    #     | l0 1             |          |     m1 r1             |
    #     |    l1 1          |          |        m2 r2          |
    #     |     : : :        |          |         :  :  :       |
    #     |           ln-2 1 |          |                  mn-1 |
    #
    #  ri = bi -> not necessary 
    
    m[0] = a[0]
    y[0] = d[0]
    
    for i in range(1, N):
       l = c[i-1]  / m[i-1]
       m[i] = a[i] - l * b[i-1]
       y[i] = d[i] - l * y[i-1]
            
    # step 2: Backward substitutions (U*x=y, for x)       
    x[N-1] = y[N-1] / m[N-1]
    for i in range(N-2, -1, -1):
        x[i] = (y[i] - b[i] * x[i+1]) / m[i]
       
    return x