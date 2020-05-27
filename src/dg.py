import numpy as np
import scipy.special
import math

class DG:
    
    def __init__(self, case):
        discretization_strategy = case["solver"]["spatial_discretization"]
        
        #TODO


class Line:

    @staticmethod
    def set_nodes_1d(N, vertices):
        """ 
        Sets N+1 nodes in equispaced positions using the vertices indicated
        by vx.
        """
        
        K = vertices.shape[1]
        x = np.zeros((N+1, K))
        for k in range(K):
            for i in range(N+1):
                x[i,k] = i * (vertices[1,k] - vertices[0,k]) / N + vertices[0,k]
                
        return x

    @staticmethod
    def node_indices_1d(N):
        """
        Generates number of node Indices for order N.
        
        >>> Line.node_indices_1d(1)
        array([[1, 0],
               [0, 1]])
            
        >>> Line.node_indices_1d(2)
        array([[2, 0],
               [1, 1],
               [0, 2]])
        """
        Np = N+1
        nId = np.zeros([Np, 2])
        for i in range(Np):
            nId[i] = [N-i, i]     
        return nId.astype(int)
        
    @staticmethod
    def jacobi_gauss_lobatto(alpha, beta, N):
        """
        Compute the order N Gauss Lobatto quadrature points, x, associated
        with the Jacobi polynomial.
        
        >>> Line.jacobi_gauss_lobatto(0.0, 0.0, 1)
        array([-1.,  1.])
        
        >>> Line.jacobi_gauss_lobatto(0,0,3)
        array([-1.       , -0.4472136,  0.4472136,  1.       ])

        >>> Line.jacobi_gauss_lobatto(0,0,4)
        array([-1.        , -0.65465367,  0.        ,  0.65465367,  1.        ])
        
        """
        if N==0:
            return np.array([0.0])
        if N==1:
            return np.array([-1.0, 1.0])
        if N>1:
            x, w = scipy.special.roots_jacobi(N-1, alpha+1, beta+1)
            return np.concatenate(([-1.0], x, [1.0]))
        
        raise ValueError('N must be positive.')

    @staticmethod
    def jacobi_gauss(alpha, beta, N):
        """
        Compute the order N Gauss quadrature points, x, 
        and weights, w, associated with the Jacobi 
        polynomial, of type (alpha,beta) > -1 ( <> -0.5).
        >>> s1 = Line.jacobi_gauss(2,1,0)
        >>> s2 = [-0.2,  2]
        >>> np.allclose(s1,s2)
        True
        >>> s1 = Line.jacobi_gauss(2,1,1)
        >>> s2 = [([-0.54691816,  0.26120387]), ([0.76094757, 0.57238576])]
        >>> np.allclose(s1,s2)
        True
        >>> s1 = Line.jacobi_gauss(2,1,2)
        >>> s2 = [([-0.70882014, -0.13230082,  0.50778763]), ([0.39524241,  0.72312171,  0.21496922])]
        >>> np.allclose(s1,s2)
        True
        """
        x = np.zeros(N)
        w = np.zeros(N)
        if N==0:
            x = -(alpha-beta)/(alpha+beta+2)
            w = 2
            return [x,w]

        # Form symmetric matrix from recurrence.
        J = np.zeros([N+1,N+1])
        h1 = np.zeros(N+1)
        aux = np.zeros(N)
        
        for i in range(N):
            aux[i] = 1+i

        for i in range(N+1):
            h1[i] = 2*i+alpha+beta

        J = np.diag(-0.5*(alpha**2-beta**2)/(h1+2)/h1) \
            + np.diag(2/(h1[0:N]+2) \
            * np.sqrt(aux*(aux+alpha+beta)*(aux+alpha) \
            * (aux+beta)/(h1[0:N]+1)/(h1[0:N]+3)),1)

        eps = np.finfo(np.float).eps

        if (alpha+beta < 10*eps):
            J[0,0] = 0.0

        J = J+np.transpose(J)

        [x,V] = np.linalg.eig(J) 
#        x = np.diag(D)
#        vt = np.transpose(V[0,:])
        w = V[0,:]**2*2**(alpha+beta+1)/(alpha+beta+1)*scipy.special.gamma(alpha+1) \
            * scipy.special.gamma(beta+1)/scipy.special.gamma(alpha+beta+1)

        return [x,w]




    @staticmethod
    def jacobi_polynomial(r, alpha, beta, N):
        """
        Evaluate Jacobi Polynomial
        
        >>> r = Line.jacobi_gauss_lobatto(0,0,1)
        >>> Line.jacobi_polynomial(r, 0, 0, 1)
        array([-1.22474487,  1.22474487])

        >>> r = Line.jacobi_gauss_lobatto(0,0,2)
        >>> Line.jacobi_polynomial(r, 0, 0, 2)
        array([ 1.58113883, -0.79056942,  1.58113883])

        >>> r = Line.jacobi_gauss_lobatto(0,0,3)
        >>> Line.jacobi_polynomial(r, 0, 0, 3)
        array([-1.87082869,  0.83666003, -0.83666003,  1.87082869])
        
        >>> r = Line.jacobi_gauss_lobatto(0,0,4)
        >>> Line.jacobi_polynomial(r, 0, 0, 4)
        array([ 2.12132034, -0.90913729,  0.79549513, -0.90913729,  2.12132034])
        
        """
        PL = np.zeros([N+1,len(r)]) 
        # Initial values P_0(x) and P_1(x)
        gamma0 = 2**(alpha+beta+1) \
                / (alpha+beta+1) \
                * scipy.special.gamma(alpha+1) \
                * scipy.special.gamma(beta+1) \
                / scipy.special.gamma(alpha+beta+1)
        PL[0] = 1.0 / math.sqrt(gamma0)
        if N == 0:
        #    return PL.transpose()
            return PL[0]
        gamma1 = (alpha+1.) * (beta+1.) / (alpha+beta+3.) * gamma0
        PL[1] = ((alpha+beta+2.)*r/2. + (alpha-beta)/2.) / math.sqrt(gamma1)
        
        if N == 1:
#            return PL.transpose()
            return PL[1]
        # Repeat value in recurrence.
        aold = 2. / (2.+alpha+beta) \
            * math.sqrt( (alpha+1.)*(beta+1.) / (alpha+beta+3.))

        # Forward recurrence using the symmetry of the recurrence.
        for i in range(N-1):
            h1 = 2.*(i+1.) + alpha + beta
            anew = 2. / (h1+2.) \
                * math.sqrt((i+2.)*(i+2.+ alpha+beta)*(i+2.+alpha)*(i+2.+beta) \
                            / (h1+1.)/(h1+3.))
            bnew = - (alpha**2 - beta**2) / h1 / (h1+2.)
            PL[i+2] = 1. / anew * (-aold * PL[i] + (r-bnew) * PL[i+1])
            aold = anew

        return PL[N]

    @staticmethod
    def vandermonde_1d(N, r):
        """
        Initialize Vandermonde matrix
        >>> r = Line.jacobi_gauss_lobatto(0,0,2)
        >>> Line.vandermonde_1d(2,r)
        array([[ 0.70710678, -1.22474487,  1.58113883],
               [ 0.70710678,  0.        , -0.79056942],
               [ 0.70710678,  1.22474487,  1.58113883]])
        """
        res = np.zeros([len(r), N+1])
        
        for j in range(N+1):
            res[:,j] = Line.jacobi_polynomial(r, 0, 0, j)
            
        return res

    @staticmethod
    def jacobi_polynomial_grad(r, alpha, beta, N):
        """
        Evaluate the derivative of the Jacobi pol. of type (alpha,beta) > -1
        at points r for order N
        >>> r = Line.jacobi_gauss_lobatto(0,0,1)
        >>> Line.jacobi_polynomial_grad(r,0,0,1)
        array([1.22474487, 1.22474487])
        >>> r = Line.jacobi_gauss_lobatto(0,0,3)
        >>> Line.jacobi_polynomial_grad(r,0,0,3)
        array([11.22497216,  0.        ,  0.        , 11.22497216])
        """
        
        dP = np.zeros([len(r)])

        if N == 0:
            return dP
           
        jPol = Line.jacobi_polynomial(r,alpha+1,beta+1,N-1)
        for i in range(len(r)):
            dP[i] = math.sqrt(N*(N+alpha+beta+1))*jPol[i]
        return dP    

    @staticmethod
    def vandermonde_1d_grad(N,r):
        """
        Initialize the gradient of the modal basis (i) at (r)
        at order (N)
        >>> r = Line.jacobi_gauss_lobatto(0,0,2)
        >>> Line.vandermonde_1d_grad(2,r)
        array([[ 0.        ,  1.22474487, -4.74341649],
               [ 0.        ,  1.22474487,  0.        ],
               [ 0.        ,  1.22474487,  4.74341649]])
        """

        DVr =  np.zeros([len(r),N+1])
        for i in range(N+1):
            DVr[:,i] = Line.jacobi_polynomial_grad(r,0,0,i)
        return DVr

    @staticmethod
    def differentiation_matrix_1d(N,r,V):
        """
        Initialize the (r) differentiation matrices
        of the interval evaluated at (r) at order N
        V is the 1d Vandermonde matrix
        >>> r = Line.jacobi_gauss_lobatto(0,0,2)
        >>> V = Line.vandermonde_1d(2,r)
        >>> Line.differentiation_matrix_1d(2,r,V)
        array([[-1.5,  2. , -0.5],
               [-0.5,  0. ,  0.5],
               [ 0.5, -2. ,  1.5]])

        >>> r = Line.jacobi_gauss_lobatto(0,0,3)
        >>> V = Line.vandermonde_1d(3,r)
        >>> A1 = Line.differentiation_matrix_1d(3,r,V)
        >>> A2 = ([[-3.00000000e+00,  4.04508497e+00, -1.54508497e+00,  5.00000000e-01], \
                   [-8.09016994e-01, -4.05396129e-16,  1.11803399e+00, -3.09016994e-01], \
                   [ 3.09016994e-01, -1.11803399e+00,  6.28036983e-16,  8.09016994e-01], \
                   [-5.00000000e-01,  1.54508497e+00, -4.04508497e+00,  3.00000000e+00]])
        >>> np.allclose(A1,A2)
        True
        """
        Vr = Line.vandermonde_1d_grad(N,r)
        Vinv = np.linalg.inv(V)
        Dr = np.matmul(Vr,Vinv)
        return Dr

    @staticmethod
    def surface_integral_dg(N,V):
        """
        Compute surface integral term in DG formulation
        >>> r = Line.jacobi_gauss_lobatto(0,0,2)
        >>> V = Line.vandermonde_1d(2,r)
        >>> Line.surface_integral_dg(2,V)
        array([[ 4.5 ,  1.5 ],
               [-0.75, -0.75],
               [ 1.5 ,  4.5 ]])

        >>> r = Line.jacobi_gauss_lobatto(0,0,3)
        >>> V = Line.vandermonde_1d(3,r)
        >>> Line.surface_integral_dg(3,V)
        array([[ 8.        , -2.        ],
               [-0.89442719,  0.89442719],
               [ 0.89442719, -0.89442719],
               [-2.        ,  8.        ]])
        """
        # Nfaces, Nfp and Np are defined as global variables
        Nfaces = 1
        Nfp = 2
        Np = N+1

        Emat = np.zeros([Np,Nfaces*Nfp])
        Emat[0,0] = 1.0
        Emat[Np-1,1] = 1.0

        Vtrans = np.transpose(V)
        Vi = np.matmul(Vtrans,Emat)
        Lift = np.matmul(V,Vi)
        return Lift

    @staticmethod
    def normals_1d(K):
        """
        Compute outward pointing normals at element faces
        >>> Line.normals_1d(4)
        array([[-1., -1., -1., -1.],
               [ 1.,  1.,  1.,  1.]])
        """
        # K is the number of elements, derived from the grid info
        # Nfaces and Nfp are defined as global variables
        Nfaces = 1
        Nfp = 2
        nx = np.zeros([Nfp*Nfaces,K])
        nx[0,:] = -1.0
        nx[1,:] = 1.0
        return nx 

    @staticmethod
    def filter_1d(N,Nc,s,V):
        """
        Initialize 1D filter matrix of size N.
        Order of exponential filter is (even) s with cutoff at Nc;
        >>> r = Line.jacobi_gauss_lobatto(0,0,3)
        >>> V = Line.vandermonde_1d(3,r)
        >>> Line.filter_1d(3,1,1,V)
        array([[ 0.3333842 ,  0.9756328 , -0.14240119, -0.1666158 ],
               [ 0.19512656,  0.66667684,  0.16667684, -0.02848024],
               [-0.02848024,  0.16667684,  0.66667684,  0.19512656],
               [-0.1666158 , -0.14240119,  0.9756328 ,  0.3333842 ]])
        """

        s_even = 2*s
        Fdiagonal = np.ones(N+1)
        alpha = -np.log(np.finfo(np.float).eps)

        for i in range(Nc,N+1):
            Fdiagonal[i] = np.exp(-alpha*((i-Nc)/(N-Nc))**s_even)

        # F = V*diag(Fdiagonal)*Vinv
        Fi = np.matmul(V,np.diag(Fdiagonal))
        Vinv = np.linalg.inv(V)
        F = np.matmul(Fi,Vinv)
        return F

if __name__ == '__main__':
    import doctest
    doctest.testmod()