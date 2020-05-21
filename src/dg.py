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
    def jacobi_polynomial(r, alpha, beta, N):
        """
        Evaluate Jacobi Polynomial
        
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
            return PL.transpose()
        
        gamma1 = (alpha+1.) * (beta+1.) / (alpha+beta+3.) * gamma0
        PL[1] = ((alpha+beta+2.)*r/2. + (alpha-beta)/2.) / math.sqrt(gamma1)
        
        if N == 1:
            return PL.transpose()

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
        """
        res = np.zeros([len(r), N+1])
        
        for j in range(N+1):
            res[j] = Line.jacobi_polynomial(r, 0, 0, j)
            
        return res.transpose()

if __name__ == '__main__':
    import doctest
    doctest.testmod()