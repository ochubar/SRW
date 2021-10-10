###############################################################################################
# Algorithms solving for the eigensystem
# including: 
#     Eigensolver for symmetric/Hermitian matrices
#     Principle component analysis (PCA)
#
# Eigensolvers: 
#     SciPy (LAPACK): standard or generalized eigenvalue problem 
#	  for a complex Hermitian or real symmetric matrix, 
#	  work size O(N^2) (N: dimension of the matrix)
#     SciPy sparse (ARPACK): find eigenvalues and eigenvectors of 
#         the real symmetric matrix or complex Hermitian matrix, 
#	  work size O(N)
#     PRIMME (multi-method): standard or generalized eigenvalue problem 
#	  for a sparse complex Hermitian or real symmetric matrix, 
#	  with preconditioning, memory-efficient
#
# Reference:
#     https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html
#     https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html
#     https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html
#     http://www.cs.wm.edu/~andreas/software/
#
# Ruizi Li, May 2021
##############################################################################################

#import sys #OC04082021 (commented-out)
from time import time
#from sklearn.decomposition import PCA as PCA #OC04082021 (moved to try-except)
try:
    from scipy.linalg import eigh as largest_eigh
    from scipy.sparse.linalg.eigen.arpack import eigsh as largest_eigsh
except:
    print("UtiMathEigen WARNING: SciPy unavailable; to install: pip install scipy")
try:
    from primme import eigsh as prmm_eigsh
    PRIMME = True
except:
    PRIMME = False
try:
    from sklearn.decomposition import PCA as PCA
except:
    print("UtiMathEigen WARNING: sklearn unavailable; to install: pip install -U scikit-learn")
import math
import numpy as np


#class Linear_alg:
class UtiMathEigen: #OC12052021
#class utilMathEigen:
    def __init__(self, dim=None, fptype = 1, alg='PM', lower=False, rescale=True):
        """
        Initiate the object
        :param dim: rank of the matrix
        :param fptype: precision, 1(2) for single(double)
        :param alg: eigensolver algorithm, 'SP' -- SciPy, 'SPS' -- SciPy Sparse, 'PM' -- Primme 
        :param lower: use only the lower triangle part of the matrix for (SP) eigensolver; valid only when alg='SP'
        """
        if PRIMME:
            self._MODELIST = ("SP", "SPS", "PM") # SciPy; SciPy Sparse; Primme
        else:
            self._MODELIST = ("SP", "SPS")
            print("UtiMathEigen WARNING: PRIMME unavailable; to install: pip install primme") #OC12052021
            #print("utilMathEigen WARNING: PRIMME unavailable; to install: pip install primme")
            print("For more information, please visit https://github.com/primme/primme")
        self.dim = dim
        self.fptype = fptype
        self._eigenVal = None
        self._eigenVecL = None
        self._eigenVecR = None
        if alg not in self._MODELIST:
            print("UtiMathEigen WARNING: undefined / unsupported algorithm {:}, expected {:}".format(alg, self._MODELIST)) #OC12052021
            #print("utilMathEigen WARNING: undefined algorithm {:}, expected {:}".format(alg, self._MODELIST))

            #self.alg = None #OC19062021 (commented-out)
            print('Setting algorithm to \'' + self._MODELIST[0] + '\'') #OC19062021
            self.alg = self._MODELIST[0] #OC19062021
            #self.alg = self._MODELIST[len(self._MODELIST) - 1] #OC19062021
        else:
            self.alg = alg
        self.lower = lower
        self.rescale = rescale
        return


    def eigen_left(self, data, n_modes=None, alg=None, queue=None, pid=None):
        """
        Solve for the eigensystem of input matrix w. largest N eigenvalues and column eigenvectors 
        :param data: input real or complex matrix with dimension (self.dim, self.dim)
        :param n_modes: number of eigenpairs (w. largest eigenvalues) of interest
        :param alg: eigensolver algorithm
        :param queue: list to append the generated eigenpairs to, default not applicable
        :param pid: predefined ID of the generated eigenpairs, applies only when param queue is not None
        :return: (queue is None) array of eigenvalues, array of eigenvectors (normalized) with dimension (rank, N), order of eigenpairs (True in descending order of eigenvalues, False o.w.)
        """
        self._eigenVal = None
        self._eigenVecL = None
        decr = False #OC29062021
        #decr = None
        
        if alg is None:
            alg = self.alg
        if alg not in self._MODELIST:
            print("UtiMathEigen WARNING: undefined / unsupported algorithm {:}, expected {:}".format(alg, self._MODELIST)) #OC29062021
            #print("Unknown eigensolver {:}\n Supported algorithms: {:}\n".format(alg, self._MODELIST))

            if queue is not None:
                queue.put( [ None, None, decr, pid ] )
                return None, None #OC29062021 (to make same return format in all cases)
                #return
            else:
                print('Setting algorithm to \'' + self.alg + '\'') #OC29062021
                alg = self.alg #OC29062021
                #return None, None, decr

        if n_modes is None:
            n_modes = self.dim
        start = time()
        if alg == 'PM':
            self._eigenVal, self._eigenVecL = prmm_eigsh(np.array(data), k=n_modes)
            decr = True #OC29062021 (based on guess)
            elapsed = round(time() - start, 3) #OC28062021
            #elapsed = (time() - start)
            print("Primme eigsh elapsed time: {:} s".format(elapsed)) #OC28062021
            #print("Primme eigsh elapsed time: {:}".format(elapsed))
        if alg == 'SP':
            self._eigenVal, self._eigenVecL = largest_eigh(np.array(data), lower=self.lower, eigvals=(self.dim-n_modes, self.dim-1))
            decr = False
            elapsed = round(time() - start, 3) #OC28062021
            #elapsed = (time() - start)
            print("SciPy eigh elapsed time: {:} s".format(elapsed)) #OC28062021
            #print("SciPy eigh elapsed time: {:}".format(elapsed))
        if alg == 'SPS':
            self._eigenVal, self._eigenVecL = largest_eigsh(data, n_modes, which='LM')
            decr = True #OC29062021 (based on tests made on Windows)
            elapsed = round(time() - start, 3) #OC28062021
            #elapsed = (time() - start)
            print("SciPy sparse eigh elapsed time: {:} s".format(elapsed)) #OC28062021
            #print("SciPy sparse eigh elapsed time: {:}".format(elapsed))
        if self.rescale:

            #print(self._eigenVal) #DEBUG #OC20062021
            self._eigenVecL *= np.sqrt(np.abs(self._eigenVal))
        
        #OC29062021: test actual order of eigenvalues / eigenvectors
        if(self._eigenVal is not None):
            if(len(self._eigenVal) >= n_modes):
                decr = True if(abs(self._eigenVal[0]) > abs(self._eigenVal[n_modes-1])) else False

        if decr:
            self._eigenVecL = np.array(self._eigenVecL.T)
        else:
            self._eigenVecL = np.array(np.flip(self._eigenVecL.T, axis=0))
            self._eigenVal = np.flip(self._eigenVal)
        if queue is not None:
            queue.put([self._eigenVal, self._eigenVecL, decr, pid])
            return
        else:
            return self._eigenVecL, self._eigenVal


    def pca(self, data, n_modes=None):
        """
        PCA algorithm
        :param data: input matrix
        :param n_modes: number of components of interest
        :return: n_modes singular values
        """
        if n_modes is None:
            n_modes = 'mle'
        svd_solver = 'auto'
        pca = PCA(n_components=n_modes, svd_solver=svd_solver)
        pca.fit(data)
        print("PCA: variance ratio: ")
        print(pca.explained_variance_ratio_) 
        return pca.singular_values_


    def angle2eigenvec_left(self, v, indxs = 0, indxe = None):
        """
        Geometric angle between a column vector and a subset of the left eigenspace
        :param v: array of column vector
        :param indxs: left-end index of the eigenvector in the eigenspace
        :param indxe: right-end index of the eigenvector in the eigenspace
        :return: angle between the vector and the eigenspace
        """
        if len(v) != len(self._eigenVecL[indxs]):
            raise Exception("Unmatched vector length of {:}, expecting {:}".format(len(v), len(self._eigenVecL[indxs])))
        if indxe is None:
            indxe = len(self._eigenVecL)
        v = v/np.linalg.norm(v)
        ev = self._eigenVecL.T[indxs:indxe].T
        dot = np.matmul(ev, np.matmul(v.conj(), ev))
        return math.asin(np.linalg.norm(v-dot))


    def reconst_mat_left(self, indxs = 0, indxe = None):
        """
        Reconstruct the matrix from a subset of the left eigenspace
        :param indxs: left-end index of the eigenvector in the eigenspace
        :param indxe: right-end index of the eigenvector in the eigenspace
        :return: reconstructed matrix
        """
        if indxe is None:
            indxe = self._eigenVecL.shape[1]
        print("Shape of eigenVec is {:}".format(self._eigenVecL.shape))
        ev = self._eigenVecL.T[indxs:indxe].T
        ed = np.diag(self._eigenVal[indxs:indxe])
        return np.matmul(np.matmul(ev, ed), ev.conj().T)

