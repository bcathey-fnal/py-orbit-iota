##############################################################################
"""
Module. Includes the Linear_Coupled_Lattice class which can be used for linear
analysis for strongly coupled lattices. Based on
https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.2.074001
We obtain the element matrices using  TEAPOT_MATRIX_Lattice and then apply
the analysis discussed in
https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.2.074001 to
generate twiss functions of the normal modes. This analysis closely matches
with the code BMAD. https://www.classe.cornell.edu/bmad/
"""
##############################################################################
# Imports
from bunch import Bunch # For testing
from orbit.matrix_lattice import BaseMATRIX # For extracting matrices
from orbit.teapot import TEAPOT_MATRIX_Lattice # Matrix lattice
import numpy as np # For linear algebra

##############################################################################
# Constants
S = np.zeros((4,4), dtype='complex128') # Symplectic matrix
S[0,1] = S[2,3] = 1.0; S[1,0] = S[3,2] = -1.0
one_over_sqrt2 = 2**(-0.5) # Useful constant

##############################################################################
class Linear_Coupled_Lattice(TEAPOT_MATRIX_Lattice):
    """
    A subclass for the TEAPOT_MATRIX_Lattice implementing strong coupling.
    """
    def __init__(self, teapot_lattice, bunch, name = None):
        """
        Constructor. Creates a linear coupled lattice object and performs
        normal mode analysis.
        """
        # Initialize the teapot matrix lattice
        TEAPOT_MATRIX_Lattice.__init__(self, teapot_lattice, bunch, name)

        # Analyze optics
        self.analyze()

    def analyze(self):
        """
        Perform linear analysis of the lattice.
        """
        # Obtain the 1-turn map at lattice start
        T = to_numpy_matrix(self.getRingMatrix())
        s = 0.0 # The s position at the start of the lattice

        # Iterate through all nodes and find the one-turn matrices
        # at the position of each nodes.
        for matrixNode in self.getNodes():
            # Check that the node is of the type BaseMATRIX
            if(isinstance(matrixNode,BaseMATRIX) == True):
                # Extract the node length and matrix
                node_l = matrixNode.getLength()
                node_matrix = to_numpy_matrix(matrixNode.getMatrix())
                N = action_angle_matrix(T) # Normalization matrix
                # Obtain the similarity transformation matrix
                V, U, Vinv = edwards_teng_transformation(T)
                Tfl = contruct_uncoupled_floquet(U)
                # Add the s position and various parameters to the nodes
                matrixNode.addParam('s', s)
                matrixNode.addParam('T', T)
                matrixNode.addParam('N', N)
                matrixNode.addParam('V', V)
                matrixNode.addParam('U', U)
                matrixNode.addParam('Vinv', Vinv)
                matrixNode.addParam('Nfl', np.dot(Tfl, Vinv))

                # Update the s position and the 1-turn matrix
                s += node_l
                T = np.dot(node_matrix, np.dot(T, np.linalg.inv(node_matrix)))

    def getRingTwissData(self):
        """
        Obtain the normal mode twiss data of the whole lattice.
        """
        cnt_nodes = 0 # Count the number of nodes in the lattice
        for matrixNode in self.getNodes():
            cnt_nodes += 1

        data = np.zeros((cnt_nodes,7)) # Create the data array

        # Iterate over all elements
        i = 0 # Counter
        for matrixNode in self.getNodes():
            # Get the un-coupled one-turn matrices
            U = matrixNode.getParam('U')
            UA = U[0:2, 0:2].copy()
            UB = U[2:4, 2:4].copy()
            # Extract the normal mode Twiss parameters
            beta, alfa = extract_twiss_parameters(UA)
            betb, alfb = extract_twiss_parameters(UB)

            # Extract the inverse of the decoupling matrix of this node
            Vinv = matrixNode.getParam('Vinv')

            #print('============================================================')
            #print('{}: s = {}'.format(matrixNode.getName(),matrixNode.getParam('s')))

            # Extract the phase advance from the last node to this node
            if i > 0: # Only do this starting from node 1
                # Extract the transfer matrix for this node
                M = to_numpy_matrix(matrixNode.getMatrix())
                # Calculate the un-coupled transfer matrices
                Munc = np.dot(Vinv, np.dot(M_prev, V_prev))
                # Extract the eigen values and eigenvectors of the un-coupled
                # matrices
                muA = extract_phase_advance(data[i-1,2], beta, data[i-1,3],
                                            alfa, Munc[0:2, 0:2])
                muB = extract_phase_advance(data[i-1,5], betb, data[i-1,6],
                                            alfb, Munc[2:4, 2:4])

            # Finally save the data 
            data[i, 0] = matrixNode.getParam('s') # Put s position
            data[i, 2] = beta; data[i, 3] = alfa
            data[i, 5] = betb; data[i, 6] = alfb
            if i > 0:
                data[i, 1] = data[i-1, 1] + muA
                data[i, 4] = data[i-1, 4] + muB

            i += 1 # Update the counter
            # Save the decoupling matrix for this node
            V_prev = matrixNode.getParam('V')
            # Extract the transfer matrix for this node
            M_prev = to_numpy_matrix(matrixNode.getMatrix())

        return data

    def getFloquet4D(self):
        """
        Get the Floquet matrix for 4D dynamics at lattice entrance
        """
        R = np.zeros((4,4)) # Rotation matrix in normalized phase space
        # One turn matrix at the entrance of the lattice
        M = to_numpy_matrix(self.getRingMatrix())
        wM, vM = np.linalg.eig(M) # Get transverse eigenvalues
        # Fill in the rotation matrix
        R[0:2,0:2] = [[np.real(wM[0]), np.imag(wM[0])],
                      [np.imag(wM[1]), np.real(wM[1])]]
        R[2:4,2:4] = [[np.real(wM[2]), np.imag(wM[2])],
                      [np.imag(wM[3]), np.real(wM[3])]]
        # Get the eigen vectors for the rotation matrix.
        # Eigenvalues are exactly same by construction.
        wR, vR = np.linalg.eig(R)
        # Finally generate Tf such that M = inv(Tf) @ R @ Tf
        Tf = np.dot(vR, np.linalg.inv(vM))
        #Tfinv = np.dot(vM, np.linalg.inv(vR))
        return np.real(Tf) # Clip the imaginary parts, which are ~0

##############################################################################
def to_numpy_matrix(orbit_matrix, size=4):
    """
    Convert from the PyORBIT matrix format to numpy format
    """
    a = np.zeros((size,size)) # Declare an empty numpy matrix
    # Populate the matrix
    for i in range(size):
        for j in range(size):
            a[i,j] = orbit_matrix.get(i,j)
    return a

def symplectic_conjugate(C):
    """
    Calculate symplectic conjugate of a 2x2 matrix.
    """
    if C.shape[0] != 2 or C.shape[1] != 2:
        raise ValueError('Argument must be a 2x2 matrix.')
    Cplus = np.zeros((2,2), C.dtype)
    Cplus[0,0] = C[1,1]
    Cplus[0,1] = -C[0,1]
    Cplus[1,0] = -C[1,0]
    Cplus[1,1] = C[0,0]
    return Cplus

def eigen_analysis(T):
    """
    Calculate eigenvalues and eigenvectors of the node matrix. Order
    according to the order section 21.2 in
    https://www.classe.cornell.edu/bmad/bmad-manual-2022-09-07.pdf
    """
    w_raw, v_raw = np.linalg.eig(T) # First obtain eigenvalues
    v = np.zeros((4,4), v_raw.dtype) # Allocate space for the eigenvectors
    w = np.zeros((4,), w_raw.dtype) # Allocate space for eigenvalues
    for p in range(2): # Iterate over planes
        # Find which eigenvectors have the maximum contribution to the
        # position in the plane with index p
        kmax = np.argmax(np.abs(v_raw[2*p,:]))

        vsv = np.dot(np.conj(v_raw[:,kmax]), np.dot(S, v_raw[:,kmax]))
        alpha = 1.0/np.sqrt(np.abs(vsv.imag)) # Scale factor

        if vsv.imag >= 0.0: # Positive values corresponds to vectors 0, 2
            v[:,2*p] = alpha*v_raw[:,kmax]
            w[2*p] = w_raw[kmax]
            # Find the pair!
            for kmax2 in range(4):
                # The pair eigenvalue must be the conjugate of the other one
                if abs(w_raw[kmax].imag + w_raw[kmax2].imag) < 1e-10 \
                        and kmax != kmax2:
                    v[:,2*p+1] = alpha*v_raw[:,kmax2]
                    w[2*p+1] = w_raw[kmax2]
                    break
        else: # Negative values corresponds to vectors 1, 3
            v[:,2*p+1] = alpha*v_raw[:,kmax]
            w[2*p+1] = w_raw[kmax]
            # Find the pair!
            for kmax2 in range(4):
                # The pair eigenvalue must be the conjugate of the other one
                if abs(w_raw[kmax].imag + w_raw[kmax2].imag) < 1e-10 \
                        and kmax != kmax2:
                    v[:,2*p] = alpha*v_raw[:,kmax2]
                    w[2*p] = w_raw[kmax2]
                    break
    return w, v # return the eigenvalue and vectors

def action_angle_matrix(T):
    """
    Calculate normalization matrix. See Eq. 21.33 in
    https://www.classe.cornell.edu/bmad/bmad-manual-2022-09-07.pdf
    """
    w, v = eigen_analysis(T)

    qa = np.angle(w[0])/(2.0*np.pi)
    qb = np.angle(w[2])/(2.0*np.pi)

    # Create the normalization matrix
    N = np.zeros((4,4), v.dtype)
    # All matrix elements are real
    N[0,:] = one_over_sqrt2*np.transpose(v[:,0]+v[:,1])
    N[1,:] = -1j*one_over_sqrt2*np.transpose(v[:,0]-v[:,1])
    N[2,:] = one_over_sqrt2*np.transpose(v[:,2]+v[:,3])
    N[3,:] = -1j*one_over_sqrt2*np.transpose(v[:,2]-v[:,3])

    return N

def edwards_teng_transformation(T):
    """
    Compute the similarity transformation T = VUV^{-1} using the Edwards
    and Teng formalism.
    """
    # First extract the important matrices
    M = T[0:2, 0:2].copy()
    N = T[2:4, 2:4].copy()
    m = T[0:2, 2:4].copy()
    n = T[2:4, 0:2].copy()
    H = m + symplectic_conjugate(n)

    # Evaluate important quantities for the similarity transformation
    TrMminusN = np.trace(M-N)
    detH = np.linalg.det(H)
    gamma = np.sqrt(0.5+0.5/np.sqrt(1.0 + 4.0*detH/(TrMminusN*TrMminusN)))
    C = -H*np.sign(TrMminusN)/(gamma*np.sqrt(TrMminusN*TrMminusN+4.0*detH))

    # Construct the transformation matrices
    V = gamma*np.eye(4)
    V[0:2, 2:4] = C
    V[2:4, 0:2] = -symplectic_conjugate(C)
    Vinv = gamma*np.eye(4)
    Vinv[0:2, 2:4] = -C
    Vinv[2:4, 0:2] = symplectic_conjugate(C)
    U = np.dot(Vinv,np.dot(T,V))

    return V, U, Vinv

def extract_twiss_parameters(U):
    """
    Extract twiss parameters from the normal mode matrix.
    """
    cos_theta = 0.5*(U[0,0] + U[1,1]) # Cosine of the phase advance
    abs_sin_theta = (1.0 - cos_theta*cos_theta)**0.5 # Absolute value of sine
    beta = abs(U[0,1])/abs_sin_theta # Now get the beta function
    sin_theta = np.sign(U[0,1])*abs_sin_theta # get the sign of the sine value
    alpha = 0.5*(U[0,0] - U[1,1])/sin_theta # Obtain value of alpha
    return beta, alpha

def contruct_uncoupled_floquet(U):
    """
    Construct the 4D Floquet matrix for the uncoupled transverse modes
    to convert to normalized coordinates.
    """
    UA = U[0:2, 0:2].copy()
    UB = U[2:4, 2:4].copy()
    # Extract the normal mode Twiss parameters
    beta, alfa = extract_twiss_parameters(UA)
    betb, alfb = extract_twiss_parameters(UB)
    Tfl = np.zeros((4,4)) # Create a Floquet matrix
    # Mode a
    Tfl[0,0] = beta**(-0.5)
    Tfl[1,0] = alfa*beta**(-0.5)
    Tfl[1,1] = beta**0.5
    # Mode b
    Tfl[2,2] = betb**(-0.5)
    Tfl[3,2] = alfb*betb**(-0.5)
    Tfl[3,3] = betb**0.5
    return Tfl # Return the Floquet matrix


def extract_phase_advance(beta1, beta2, alfa1, alfa2, M):
    """
    Extract phase advance from a transport matrix given the initial
    and final Twiss parameters.
    """
    # M_{12} = \sqrt{\beta_1 \beta_2} \sin \mu
    sin_mu = M[0,1]*(beta1*beta2)**(-0.5)
    # M_{11} = \sqrt{\frac{\beta_2}{\beta_1}} (\cos \mu + \alpha_1 \sin \mu)
    cos_mu_1 = M[0,0]*(beta1/beta2)**0.5 - alfa1*sin_mu
    # M_{11} = \sqrt{\frac{\beta_1}{\beta_2}} (\cos \mu - \alpha_2 \sin \mu)
    # cos_mu_2 = M[1,1]*(beta2/beta1)**0.5 + alfa2*sin_mu
    mu = np.angle(cos_mu_1 + 1j*sin_mu)/(2.0*np.pi) # Extract the phase advance
    return mu
