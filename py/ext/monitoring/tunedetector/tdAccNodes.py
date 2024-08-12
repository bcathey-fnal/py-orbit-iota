"""
Module. Tune detector node class.
"""
# import the function that finalizes the execution
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
#from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

#import the base monitoring AccNode class
from ext.monitoring.monitoringAccNodes import Monitoring_Base_AccNode

import numpy as np # For matrix calculations

class TuneDetector_AccNode(Monitoring_Base_AccNode):
    """
    The subclass of the Monitoring_Base_AccNode which is in turn a subclass of AccNodeBunchTracker
    """
    def __init__(self, monitoring_calc, name = "no name"):			
        """
        Constructor. Creates the tunedetector node element.
        """
        Monitoring_Base_AccNode.__init__(self, monitoring_calc, "TuneDetector", name)
        self.setTmode(np.eye(4), np.zeros((6,6)))
        self.isfirst = False

    def setTmode(self, Vinv, T):
        """
        Set the transform matrix for the node. By default it's identity.
        
        Parameters:
        - Vinv: 4x4 matrix to transform into normal mode coordinates.
        - T: 6x6 one-turn matrix of the linear lattice at the position.

        Returns: None
        """
        # Calculate the dispersion vector from the one-turn matrix
        Dvec = np.dot(np.linalg.inv(np.eye(4) - T[:4,:4]), T[0:4,5])
        Tmode = np.zeros((4,5)) # prepare a 4x5 matrix
        Tmode[:,:4] = Vinv # Assign the Floquet matrix
        Tmode[:,4] = Dvec*0.0 # Assign the dispersion vector
        self.Tmode = Tmode.tolist() # Convert the numpy array to a list

    def isFirst(self, isfirst):
        """
        Flag whether this is the first tune detector node in the lattice.
        """
        self.isfirst = isfirst

    def __str__(self):
        return '{}: TuneDetector, l = 0, Tmode = {}, enable={}'.format(
                self.getName(), self.Tmode, self.enable)
				
    def track(self, paramsDict):
        """
        Track the bunch through the tune detector
        """
        if(self.switcher != True): return
        bunch = paramsDict["bunch"]
        self.monitoring_calc.trackBunch(bunch, self.Tmode, self.isfirst)
		
