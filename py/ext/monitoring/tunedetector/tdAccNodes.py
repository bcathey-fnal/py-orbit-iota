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
        self.Tmode = np.eye(4).tolist()
        self.isfirst = False

    def setTmode(self, Tmode):
        """
        Set the transform matrix for the node. By default it's identitity.
        """
        self.Tmode = Tmode.tolist()

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
		
