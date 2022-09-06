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
        self.TFloquet = np.eye(4).tolist()

    def setTFloquet(self, TFloquet):
        #pass
        self.TFloquet = TFloquet.tolist()

    def __str__(self):
        return '{}: TuneDetector, l = 0, TFloquet = {}'.format(self.getName(), self.TFloquet)
				
    def track(self, paramsDict):
        """
        Track the bunch through the tune detector
        """
        if(self.switcher != True): return
        bunch = paramsDict["bunch"]
        self.monitoring_calc.trackBunch(bunch, self.TFloquet)
		
