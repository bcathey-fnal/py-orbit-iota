"""
Module. Electron cooling node class.
"""
# import the function that finalizes the execution
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
#from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

#import the base SC AccNode class
from ext.scattering.scatteringAccNodes import Scattering_Base_AccNode

class ECooler_AccNode(Scattering_Base_AccNode):
    """
    The subclass of the Scattering_Base_AccNode which is in turn a subclass of AccNodeBunchTracker
    """
    def __init__(self, scattering_calc, name = "no name"):			
        """
        Constructor. Creates the ecooler node element.
        """
        Scattering_Base_AccNode.__init__(self, scattering_calc, name)
        self.setType("ECooler")
        self.xc = self.yc = 0.0 # The default offset is 0!

    def setoffset(self, xc, yc):
        self.xc = xc
        self.yc = yc

    def __str__(self):
        return '{}: ECooler, l = {}, xc = {}, yc = {}'.format(self.getName(),
            self.scattering_length, self.xc, self.yc)
				
    def track(self, paramsDict):
        """
        Track the bunch through the ecooler
        """
        if(self.switcher != True): return
        bunch = paramsDict["bunch"]
        self.scattering_calc.trackBunch(bunch, self.scattering_length, self.xc, self.yc)
		
