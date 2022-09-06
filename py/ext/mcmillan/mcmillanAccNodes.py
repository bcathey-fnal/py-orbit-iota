"""
Module. Includes abstract classes for all types of McMillan accelerator nodes.
Based on SC_Base_AccNode
"""
# import the function that finalizes the execution
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

class McMillan_AccNode(AccNodeBunchTracker):
    """
    The subclass of the AccNodeBunchTracker class. It is a base class for monitoring nodes. It uses the c++ McMillan class.
    """
    def __init__(self, mcmillan_calc, length, ks, name = "no name"):			
        """
        Constructor. Creates the mcmillan accelerator node element.
        """
        AccNodeBunchTracker.__init__(self,name)
        self.setType("McMillan")
        self.setLength(length) # Length of thick lens
        self.switcher = True # Switches the calculation on or off
        self.addParam('ks', ks) # Solenoid field
        self.mcmillan_calc = mcmillan_calc

    def __str__(self):
        return '{}: {}'.format(self.getName(), self.getType())

    def isRFGap(self):
        """
        Returns False. The RF Gap node returns True.
        """
        return False
		
    def setCalculationFlag(self, switcher):
        """
        Sets the boolean parameter that define if the calculations will be performed. True of False.
        """
        self.switcher = switcher
		
    def getCalculationFlag(self):
        """
        Returns the boolean parameter that define if the calculations will be performed. True of False.
        """
        return 	self.switcher	
		
    def getCalculator(self):
        """
        Returns the scattering calculator.
        """
        return self.mcmillan_calc
		
    def trackDesign(self, paramsDict):
        """
        This method is for Linac Nodes compatibility. It is empty and should not be used for anything else.
        """
        pass	
		
    def track(self, paramsDict):
        """
        Track the bunch through the mcmillan lens
        """
        if(self.switcher != True): return
        bunch = paramsDict["bunch"]
        ks = self.getParam('ks')
        track_length = self.getLength(self.getActivePartIndex())
        #print('Tracking through {}:{} l = {}, ks = {}'.format(self.getName(),
        #    self.getActivePartIndex(), track_length, ks))
        self.mcmillan_calc.trackBunch(bunch, track_length, ks)
