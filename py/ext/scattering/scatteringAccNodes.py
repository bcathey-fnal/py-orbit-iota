"""
Module. Includes abstract classes for all types of scattering accelerator nodes.
Based on SC_Base_AccNode
"""
# import the function that finalizes the execution
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

class Scattering_Base_AccNode(AccNodeBunchTracker):
    """
    The subclass of the AccNodeBunchTracker class. It is a base class for scattering nodes. It uses the c++ scattering calculators.
    """
    def __init__(self, scattering_calc, name = "no name"):			
        """
        Constructor. Creates the scattering accelerator node element.
        """
        AccNodeBunchTracker.__init__(self,name)
        self.setType("Scattering_Base")
        self.scattering_length = 0.
        self.switcher = True
        self.scattering_calc = scattering_calc
                
    def isRFGap(self):
        """
        Returns False. The RF Gap node returns True.
        """
        return False
		
    def setLengthOfScattering(self, scattering_length):
        """
        Defines the path length that will be used in scattering kick calculations.
        """
        self.scattering_length = scattering_length
		
    def getLengthOfScattering(self):
        """
        Returns the path length that will be used in scattering kick calculations.
        """		
        return self.scattering_length
		
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
        return self.scattering_calc
		
    def trackDesign(self, paramsDict):
        """
        This method is for Linac Nodes compatibility. It is empty and should not be used for Space Charge calculations.
        """
        pass	
		
    def track(self, paramsDict):
        """
        It is tracking the bunch through the scattering calculator. Each subclass should implement this method.
        """
        pass

