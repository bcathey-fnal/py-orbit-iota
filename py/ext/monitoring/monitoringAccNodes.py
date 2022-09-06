"""
Module. Includes abstract classes for all types of monitor accelerator nodes.
Based on SC_Base_AccNode
"""
# import the function that finalizes the execution
#from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

class Monitoring_Base_AccNode(AccNodeBunchTracker):
    """
    The subclass of the AccNodeBunchTracker class. It is a base class for monitoring nodes. It uses thee c++ monitor classes.
    """
    def __init__(self, monitoring_calc, typename, name = "no name"):			
        """
        Constructor. Creates the scattering accelerator node element.
        """
        AccNodeBunchTracker.__init__(self,name)
        self.setType(typename)
        self.monitoring_length = 0. # Something for future?? Striplines maybe...
        self.switcher = True
        self.monitoring_calc = monitoring_calc

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
        return self.monitoring_calc
		
    def trackDesign(self, paramsDict):
        """
        This method is for Linac Nodes compatibility. It is empty and should not be used for anything else.
        """
        pass	
		
    def track(self, paramsDict):
        """
        Track the bunch through the monitor
        """
        if(self.switcher != True): return
        bunch = paramsDict["bunch"]
        self.monitoring_calc.trackBunch(bunch)
