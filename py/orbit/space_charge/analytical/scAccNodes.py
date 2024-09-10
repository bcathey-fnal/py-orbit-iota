"""
Module. Includes classes for the analytical space charge calculator nodes.
"""

import sys
import os
import math
import numpy as np
import orbit_mpi
import time
from orbit_mpi import mpi_datatype
from orbit_mpi import mpi_op

# import the function that finalizes the execution
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

#import the base DirectForce AccNode class
from orbit.space_charge.scAccNodes import SC_Base_AccNode

class SCanalyticalAccNode(SC_Base_AccNode):
	"""
	The subclass of the AccNodeBunchTracker class. It uses SpaceChargeFrozen wrapper for the c++ space charge calculator.
	"""
	def __init__(self, sc_calculator, lattice_functions, name = "no name"):	
		"""
		Constructor. Creates the SC accelerator node element.
		"""
		SC_Base_AccNode.__init__(self, sc_calculator, name)
		self.setType("FrozenSC")
	
	def set_lattice_functions(self, lattice_functions):
		self.lattice_functions = lattice_functions
		
	def track(self, paramsDict):
		"""
		It is tracking the bunch through the Space Charge calculator.
		"""
		bunch = paramsDict["bunch"]
		beta_x = self.lattice_functions['betax']
		beta_y = self.lattice_functions['betay']
		eta_x  = self.lattice_functions['etax']
		eta_y  = self.lattice_functions['etay']
		co_x   = self.lattice_functions['orbitx']
		co_y   = self.lattice_functions['orbity']
		"""
		print 'debug: beta_x =', beta_x
		print 'debug: beta_y =', beta_y
		print 'debug: eta_x  =', eta_x
		print 'debug: eta_y  =', eta_y
		print 'debug: co_x   =', co_x
		print 'debug: co_y   =', co_y
		"""
		self.sc_calculator.setLatticeParameters(beta_x, beta_y, eta_x, eta_y, co_x, co_y)
		self.sc_calculator.trackBunch(bunch, self.sc_length)		
