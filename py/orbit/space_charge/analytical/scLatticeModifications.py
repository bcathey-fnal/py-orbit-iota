"""
Module. Includes functions that will modify the accelerator lattice by inserting the SC accelerator nodes.
"""

# import acc. nodes
from scAccNodes import SCanalyticalAccNode

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import the general SC lattice modification function 
from orbit.space_charge.scLatticeModifications import setSC_General_AccNodes



def setSCanalyticalAccNodes(lattice, sc_path_length_min, space_charge_calculator):
	"""
	It will put a set of frozenSC_AccNodes into the lattice as child nodes of the first level accelerator nodes.
	The SC nodes will be inserted at the beginning of a particular part of the first level AccNode element.
	The distance between SC nodes should be more than sc_path_length_min	The function will return the array of SC nodes as a convenience for the user.
	"""
	accNodes = lattice.getNodes()
	if(len(accNodes) == 0): return
	#-----------------------------------------------
	# nodes_arr[(accNode, part_index, position, path_length)] 
	#-----------------------------------------------
	nodes_arr = []
	length_total = 0.
	running_path = 0.
	rest_length = 0.
	for i, accNode in enumerate(accNodes):
		nParts = accNode.getnParts()
		for ip in range(nParts):
			part_length = accNode.getLength(ip)
			#if(part_length > 1.0):
			#	print "Warning! Node ",accNode.getName(), " has length ", part_length, "m which is greater than 1 m.  Space charge algorithm may be innacurate!"
			if(running_path > sc_path_length_min):
				# the lattice functions at the entrance of the present node are the ones at the exit of the previous node
				lattice_functions_dict = accNodes[i-1].getParamsDict()
				# if the lattice does not provide closed orbit or vertical dispersion
				err_mess = []
				opticsparameters2check = ['etay', 'orbitx', 'orbity']
				for param in opticsparameters2check:
					try:    
						accNodes[i-1].getParamsDict()[param]
					except: 
						accNodes[i-1].addParam(param, 0.)
						err_mess.append(param)
				if err_mess:
					print 'WARNING: Lattice node', accNodes[i-1].getParamsDict()['node_index'], 'does not provide:', err_mess
				nodes_arr.append((accNode, lattice_functions_dict, ip, length_total, running_path))
				running_path = 0.
			running_path += part_length
			length_total += part_length
	if(len(nodes_arr) > 0):
		rest_length = length_total - nodes_arr[len(nodes_arr) - 1][-2] # this index must be -2 since we added the lattice_functions_dict!!
	else:
		rest_length = length_total
	# now the first SC node in the beginning of the lattice
	# if the lattice does not provide closed orbit or vertical dispersion				
	err_mess = []
	opticsparameters2check = ['etay', 'orbitx', 'orbity']
	for param in opticsparameters2check:
		try:    
			accNodes[-1].getParamsDict()[param]
		except: 
			accNodes[-1].addParam(param, 0.)
			err_mess.append(param)
	if err_mess:
		print 'WARNING: Lattice node', accNodes[-1].getParamsDict()['node_index'], 'does not provide:', err_mess
	nodes_arr.insert(0,(accNodes[0],accNodes[-1].getParamsDict(),0,0.,rest_length))
	#---------------------------------------------------
	# Now we put all SC nodes as a childeren of accNodes
	#---------------------------------------------------
	scNodes_arr = []
	for inode in range(len(nodes_arr)-1):
		(accNode, lattice_functions_dict, part_index, position, path_length) = nodes_arr[inode]
		(accNodeNext, lattice_functions_dictNext, part_indexNext, positionNext, path_lengthNext) = nodes_arr[inode+1]
		scNode = SCanalyticalAccNode(space_charge_calculator,accNode.getName()+":"+str(part_index)+":")
		scNode.setLengthOfSC(path_lengthNext)
		scNode.setName(scNode.getName()+"frozen_SC")
		scNode.set_lattice_functions(lattice_functions_dict)
		scNodes_arr.append(scNode)
		accNode.addChildNode(scNode,AccNode.BODY,part_index,AccNode.BEFORE)
	#set the last SC node
	(accNode, lattice_functions_dict, part_index, position, path_length) = nodes_arr[len(nodes_arr)-1]
	scNode = SCanalyticalAccNode(space_charge_calculator,accNode.getName()+":"+str(part_index)+":")
	scNode.setLengthOfSC(rest_length)
	scNode.set_lattice_functions(lattice_functions_dict)
	scNodes_arr.append(scNode)
	accNode.addChildNode(scNode,AccNode.BODY,part_index,AccNode.BEFORE)
	
	return scNodes_arr
