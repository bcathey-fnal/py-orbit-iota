"""
Simple test script to check that MPI is functioning properly. Calling the
script should print the number of nodes.
"""

import orbit_mpi # Proper node identification for parallel operations

def get_mpi_rank():
    """
    Get total number of nodes and the rank of this MPI node.

    Parameters: None

    Returns:
    - node_mpi_rank: The rank of this node.
    - mpi_size: Total number of nodes in this job.
    """
    node_mpi_rank = 0  # By default this is primary node
    mpi_size = 1 # By default there is only one node
    mpi_init = orbit_mpi.MPI_Initialized() # Is MPI initialized
    comm = orbit_mpi.mpi_comm.MPI_COMM_WORLD # Get the communication world
    if mpi_init:
        node_mpi_rank = orbit_mpi.MPI_Comm_rank(comm)
        mpi_size = orbit_mpi.MPI_Comm_size(comm)
    return node_mpi_rank, mpi_size # Return the rank and size

if __name__ == "__main__":
    node_mpi_rank, mpi_size = get_mpi_rank()
    if node_mpi_rank == 0:
        print("PyORBIT started on {} nodes.".format(mpi_size))

