# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:44:14 2026

@author: wiewe372
"""

from mpi4py import MPI

# Initialize communicator
comm = MPI.COMM_WORLD

rank = comm.Get_rank()
size = comm.Get_size()

# Each process contributes its rank
local_value = rank

# Reduce all values to rank 0 (sum them)
total_sum = comm.reduce(local_value, op=MPI.SUM, root=0)

# Only rank 0 prints the result
if rank == 0:
    print(f"Sum of ranks (0 to {size-1}) is: {total_sum}")

