# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 16:07:23 2026

@author: wiewe372
"""

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(f"Hello from rank {rank} out of {size}")