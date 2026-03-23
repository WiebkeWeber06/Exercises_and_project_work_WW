#!/usr/bin/env python

class Fish:
    def __init__(self):
        """Constructor for this class."""
        self.members = ["Shark", "Salmon", "Goldfish"]

    def printMembers(self):
        print("Printing members of the Fish class")
        for member in self.members:
            print(f"\t{member}")
