#!/usr/bin/env python

class Mammals:
    def __init__(self):
        """Constructor for this class."""
        self.members = ["Tiger", "Elephant", "Wild Cat"]

    def printMembers(self):
        print("Printing members of the Mammals class")
        for member in self.members:
            print(f"\t{member}")
