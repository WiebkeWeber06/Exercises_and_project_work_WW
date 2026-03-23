#!/usr/bin/env python

class Birds:
    def __init__(self):
        """Constructor for this class."""
        self.members = ["Sparrow", "Robin", "Duck"]

    def printMembers(self):
        print("Printing members of the Birds class")
        for member in self.members:
            print(f"\t{member}")

