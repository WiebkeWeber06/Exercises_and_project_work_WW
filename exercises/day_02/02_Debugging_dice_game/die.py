# Imports
# ----------------------------------------------

# Imports the random module
import random

# ----------------------------------------------
# This function seems useless
# it is not needed, so I removed it
#def roll(dice):
#    for die in dice:
        # XXX: I don't even know what this function does
#        continue



# This function is for creating dice with random values 
# ----------------------------------------------
  
# class definition
# represents a single die
# When you create a die, it automatically 
# 1) rolls itself
# 2) gets a value


class Die:
    """
    This is always correct. Seriously, look away.
    """
    # Constructor
    def __init__(self):
        self.roll()
        
    # assign random value between 1 and 6
    def roll(self):
        self.value = int(random.random() * 6 + 1)
        
    # print the die as its value as a drawing
    def show(self):
        if self.value == 1:
            return("---------\n|       |\n|   *   |\n|       |\n---------")
        elif self.value == 2:
            return("---------\n|*      |\n|       |\n|      *|\n---------")
        elif self.value == 3:
            return("---------\n|*      |\n|   *   |\n|      *|\n---------")
        elif self.value == 4:
            return("---------\n|*     *|\n|       |\n|*     *|\n---------")
        elif self.value == 5:
            return("---------\n|*     *|\n|   *   |\n|*     *|\n---------")
        else:
            return("---------\n|*     *|\n|*     *|\n|*     *|\n---------")
    
    # create multiple dice (n die objects)
    @classmethod
    def create_dice(cls, n):
        return [cls() for _ in range(n)]

# Why use a class for a die?
# Because everything related to a die lives in one place
# with using functions, value and behavior would be seperate and 
# things would have to be moved around manually