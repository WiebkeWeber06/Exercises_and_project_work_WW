
# Imports
# ----------------------------------------------

# imports the die.py file and allows for usage of class "Die"
from .die import Die

# imports the useless utils.py file
#from .utils import i_just_throw_an_exception

# imports the sys package, which allows us to use sys.exit()
import sys


# Define class GameRunner
# ----------------------------------------------

# contains the program which is run for the game

# The flow is:
# 1) create a GameRunner
# 2) create 5 dice
# 3) show the dice
# 4) ask for a guess
# 5) compare the guess to the sum
# 6) count consecutive correct answers with c
# 7) after each round, ask whether to play again
# 8) stop if the player enters n
# 9) win if c == 6


class GameRunner:

    def __init__(self):
        self.dice = Die.create_dice(5)
        self.reset()

    def reset(self):
        self.round = 1
        self.wins = 0
        self.losses = 0
        
# we got problems here I think
# replace the 1 with die.value
    def answer(self):
        return sum(die.value for die in self.dice)

    @classmethod
    def run(cls):
        # move runner outside the loop because otherwise,
        # it is impossible to win
        runner = cls()
        c = 0
        
        while True:
            print(f"Round {runner.round} \n")
            
            for die in runner.dice:
                die.roll()
                print(die.show())

            # Get user input
            while True:

                guess = input("What is the sum?: ")
                
            # Exclude everything which is not integer    
                try:
                    guess = int(guess)
                    print("\n")
                except ValueError:
                    print("Invalid input, please try again")
                    continue
                
            # Exclude numbers which are out of range    
                if guess not in range(5, 31):
                    print("Value can't be more than 30 and not below 5 "
                          + "with 5 dice, " +
                          "the math is not mathing")
                    continue
                
            # continue with game if correct input    
                else:
                    break
                    
            # guess is right    
            if guess == runner.answer():
                print("Congrats, you can add like a 5 year old... \n")
                runner.wins += 1
                c += 1
            
            # guess is wrong
            else:
                print("Sorry that's wrong")
                print("The answer is: {}".format(runner.answer()))
                print("Like seriously, how could you mess that up \n")
                runner.losses += 1
                c = 0
                
            print("Wins: {}/6, Losses: {}".format(runner.wins, runner.loses))
            
            # 6 consecutive rounds have been won
            if c == 6:
                print("You won... Congrats...")
                print("The fact it took you so long is pretty sad")
                break
            
# What is the sense of this?
# Change input text from Y -> y
# asks if one wants to continue the round

            while True:
                print(f"Finished round {runner.round} \n")
                runner.round += 1
                prompt = input(f"Would you like to continue with round {runner.round}? [y/n]: ")
        # I think we got problems here too
        # if yes, then exit this loop and go back to the first loop
                if prompt == 'y':
                    break
        # if no, then terminate the whole program        
                elif prompt == 'n':
                    print("You chose to give up already. What was I expecting?")
        # add break to terminate the program
                    sys.exit()
        # add a case for typo, so we repeat the loop
        # and request a different input            
                else:
                    print("Wrong input, did you do a typo? Please try again.")
                    continue                    
