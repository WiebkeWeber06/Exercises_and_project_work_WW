# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:20:32 2026

@author: wiewe372
"""

class Person:
    def __init__(self, firstname, lastname):
        self.firstname = firstname
        self.lastname = lastname

    def get_full_name(self):
        return f"{self.firstname} {self.lastname}"
    
class Student(Person):
    def __init__(self, firstname, lastname, subject):
        super().__init__(firstname, lastname)  # initialize Person
        self.subject = subject

    def printNameSubject(self):
        print(f"{self.get_full_name()}, {self.subject}")    
        
        
me = Student('Benedikt', 'Daurer', 'physics')
me.printNameSubject()        


class Teacher(Person):
    def __init__(self, firstname, lastname, course):
        super().__init__(firstname, lastname)
        self.course = course

    def printNameCourse(self):
        print(f"{self.get_full_name()}, {self.course}")
        
teacher = Teacher('Jan', 'Dalton', 'Bacteria in Yellowstone')
teacher.printNameCourse()        