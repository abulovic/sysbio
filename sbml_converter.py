#!/usr/bin/python

#tryout the libsbml parser to see how difficult it is easy to get user input from sbml

from libsbml import *
import sys

if __name__ == "__main__":
    reader = SBMLReader()
    document = reader.readSBML("")
    if document.getNumErrors() != 0:
        print "Something went wrong!"
        sys.exit()
