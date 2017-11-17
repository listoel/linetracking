# -*- coding: utf-8 -*-

########################################################################
#                                                                      #
#       Other tools for MAD-X-like tracking in python.                 #
#       Author: L.S. Stoel                                             #
#       Version 0.1 - Work in progress                                 #
#                                                                      #
########################################################################

##################################################
#                                                #
#   lineprint                                    #
#                                                #
##################################################

def lineprint(line):
    dist = 0.0
    for element in line:
        dist += element.len
        print(element.name, '{0:.4f}'.format(dist))
    return
