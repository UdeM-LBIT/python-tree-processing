from LabelGTCRec import LabelGTC

import os
import sys

import copy

from ..TreeLib import *
from ..TreeLib import TreeUtils, TreeClass

s = TreeClass("((A,B),C);")

print("\n\n\n\n")
print("**********************************************************************************")
print("*                                       INPUT                                    *")
print("**********************************************************************************")
print("\n\n\n\n")

seuil = 0.6

print("------THRESHOLD------")
print("\n")
print(seuil)
print("\n\n\n\n")

print("-----SPECIES TREE-----")
print("\n")
print(s)
print("\n")
print("\n\n\n\n")


g = TreeClass("((((a_A,x_B)0.8,(b_B,e_C)0.8),y_C)0.1,((i_B,k_A)0.1,((c_C, j_A)0.1,(d_C,(g_B,h_A)0.8)0.1)0.1)0.1)0.8;")

cst = [TreeClass("(((a_A,x_B),(b_B,e_C)),y_C);"), TreeClass("(c_C,j_A);"), TreeClass("(d_C,(g_B,h_A));"), TreeClass("(i_B,k_A);")]

print("-----COVERING SET OF TREES-----")
print("\n")
for tree in cst:
    print(tree)
    print("\n")
print("\n\n\n\n")

lgtc = LabelGTC(s, g, cst, seuil)

print("-----GENES TREE-----")
print("\n")
print(lgtc.getGenesTree().get_ascii(show_internal=True))
print("\n")
print("\n\n\n\n")

print("**********************************************************************************")
print("*                                 ALGORITHM BEGIN                                *")
print("**********************************************************************************")
print("\n\n\n\n")

print("CURRENT TREE BEING ANALIZED")
print(lgtc.getGenesTree().get_ascii(show_internal=True))
print("\n")

lgtc.mergeResolutions()
print("**********************************************************************************")
print("*                            ALGORITHM ENDED SUCCESFULLY                         *")
print("**********************************************************************************")
print("\n\n\n\n")