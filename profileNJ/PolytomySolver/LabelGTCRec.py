"""
Author: Kylian Berrebi
Date: 06/2017
Functions for the LabelGTC resolution
"""

"""
Expected input:
    - A species tree
    - A gene tree
    - A covering set of trees (list of TreeClass trees)
    - A support for each (non-leaf) node that represents the confidence in the subtree rooted at this node
    - A threshold value (from 0 to 1) from wich the subtree architecture should be trusted
"""

import os
import sys

from sys import stdout

import re

import time

import copy

from PolySolver import PolytomySolver
from PolySolver import GeneTreeSolver

from ..TreeLib import *
from ..TreeLib import TreeUtils, TreeClass



class LabelGTC:
    
    """
    Usage example :
    
    speciesTree = TreeClass("...")
    genesTree = TreeClass("...")
    coveringSetTree = [TreeClass("..."), TreeClass("..."), ...]
    threshold = 0.8
    
    lgtc = LabelGTC(speciesTree, genesTree, coveringSetTree, threshold)
    
    lgtc.mergeResolutions()
    """
    
    
    def __init__(self, speciesTree, genesTree, covSetTree, threshold):
        
        self.speciesTree = speciesTree
        
        #self.speciesTree.label_internal_node()
        
        self.genesTree = genesTree
        
        self.genesTree.label_internal_node()
        
        self.covSetTree = covSetTree
        self.threshold = threshold
        
        
     
    def getSpeciesTree(self):
        return self.speciesTree
        
    
    def getGenesTree(self):
        return self.genesTree
        
    
    def setThreshold(self, newThreshold):
        self.threshold = newThreshold
        

    
    def checkCovSetTree(self):
        """check if the covering set of tree is contained by the tree of genes, and add the label cst to each node of the tree of genes that is also a node of a tree from the covering set of tree (except the root node of each tree) """
        
        cpt = 0
        #list containing all the leaves of each tree from the set of cevering trees
        leaves_list_cst = []
        
        #list of all the leaves of the tree of genes
        leaves_list_gt = self.genesTree.get_leaf_names()
            
        
        for subtree in self.covSetTree:
            
            if subtree.is_leaf():
                leaves_list_cst.append(subtree.name)
            
            else:
                leaves_list_cst.extend(subtree.get_leaf_names())
            
            #searching in the tree of genes a subtree identical to the current subtree
            for g_node in self.genesTree.traverse("postorder"):
                
                if not g_node.has_feature('cst'):
                    g_node.add_features(cst=0)
                
                #testing if the tree rooted at the current node has the same topology than the current subtree
                if g_node.has_same_topo(subtree):
                    cpt += 1
                    if g_node.cst == 0:
                        g_node.add_features(cst=2)
                    
                    #add the cst feature to all descendants of the current node
                    for child in g_node.get_descendants():
                        child.add_features(cst=1)
                        
        #all the subtrees have been found in the tree of genes
        return (cpt==len(self.covSetTree) and len(leaves_list_cst) == len(leaves_list_gt) and set(leaves_list_cst)==set(leaves_list_gt))
    
    

    def binaryLabeling(self):
        """Binarization of the support for each node according to the threshold"""
        
        #checking if the covering set of tree is conform with the tree of genes
        if self.checkCovSetTree():
        
            for g_node in self.genesTree.traverse("levelorder"):
                
                    #adding the binary feature
                    if g_node.support >= self.threshold:
                        g_node.add_features(binconfidence = 1)
                    else:
                        g_node.add_features(binconfidence = 0)
        
        else:
            raise Exception("The covering set of tree is not conform with the tree of genes")
            
    
    
    def largerCSE(self):
        """Identify the larger covering set of edges such that each edge of this set has no ancestral edge labelled 1, by adding lcse feature to the concerned nodes"""
        
        #To know when the loop should be stopped TODO (attribute)
        leaves_to_compare = []
        all_leaves = self.genesTree.get_leaf_names()
       
        for g_node in self.genesTree.traverse("levelorder"):
            
            children = g_node.get_children()
            
            #Searching for potential children for the set, and adding to them a feature indicating that they are in the set
            for child in children:
                
                if (child.binconfidence == 1 and child.cst != 1):
                    child.add_features(lcse=1)
                    leaves_to_compare += child.get_leaf_names()
        
                elif child.cst == 2:
                    child.add_features(lcse=1)
                    leaves_to_compare += child.get_leaf_names()
            
            #Stopping the loop when the larger covering set of edges is found        
            if set(all_leaves).issubset(set(leaves_to_compare)):
                break
        
        
    
    def globalProcessing(self):
        """Calling LabelGTC recursively on concerned subtrees of the tree of genes (processing the global case)"""
        
        print("\n\n")
        print("************************************************* NEW INSTANCE **********************************************")
        print("\n\n")
        
        #Computing first the larger covering set of tree
        self.largerCSE()
        
        #List of leaves to compare with sub_leaves in order to know when the loop should be stopped TODO (attribute)
        toCompare_leaves = self.genesTree.get_leaf_names()
        sub_leaves = []
        
        #Storing the trees that will be computed recursively
        rec_trees = []
        
        cpt = 0
        cpt_bis = 0
              
        #Applying recursively the LabelGTC algorithm
        for g_node in self.genesTree.traverse("levelorder"):
            
            #Testing if all the covering set of edges has been checked
            if sub_leaves == toCompare_leaves:
                break
            
            #Testing only the subtrees of the covering set of edges
            if g_node.has_feature('lcse') and (not g_node.is_root()) and (not g_node.has_feature('root')):
                sub_leaves += g_node.get_leaf_names()
                
                #Applying recursively the LabelGTC algorithm to the nodes that have a 1 confidence in the covering set of edges
                if g_node.binconfidence == 1 and g_node.cst == 0:
                    #The covering set of tree for the new instance of LabelGTC
                    cst_subtree = []
                    
                    #Building the covering set of tree for the new instance
                    for tree in self.covSetTree:
                        
                        tree_leaf_names = tree.get_leaf_names()
                        g_node_leaf_name = g_node.get_leaf_names()
                        included = True
                        
                        for element in tree_leaf_names:
                            
                            if not element in g_node_leaf_name:
                                included = False
                                break
                                
                        if included:
                            cst_subtree.append(tree)
                            
                    print("CURRENT TREE BEING ANALIZED:")
                    print(g_node.get_ascii(show_internal=True))
                    print("\n")
                    #This node is the root of the tree of genes of the new instance      
                    g_node.add_features(root=1)
                    g_node.detach()
                    #New instance with the current subtree and the reduced covering set of tree (limited to the subtree)
                    lgtc = LabelGTC(self.speciesTree, g_node, cst_subtree, self.threshold)
                    lgtc.mergeResolutions()
                    
        
        #Creating a reduced version of the tree of genes by contracting each subtree that has a root node with label lcse and a binconfidence equals to 1, to a single leaf node
        for g_node in self.genesTree.traverse("levelorder"):
            
            cpt += 1
            
            if g_node.has_feature('lcse') and cpt > 1 and not g_node.has_feature('treated') and not g_node.is_leaf():
                built = self.genesTree&g_node.name
                #Storing the subtree to attach it back later when modified
                rec_trees.append(built)
                rebuilt = TreeClass()
                built.insert_node_between(built.up, rebuilt)
                rebuilt.add_features(root=1)
                rebuilt.add_features(treated=1)
                built.up.name = built.name 
                built.detach()
                
        print("\n\n")        
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("\n")           
        print("Using minSGT algorithm on the current node")
        print("\n")
        print("-----------CONTRACTED TREE-----------")              
        print(self.genesTree.get_ascii(show_internal=True))
        print("-------------------------------------")
        print("\n")
        
        #Attaching back the previously contracted subtrees now resolved   
        for g_node in self.genesTree.traverse("levelorder"):
            
            if cpt_bis > len(rec_trees) - 1:
                break
                
            elif g_node.has_feature('treated'):
                #rec_trees[cpt_bis].add_features(binconfidence=0)
                tobuild = self.genesTree&g_node.name
                parent = tobuild.up
                rec_trees[cpt_bis].add_features(treated=1)
                parent.add_child(rec_trees[cpt_bis])
                tobuild.delete()       
                cpt_bis += 1
                
        print("-----------REBUILT TREE--------------")  
        print(self.genesTree.get_ascii(show_internal=True))
        print("-------------------------------------")
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")       
        print("\n\n\n\n")                  
    
    def polyRes(self):
        """Using PolytomySolver Algorithm"""
        
        
        dupcost = 1
        losscost = 1
        
        self.genesTree.set_species()
        
        lcamap = TreeUtils.lcaMapping(self.genesTree, self.speciesTree, multspeciename=False)
        print(self.speciesTree.write(features=[]))
        print(type(self.genesTree))
        print(self.genesTree)
        for node in self.genesTree.traverse("levelorder"):
            print(node.features)
        s_images = TreeUtils.getImageTreeNode(
            self.genesTree, self.speciesTree, self.genesTree.lcaMap)

        #print "COST=", ps.getTableValue(self.speciesTree, 1)
        
        gts = GeneTreeSolver(self.genesTree, self.speciesTree, lcamap, dupcost, losscost)
        
        gts.labelInternalNodes(self.genesTree)
        gts.labelInternalNodes(self.speciesTree)
        
        gts.use_dp = False
        if dupcost != losscost:
            gts.use_dp = True

        r = [x + ";" for x in gts.solvePolytomies(1000)]

        print "NBSOLS=", len(r)
        print r
    
    
    def mergeResolutions(self):
        """Using the different kind of resulutions according to the labelling of the tree of genes"""
        
        onlyLeaves = True
        polyResCompatible = True
        minTRSCompatible = True
        minSGTCompatible = True
        cpt = 0
        
        self.binaryLabeling()
        
        #Testing the case where the covering set of trees is only composed by leaves
        for subtree in self.covSetTree:
            if not subtree.is_leaf():
                onlyLeaves = False
                break
        
        if onlyLeaves:
            print("-------> Using M-PolyRes algorithm")
        
        #Testing global case
        else:
            
            for g_node in self.genesTree.traverse("levelorder"):
                
                if not (g_node.is_root() or g_node.has_feature('root')):
                #Processing non-terminal edges
                    if g_node.cst == 0:
                    
                        if g_node.binconfidence == 1:
                            polyResCompatible = False
                            minSGTCompatible = False
                            cpt += 1
                            
                        elif g_node.binconfidence == 0:
                            minTRSCompatible = False
                            cpt += 1
                            
                    #Processing terminal edges    
                    elif g_node.cst == 2:
                    
                        if g_node.binconfidence == 0:
                            polyResCompatible = False
                            cpt += 1
                            
                        elif g_node.binconfidence == 1:
                            minTRSCompatible = False
                            minSGTCompatible = False
                            cpt += 1
                            
                    if (not polyResCompatible) and (not minTRSCompatible) and (not minSGTCompatible):
                        break
        
            if polyResCompatible:
                print("-------> Using polyRes algorithm")
                self.polyRes()
            
            if minTRSCompatible and cpt > 2:
                print("-------> Using minTRS algorithm")
                
            if minSGTCompatible:
                print("-------> Using minSGT algorithm")
            
            if not (polyResCompatible or minTRSCompatible or minSGTCompatible):
                print(" -------> Using globalProcessing")
                self.globalProcessing()
                
           
                
            
                
        
        
        
        



    