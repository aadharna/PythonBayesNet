# CISC 6725 AI; Fordham Univ AI 2016
# 2/26/2017
# Aaron Dharna
# Inspired by Damian Lyons' implementation from AI-2016
#
# Assignment Five, Bayesian Networks
#
# The bayesian network is represented using an associative table BayesDict
# the CPT is also represented as an associative table within BayesDict
#
# BayesDict['NameofVariab'] is the node in the network
# BayesDict['NameofVariab']['cpt'] is the cpt
# cpt['Parent1']['Parent2']....['Parentn'] is how to look up the cpt
#
# if A and B are the random variables in the problem,
# a state of the world is indicated truth values for A and B. For example ['A','B'] is
# a state with both true, ['A','nB'] has A true and B false, etc. so 'n' is prefixed to
# the name to indicate the variable is false, similar to the convention in the [Russel & Norvig's AI].
#
# ------------------------------
#
# IF the BayesNet is of the correct graph structure, then Bayesian Inference will be performed
# on the class nodes.
#
# ------------------------------

import random
import numpy as np

def readBayesFiles(fileName, variableNames, BayesDict):

    Bayesfile = open(fileName, "r")

    BayesfileLines = Bayesfile.readlines()
    #print BayesfileLines

    for each_line in BayesfileLines:
        linelist = each_line.split()

        if linelist[0] == 'END':
            print "file finished"
            assignChildren(variableNames, BayesDict)
            return

        variableNames.append(linelist[0])
        readNode(linelist, BayesDict)

    Bayesfile.close()
    return

def readNode(Bayesline, nodetab):
    '''read bayesian network node from line of file'''

    nodename=Bayesline.pop(0)
    print "Reading node ",nodename
    print  Bayesline

    nodetab[nodename]={}
    parent=Bayesline.pop(0)
    #print "parent: ", parent

    nodetab[nodename]['children'] = []
    nodetab[nodename]['permutations'] = []

    if (parent == 'NONE'):
        nodetab[nodename]['numparents'] = 0
        nodetab[nodename]['prob']=float(Bayesline.pop(0))
        nodetab[nodename]['parents'] = parent
        return

    #figure out how many parents the node has
    nodetab[nodename]['parents'] = [parent]
    while (len(Bayesline)>0):
        if ( not Bayesline[0].isalpha()):
            break

        parent = Bayesline.pop(0)
        nodetab[nodename]['parents'].append(parent)
        #print "parent: ", parent
        #print nodetab[parent]
#    print "All parent list is ",nodetab[nodename]['parents']

    nodetab[nodename]['numparents']=len(nodetab[nodename]['parents'])

    #read in the CPT with this many parents
    nodetab[nodename]['cpt']={}
    ttread(nodetab[nodename]['cpt'], nodetab[nodename]['parents'], '', '', Bayesline, nodetab)
    return

def assignChildren(variableNames, BayesDict):
    #print "assigning children fn"
    for graphNode in variableNames:
        for each_attribute in BayesDict[graphNode]:
            #print graphNode, each_attribute
            if (each_attribute == 'parents'):
                for nodenames in variableNames:
                    if (nodenames in BayesDict[graphNode][each_attribute]):
                        BayesDict[nodenames]['children'].append(graphNode)
                        BayesDict[nodenames]['permutations'].append((graphNode, 'n'+graphNode))

    return


def ttread(parentcpt_dict, Randomvars, str, name, shrinkingCpt, nodetab):
    '''read a cpt from a file into the bayesian network'''
    #this if tests for the completion of a row in the CPT

    global condensed_cptlines

    if (len(Randomvars)==0):

        newval=float( shrinkingCpt.pop(0) ) # read the value
        parentcpt_dict[name]=newval
        condensed_cptlines += 1
        return

    newvars=list(Randomvars)
    x=newvars.pop(0)
    #this IF checks making nested associative list or reading a value into the list
    #which is done on the last parent of cpt[parent1]....[parentn]

    if (len(newvars)>0):
        parentcpt_dict[x]={}
        parentcpt_dict['n'+x]={} # make a nested associative list
        ttread(parentcpt_dict[x], newvars, str+x, x, shrinkingCpt, nodetab) #positive parent values
        ttread(parentcpt_dict['n'+x], newvars, str+'n'+x, 'n'+x, shrinkingCpt, nodetab) #negative parent values
        return

    ttread(parentcpt_dict, newvars, str+x, x, shrinkingCpt, nodetab)
    ttread(parentcpt_dict, newvars, str+'n'+x, 'n'+x, shrinkingCpt, nodetab)
    return


def ttlist(Randomvars, row, openfile, Bayesdict, combo_list):
    '''list all rows of the joint distribution of vars'''
    global totalcpt_lines


    #this if tests for the completion of a row
    if len(Randomvars)==0:

        prob=bayeseval(row, Bayesdict)
 #       print linenum,"Prob(",row,") = ",prob
        totalcpt_lines += 1
        combo_list.append(row)

        if (not openfile==0):
            openfile.write(str(row)+" "+str(prob)+"\n")

        return

    newvars=list(Randomvars)
    x = newvars.pop(0)
    r = list(row)
    r.append(x)

    ttlist(newvars, r, openfile, Bayesdict, combo_list) # positive case for x
    rn = list(row)
    rn.append('n'+x)
    ttlist(newvars, rn, openfile, Bayesdict, combo_list) # negative case for x
    return

def bayeseval(row, Bayesdict):
    '''evaluate case of bayesian network'''

    prob=1.0;

    for varnode in row:

        if varnode[0]=='n':
            nodename=varnode[1:len(varnode)]
            prob = prob * (1-bayesnodeeval(nodename, row, Bayesdict))
        else:
            nodename=varnode
            prob = prob * bayesnodeeval(nodename, row, Bayesdict)

    return prob

def bayesnodeeval(nodename, row, nodetab):
    '''evaluate single node in bayesian network'''

    if (nodetab[nodename]['numparents']==0):
        return nodetab[nodename]['prob']

    cpt = nodetab[nodename]['cpt']
    parentlist = nodetab[nodename]['parents']
    vlist=[]

    #extract the case for only the parents of the nodes
    for parent in row:
        if truename(parent) in parentlist:
            vlist.append(parent)
    return cpteval(cpt, vlist)

#returns the name of the passed node after cutting off the 'n'
def truename(parentName):
    if parentName[0] == 'n':
        return parentName[1:len(parentName)]

    return parentName

#DFS of the BayesDict
def cpteval(initialCpt, variableList):

    if (len(variableList) == 1):
        return initialCpt[variableList[0]]

    #remove the failed node
    name = variableList.pop(0)
    #keep searching but now from one node deeper
    return(cpteval(initialCpt[name], variableList))

def inferenceWrapper(varlist, BayesDict, combo_list, list):

    attributes = []

    for situations in combo_list:
        parentNode = situations.pop(0)
        temp = truename(parentNode)
        if BayesDict[temp]['numparents'] == 0:
            if len(BayesDict[temp]['children']) == (len(varlist) - 1):
                attributes = situations
                ans = Probability(parentNode, attributes, BayesDict)
                list.append(ans)
                print 'Probabilty({}, given:{}) -- {}'.format(parentNode, attributes, ans)
                situations.insert(0, parentNode)
            else:
                print 'Your graph is not of the correct format to perform inference.'
                return 0

    return 1

def Probability(node, children, BayesDict):

    prob = 1.0
    temp1 = 0.0
    temp2 = 0.0
    print
    for each_attribute in children:
        #print each_attribute
        if each_attribute[0] == 'n':
            temp1 = 1 - BayesDict[truename(each_attribute)]['cpt'][node]
        else:
            temp1 = BayesDict[each_attribute]['cpt'][node]

        #print temp1
        prob = prob * temp1
        #print 'prob: {}.'.format(prob)


    if node[0] == 'n':
        temp2 = 1 - BayesDict[truename(node)]['prob']
    else:
        temp2 = BayesDict[node]['prob']
    #print temp2
    prob = prob * temp2

    prob = float(prob) #/ float(sum)
    #print 'prob: ', prob

    return prob

def evaljointBayes(fname):
    global totalcpt_lines, condensed_cptlines
    totalcpt_lines = 0
    condensed_cptlines = 0

    BayesDict = {}
    varlist = []
    combo_list = []
    list = []
    MAP_class = 0
    readBayesFiles(fname, varlist, BayesDict)
    print BayesDict
    print varlist

    openfile = open(fname+'.joint.txt','w')
    ttlist(varlist, [], openfile, BayesDict, combo_list)
    openfile.close()

    #Yes, this causes the function to be run more than once.
    #HOWEVER, it is extremely easy to check if the right conditions are met
    #When they are not met is when this case comes back 0.
    #Therefore calling the function here does not cause double work to be done.
    if inferenceWrapper(varlist, BayesDict, combo_list, list) == 0:
        MAP_class = 'n/a'
    else:
        inferenceWrapper(varlist, BayesDict, combo_list, list)
        MAP_class = (combo_list[np.argmax(list)][0], list[np.argmax(list)])
    #print list

    print "------------------------------------------------------"
    print "Bayesian Network                  : ",fname
    print "Num of lines in joint distribution: ",totalcpt_lines
    print "Num of CPT lines                  : ",condensed_cptlines
    print "Compactness                       : ",float(condensed_cptlines)/float(totalcpt_lines)
    print "Num multiply operations           : ",len(varlist)*totalcpt_lines
    print "Num add operations                : ",float((len(varlist)*totalcpt_lines))/float(2)
    print "The MAP class (if appliciable) is : ",MAP_class
    print "------------------------------------------------------"


ans = raw_input("Enter the BayesFile. e.g. Bayes1.txt ")
evaljointBayes(ans)
