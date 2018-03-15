## CISC 6725 AI; Fordham Univ AI 2016

 Aaron Dharna
 
 Adapted and Expanded from Damian Lyons' implementation from AI-2016

#### Assignment Five, Bayesian Networks

----
##### Storing the graph
 The bayesian network is represented using _associative table_ table BayesDict
 the CPT is also represented as an _associative table_ within BayesDict

 BayesDict['NameofVariab'] is the node in the network
 
 BayesDict['NameofVariab']['cpt'] is the cpt
 
 cpt['Parent1']['Parent2']...['Parentn'] is how to look up the cpt
 
----------

##### Accessing the network
 if A and B are the random variables in the problem,
 a state of the world is indicated truth values for A and B. For example ['A','B'] is
 a state with both true, ['A','nB'] has A true and B false, etc. so 'n' is prefixed to
 the name to indicate the variable is false, similar to the convention in the [Russel & Norvig's AI].

 ------------------------------
##### NB:
 IF the BayesNet is of the correct graph structure, then Bayesian Inference will be performed on the class nodes.

 ------------------------------
