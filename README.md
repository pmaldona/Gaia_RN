# Gaia_RN
Program that calculates the feasibility of a reaction network to being an organization.

The program consists of a main function, which receives as input the path of an sbml file. And it returns a list type file with information from the reaction network.

## Function description

``rn.proc(f=rn.name,tot_org=T)``

### Input

``rn.name`` path of .sbml file of the reaction network  
``tot_org`` Boolean that submitted the calculation of the feasibility of being an organization

### Output

List element components:

``$sp.id`` character vector with sbml file component id  
``$sp.idn`` character vector with sbml file component idn  
``$sp.name`` character vector with sbml file component name   
``$reac`` list of all reaction extracted form the *.sbml file with their respective reactive and product stoichiometry and components involved.  
``$mp`` productive part of the stoichiometric matrix  
``$mr`` reactive part of the stoichiometric matrix  
``$nsp`` components required to become an organization.  




