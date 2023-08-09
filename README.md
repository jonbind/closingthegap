# Closing the bridge
Comprehensive data and code repository for the study on complexation of small interfering RNAs (siRNAs) with polyethyleneimine. This work aims to bridge experimental observations with simulation insights, providing a holistic understanding of the process. Includes data, analysis scripts and simulation parameters.

![Alt text](./README/complex.png)

# How to use
Download and run the Polymer Generator using this command 
	python3 ThePolymerGenerator.py
 
* The program starts with a linear bead in tree level 1, which has 2 open binding options.
* Bindings are performed depending on the pre-defined probabilities of terminatingProbability and degreeOfBranching. The two beads bound to it are assigned to tree level 2.
* Before each new tree level, the maximum polymer mass is compared with the current polymer mass (-> remainingMass). In addition, the number of currently open binding options is multiplied by the mass of the terminating bead (-> terminatingMass) to determine if polymerization termination needs to be initiated. If this is the case, all remaining binding options are connected to terminating beads and the polymerization ends.

	if (remainingMass < terminatingMass):
	initiate polymerization termination
	else:
	proceed to the next tree level

* The open binding options of the two new beads are again occupied depending on the probabilities and the new beads are assigned to tree level 3. This process continues successively until termination is initiated according to the above-mentioned condition.
* After polymerization termination, some calculations are performed on the beads, such as coordinates, bond lengths, and angles.
* Finally, it is checked whether there are directly overlapping beads in all 3 spatial dimensions (i.e., same x, y AND z). If this is the case, the polymerization starts again, otherwise it is successfully completed.
* The SMILES string and its corresponding mapping are then calculated.


* **file.gro** : Standard .gro file of your created PEI molecule

+ **file.itp** : Topology file of your created PEI molecule

- **file.map** and file.mapBead : Mapping files of the coarse grained structure
  
* **file.smiles** : SMILES expression of the created molecule
  
* **filepH.itp** : Topology files used for the complexation studies. Charges are assigned according to the chosen pH in the script
  
* **fileTitration.itp** : Topology files used for titratable Martini simulations

# Part 1: Parameterising angle and bond distribution
Here all the simulation files are bundled together so feel free to test it out

# Part 2: Simulation of titratableMartini setups
This parts needs further balall

# Part 3: Complexation of siRNA with PEI
