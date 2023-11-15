# Closing the gap
Comprehensive data and code repository for the study on complexation of small interfering RNAs (siRNAs) with polyethyleneimine. This work aims to bridge experimental observations with simulation insights, providing a holistic understanding of the process. Includes data, analysis scripts and simulation parameters.

![Alt text](./README/complex.png)

# How to use
Download and run the Polymer Generator using this command 

	python3 ThePolymerGenerator.py
 
The program starts with a linear bead in tree level 1, which has 2 open binding options.
Bindings are performed depending on the pre-defined probabilities of terminatingProbability and degreeOfBranching. The two beads bound to it are assigned to tree level 2.
Before each new tree level, the maximum polymer mass is compared with the current polymer mass (-> remainingMass). In addition, the number of currently open binding options is multiplied by the mass of the terminating bead (-> terminatingMass) to determine if polymerization termination needs to be initiated. If this is the case, all remaining binding options are connected to terminating beads and the polymerization ends.


  
    if (remainingMass < terminatingMass):  
    initiate polymerization termination  
    else:  
    proceed to the next tree level  

The open binding options of the two new beads are again occupied depending on the probabilities and the new beads are assigned to tree level 3. This process continues successively until termination is initiated according to the above-mentioned condition.
After polymerization termination, some calculations are performed on the beads, such as coordinates, bond lengths, and angles.
Finally, it is checked whether there are directly overlapping beads in all 3 spatial dimensions (i.e., same x, y AND z). If this is the case, the polymerization starts again, otherwise it is successfully completed.
The SMILES string and its corresponding mapping are then calculated.
The script completes after generating the following output files


* **file.gro** : Standard .gro file of your created PEI molecule

+ **file.itp** : Topology file of your created PEI molecule

- **file.map** and file.mapBead : Mapping files of the coarse grained structure
  
* **file.smiles** : SMILES expression of the created molecule
  
* **filepH.itp** : Topology files used for the complexation studies. Charges are assigned according to the chosen pH in the script
  
* **fileTitration.itp** : Topology files used for titratable Martini simulations

# Part 1: Parameterising angle and bond distribution
Following the guidelines provided by [Martini Website](http://cgmartini.nl/index.php/martini-3-tutorials/parameterizing-a-new-small-molecule), we undertook the parameterization of a small molecule. This process can be summarized in six essential steps all of which you can find in the subdirectory. In this example we deposited the bPEI with a 33% branching. Note: Due to file size limit of 100 MB some simulations have to be performed again to get all the necessary files. 


# Part 2: Simulation of titratableMartini setups
In alignment with the methodologies presented in the original [paper](https://pure.rug.nl/ws/portalfiles/portal/130104951/5.0014258.pdf) and the respective [tutorial](http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/556-11-titratable-martini), we expanded the scope of pH simulations to encompass molecules of varying sizes. Specifically, we examined the influence of distinct protonation states on the radius of gyration for bPEI. This folder contains the relevant files pertaining to the 5 kDa bPEI.


# Part 3: Complexation of siRNA with PEI
This section contains simulation data for the titration of a 1.5 kDa branched polyethylenimine (bPEI) molecule with a single siRNA molecule. Within this directory, you'll find a folder for the exemplary pH 5.5, along with six subfolders corresponding to different N/P ratios: 0.6, 1.2, 1.8, 3.6, 6.6, and 9.6. These simulations explore the interactions between bPEI and siRNA, and investigate the effects of varying N/P ratios on these interactions.
To run the simulations simply execute the following commands (gromacs has to be patched with PLUMED):
```bash
#!/bin/bash

# GROMACS preprocessing and simulation commands

# Minimization
gmx grompp -p top.top -c solvated.gro -f minimization.mdp -o minimization.tpr
gmx mdrun -deffnm minimization -v

# Equilibration
gmx grompp -p top.top -c minimization.gro -f equilibration.mdp -o equilibration.tpr -r minimization.gro -maxwarn 2 -n index.ndx
gmx mdrun -deffnm equilibration -v

# Adsorption
gmx grompp -p top.top -c equilibration.gro -f adsorption.mdp -o adsorption.tpr -maxwarn 1 -r equilibration.gro -n index.ndx
gmx mdrun -deffnm adsorption -v -plumed adsorption.dat

# Dynamic Simulation
gmx grompp -p top.top -c adsorption.gro -f dynamic.mdp -o dynamic.tpr -maxwarn 1 -r adsorption.gro -n index.ndx
gmx mdrun -deffnm dynamic -v -plumed plumed-rna.dat
```

#
**If you have specific questions regarding the replication of these simulations feel free to ask: jonas.binder@cup.uni-muenchen.de or benjamin.winkeljann@cup.uni-muenchen.de**
