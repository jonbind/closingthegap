# closingthebridge
Comprehensive data and code repository for the study on complexation of small interfering RNAs (siRNAs) with polyethyleneimine. This work aims to bridge experimental observations with simulation insights, providing a holistic understanding of the process. Includes data, analysis scripts and simulation parameters.

# How to use
Just download and run the ThePolymerGenerator.py file. It will output six files: file.gro, file.itp,file.map,file.mapBead,file.smiles,filepH.itp,fileTitration.itp

* **file.gro** : Standard .gro file of your created PEI molecule

+ **file.itp** : Topology file of your created PEI molecule

- **file.map** and file.mapBead : Mapping files of the coarse grained structure
  
* **file.smiles** : SMILES expression of the created molecule
  
* **filepH.itp** : Topology files used for the complexation studies. Charges are assigned according to the chosen pH in the script
  
* **fileTitration.itp** : Topology files used for titratable Martini simulations

# Part 1: Parameterising angle and bond distribution

# Part 2: Simulation of titratableMartini setups

# Part 3: Complexation of siRNA with PEI
