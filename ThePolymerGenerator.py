# -*- coding: utf-8 -*-
"""
@authors:
    code written by                                 Joshua Winkeljann 
    extensively challenged by                       Jonas Binder
    Project coordination by                         Benjamin Winkeljann
    
    
    First creation of the project on    Sep 06 2021
Last update (pH dependency) on      Sep 19 2022

Please feel free to use, adapt and modify the code for your own purposes.
We would be happy to be credited by citing:
    
    Binder et al. Coarse grained modelling of polyethylenimine using the Martini 3 forcefield


For questions simply send an email to joshi.winkeljann@gmail.com or jonas.binder@cup.uni-muenchen.de



With the PolymerGenerator you can create random polymers out of TN6d, SN4 and SN3a Beads
for further use of Coarse Grained Simulations with the Martini Software



The PolymerGenerator generates different output files, with the following endings:
    
    - .itp  (3 different version)
        - one standard file
        - one file for Titration purposes
        - one file for pH dependency
        
    - .gro
    - .map
    - .beadMap
    - .smiles
    
    
    
The coded is structured in 5 different parts, namely:
    
    - PART 1: Values to be defined by the user
    - PART 2: Classes and their corresponding functions
    - PART 3: Class independent functions
    - PART 4: Start of the Main Programm
    - PART 5: Creation of output files
    
    
       
For successful using of the PolymerGenerator consider the following aspects:
    
    - Right at the beginning of the code, the user can modify the different values
    
        - the maximumPolymerMass should be larger than 200 to generate useful polymers, but there is no upper limit known for the moment
        - the degreeOfBranching is a value between 0 and 1 as it defines a probability
        - same is true for the terminatingProbability
        
    - If you like a pH dependent generation, the pH value should be set either to 7.4 or 5.5 because therefore the corresponding 
    probabilities for charging are known (Otherwise these need to be set manually in PART5 creation of output files)
    
    - The values for the different Beads can also be modified
    Here there are a few things to consider: For the corresponding string values like name, size and typeOfBead changes should only be made if the
    user has a good overview about how the code works as there are a few hardcoded if/else conditions with these names implemented
    
    - Bond lengths between the different beads are hardcoded and can be modified in "PART3 Class independent functions" in the getBondLength() function
    
    - Same for Force constant and distance in writeDistanceAndForceConstant() function in "PART3 Class independent functions"
    and the equilibrium angle and force constant in the writeEquilibriumAngleAndForceConstant() function in "PART3 Class independent functions"
"""




import random
import math




''' PART 1: Values to be defined by the user ================================================================================================
========================================================================================================================================='''

# To be defined by the user

maximumPolymerMass = 5000

degreeOfBranching = 0.584

terminatingProbability = 0.2

# For creating a filepH.itp, the pH value needs to be set either 7.4 or 5.5
pH = 7.4

# User defined Beads: size options: tiny, small or regular
class TerminalBead():   
    def __init__(self):
        self.name = "TN6d"
        self.size = "tiny"
        self.mass = 30.05
        self.additionalBonds = 0
        self.typeOfBead = "terminal"

        self.amountOfN = 1
        self.amountOfC = 1
        self.amountOfH = 4
        
class LinearBead():
    def __init__(self):
        self.name = "SN4"
        self.size = "small"
        self.mass = 43.97
        self.additionalBonds = 1
        self.typeOfBead = "linear"

        self.amountOfN = 1
        self.amountOfC = 2
        self.amountOfH = 5

class DendriticBead():
    def __init__(self):
        self.name = "SN3a"
        self.size = "regular"
        self.mass = 56.09
        self.additionalBonds = 2
        self.typeOfBead = "dendritic"

        self.amountOfN = 1
        self.amountOfC = 3
        self.amountOfH = 6
        
# Dendritic = t, linear = s, terminal = p





''' PART 2: Classes and their corresponding functions =======================================================================================
========================================================================================================================================='''


# Bead Class with all the necessary information    
class Bead():
    
    def __init__(self,
                 name,
                 mass,
                 additionalBonds,
                 typeOfBead,
                 size,
                 identificationNumber = -1,
                 treeLevel = -1):
        
        self.name = name
        self.mass = mass
        self.typeOfBead = typeOfBead
        self.size = size
        self.additionalBonds = additionalBonds
        self.identificationNumber = identificationNumber
        self.remainingBonds = additionalBonds
        self.treeLevel = treeLevel


    # function to reduce the remaining bond amount    
    def makeBond(self):
        if self.remainingBonds > 0:
            self.remainingBonds -= 1
        
            
#-----------------------------------------------------------------------------      
   
# class to keep track of the current State of the polymerization   
class CurrentState():
    
    def __init__(self):
        self.currentPolymerMass = 0
        self.currentNumberOfBranches = 1
        self.beadList = []
        self.idList = []
        self.bondList = []
        self.angleList = []
        self.mappingList = []
        self.currentBeadNumber = 1
        
        self.amountOfTerminalBeads = 0
        self.amountOfLinearBeads = 0
        self.amountOfDendriticBeads = 0

        self.currentNumberOfC = 1
        self.currentNumberOfN = 1
        self.currentNumberOfH = 1
        
    # method to add a baet to the polymer
    def addBead(self, bead: Bead):
        
        # Reduce 1 branch because of added Bead
        self.currentNumberOfBranches -= 1
        
        # Add the additional Bonds to the resulting Polymer
        self.currentNumberOfBranches += bead.additionalBonds
        
        # Add Bead Mass to the currentMass
        self.currentPolymerMass += bead.mass
        
        
        
        
        # Add bead to the list and give it an ID, specify the type and size
        self.beadList.append(bead)
        self.idList.append(
            [bead.name,
             bead.identificationNumber,
             bead.typeOfBead,
             bead.size])
        
        # increase the bead number and therefore the id for the next bead
        self.currentBeadNumber += 1
        
             
    # use the bond list to find out the bond legths
    def appendBondLengths(self):
    
        # iterate through the list
        for i in range(len(self.bondList)):
            
            # find the ids of the bond partners
            partner1id = self.bondList[i][0]
            partner2id = self.bondList[i][1]
            
            # iterate through the idList to find the bead size to the id which defines the bond length
            for j in range(len(self.idList)):
                if partner1id == self.idList[j][1]:
                    partner1size = self.idList[j][3]
                
                if partner2id == self.idList[j][1]:
                    partner2size = self.idList[j][3]
                
            # append the bond lengths to the bond list
            self.bondList[i].append(getBondLengths(partner1size, partner2size))
    

    # method to count the residue
    def selectAndCountResidue(self, typeOfBead):
        
        if typeOfBead == "terminal":
            self.amountOfTerminalBeads += 1
            return "P"+str(self.amountOfTerminalBeads)
        
        elif typeOfBead == "linear":
            self.amountOfLinearBeads += 1
            return "S"+str(self.amountOfLinearBeads)
        
        elif typeOfBead == "dendritic":
            self.amountOfDendriticBeads += 1
            return "T"+str(self.amountOfDendriticBeads)

  
    # find out how many beads of each type are there
    def appendResidues(self):
        
        # iterate through the id list to count the different types of beads
        for i in range(len(self.idList)):
            self.idList[i].append(self.selectAndCountResidue(self.idList[i][2]))
            

    # method to create a second bond list to add the amount of bonds of each bead to the list       
    def createSecondBondList(self):
        # create new list to be filled
        newList = []
        
        # create helping variable
        count = 0
        
        # iterate through id list
        for i in range(len(self.idList)):
            
            # iterate through existing bond list
            for j in range(len(self.bondList)):
                
                if self.bondList[j][0] == i+1:     # i starts with 0, but first id is 1
                    count += 1
                    newList.append([self.bondList[j][0], self.bondList[j][1], count])
            
            # reset count variable for the next id number
            count = 0
        
        return newList
        

    def initializeFirst3BeadCoordinates(self):
        
        # append initial coordinates for the first 3 beads
        self.idList[0].append(0.0)
        self.idList[0].append(0.0)
        self.idList[0].append(0.0)
        
        
        # check if bead 2 got already terminated
        self.idList[1].append(getBondLengths(self.idList[1][3], self.idList[0][3]))
        self.idList[1].append(0.0)
        self.idList[1].append(0.0)
        
        
        # check if bead 3 got already terminated
        self.idList[2].append(- getBondLengths(self.idList[2][3], self.idList[0][3]))
        self.idList[2].append(0.0)
        self.idList[2].append(0.0)
        

    # calculate the coordinates
    def calculateCoordinates(self):
        
        # new list with bondCounts
        newBondList = self.createSecondBondList()
        
        #initialize the first three beads
        self.initializeFirst3BeadCoordinates()
        

        # iterate through every id
        for i in range(len(self.idList)):
            
            # already done manually
            if i <= 2:
                continue
            
            # find out first bond partner
            for j in range(len(self.bondList)):
                if self.bondList[j][1] == i+1:    # i starts with 0, but first id is 1
                    
                    partner1id = self.bondList[j][0]
                    
                    
                    # find out second partner
                    for y in range(len(self.bondList)):
                        
                        if partner1id == self.bondList[y][1]:
                            
                            partner2id = self.bondList[y][0]
                            
                            # get the coordinates of the former bond partners
                            partner2x = self.idList[partner2id - 1][5]
                            partner2y = self.idList[partner2id - 1][6]
                            partner2z = self.idList[partner2id - 1][7]
                            
                            partner1x = self.idList[partner1id - 1][5]
                            partner1y = self.idList[partner1id - 1][6]
                            partner1z = self.idList[partner1id - 1][7]
                            
                            
                            # find out bondNumber, because dendritics have a first and second bond which is important for the angle
                            bondNumber = newBondList[j][2]
                            
                            # find out angle (2D)
                            angle = getAngle(self.idList[partner1id - 1][2], bondNumber)   # id starts from 1, list index from 0!
                          
                            # get new Vector 
                            newXValue, newYValue, newZValue = calculateRotatedVector(partner1x - partner2x,
                                                                                     partner1y - partner2y,
                                                                                     partner1z - partner2z,
                                                                                     angle,
                                                                                     self.bondList[j][2])
                            
                            # Add the vector to the Coordinates of the bondPartner to get the new Coordinates
                            newXCoordinate = partner1x + newXValue
                            newYCoordinate = partner1y + newYValue
                            newZCoordinate = partner1z + newZValue
                            
                            self.idList[i].append(newXCoordinate)
                            self.idList[i].append(newYCoordinate)
                            self.idList[i].append(newZCoordinate)
                            
                            # get the 3D angle because of the small pertubation in z direction
                            angle = calculate3DAngle(
                                partner2x, partner2y, partner2z,
                                partner1x, partner1y, partner1z,
                                newXCoordinate, newYCoordinate, newZCoordinate)
                            
                            angleInDegree = 180 - angle*180/math.pi
                            
                            # append the information to the angle list
                            self.angleList.append([partner2id, partner1id, i+1, angleInDegree])
                            
            
    # creation of List for mapping file
    def createMappingList(self):
        # variable to count the number of all atoms
        numberOfAtoms = 0

        for bead in self.idList:
            typeOfBead = bead[2]

            # get the amount of N, C and H atoms for each bead
            amountOfN, amountOfC, amountOfH = self.getAmountOfAtoms(typeOfBead)

            # append the information to the mapping List depending on the number of atoms
            for i in range(amountOfN):
                self.mappingList.append([
                    numberOfAtoms,
                    "N"+str(self.currentNumberOfN),
                    bead[4]
                ])
                self.currentNumberOfN += 1
                numberOfAtoms += 1

            for i in range(amountOfC):
                self.mappingList.append([
                    numberOfAtoms,
                    "C" + str(self.currentNumberOfC),
                    bead[4]
                ])
                self.currentNumberOfC += 1
                numberOfAtoms += 1

            for i in range(amountOfH):
                self.mappingList.append([
                    numberOfAtoms,
                    "H" + str(self.currentNumberOfH),
                    bead[4]
                ])
                self.currentNumberOfH += 1
                numberOfAtoms += 1


    # function to return number of individual atoms
    def getAmountOfAtoms(self, typeOfBead):
        if typeOfBead == "terminal":
            amountOfN = TerminalBead().amountOfN
            amountOfC = TerminalBead().amountOfC
            amountOfH = TerminalBead().amountOfH

        elif typeOfBead == "linear":
            amountOfN = LinearBead().amountOfN
            amountOfC = LinearBead().amountOfC
            amountOfH = LinearBead().amountOfH

        elif typeOfBead == "dendritic":
            amountOfN = DendriticBead().amountOfN
            amountOfC = DendriticBead().amountOfC
            amountOfH = DendriticBead().amountOfH

        else:
            print("There's an error! Type of Bead seems to be wrong!")

        return amountOfN, amountOfC, amountOfH





    # function for debugging purposes
    def printPolymer(self):
        print(self.currentPolymerMass)
        print(self.idList)
        print(self.bondList)


#=============================================================================



''' PART 3: Class independent functions =====================================================================================================
========================================================================================================================================='''

# choose next bead based on the probabilities
def chooseBead(currentState: CurrentState,
               terminatingProbability,
               degreeOfBranching,
               treeLevel):
    
    # Check terminating probability
    if random.random() < terminatingProbability:
        # terminate the branch
        currentState.addBead(Bead(TerminalBead().name,
                                  TerminalBead().mass,
                                  TerminalBead().additionalBonds,
                                  TerminalBead().typeOfBead,
                                  TerminalBead().size,
                                  currentState.currentBeadNumber,
                                  treeLevel + 1))
        
    else:
        # select the next bead
        if random.random() < degreeOfBranching:
            currentState.addBead(Bead(DendriticBead().name,
                                      DendriticBead().mass,
                                      DendriticBead().additionalBonds,
                                      DendriticBead().typeOfBead,
                                      DendriticBead().size,
                                      currentState.currentBeadNumber,
                                      treeLevel + 1))
        else:
            currentState.addBead(Bead(LinearBead().name,
                                      LinearBead().mass,
                                      LinearBead().additionalBonds,
                                      LinearBead().typeOfBead,
                                      LinearBead().size,
                                      currentState.currentBeadNumber,
                                      treeLevel + 1))
              
            
            
#----------------------------------------------------------------------------- 



# iterate through the next tree level
def iterateThroughTreeLevel(currentState: CurrentState,
                            terminatingProbability,
                            degreeOfBranching,
                            treeLevel,
                            maximumPolymerMass):
 
    # Iterate through the bead list
    for bead in currentState.beadList:
        
        # Just work on the beads with the correct tree level
        if bead.treeLevel == treeLevel:
        
            # Just add Bonds if there are remaining Bond opportunities left
            while bead.remainingBonds > 0:
                
                # Check if there is enough open mass for new branched beads
                remainingMass = maximumPolymerMass - currentState.currentPolymerMass
                terminatingMass = currentState.currentNumberOfBranches * TerminalBead().mass
                
                if remainingMass <= terminatingMass:
                    return True
                    
            
                chooseBead(currentState, terminatingProbability, degreeOfBranching, bead.treeLevel)
            
                # Add Bond to the bond list
                currentState.bondList.append([bead.identificationNumber,
                                              currentState.currentBeadNumber - 1])
                
                # reduce the remaining bond opportunities
                bead.makeBond()
 
    
    return False
   

#----------------------------------------------------------------------------- 
             
                
def terminateOpenBonds(currentState: CurrentState, terminatingProbability, degreeOfBranching):
    # Iterate through bead list
    for bead in currentState.beadList:
        
        # Find the beads with remaining open bond opportunities
        while bead.remainingBonds > 0:
            currentState.addBead(Bead(TerminalBead().name,
                                      TerminalBead().mass,
                                      TerminalBead().additionalBonds,
                                      TerminalBead().typeOfBead,
                                      TerminalBead().size,
                                      currentState.currentBeadNumber,
                                      treeLevel + 1))
   
            # reduce the remaining bond opportunities
            bead.makeBond()
            
            # Add Bond to the bond list
            currentState.bondList.append([bead.identificationNumber,
                                          currentState.currentBeadNumber - 1])
            
    print("Polymerization successful!")


#-----------------------------------------------------------------------------

# function to determine the bond length based on the bead sizes
def getBondLengths(partner1size, partner2size):
    
    if (partner1size == "tiny" and partner2size == "tiny"):
        return 0.34
    elif (partner1size == "tiny" and partner2size == "small")or (partner1size == "small" and partner2size == "tiny"):
        return 0.365
    elif (partner1size == "tiny" and partner2size == "regular") or (partner1size == "regular" and partner2size == "tiny"):
        return 0.395
    elif (partner1size == "small" and partner2size == "small"):
        return 0.41
    elif (partner1size == "regular" and partner2size == "small") or (partner1size == "small" and partner2size == "regular"):
        return 0.43
    elif (partner1size == "regular" and partner2size == "regular"):
        return 0.47
    else:
        return 1.0


#-----------------------------------------------------------------------------

# function to determine the 2D angle
def getAngle(typeOfBead, bondNumber):    
    
    if typeOfBead == "linear":
        return 0
    elif typeOfBead == "dendritic" and bondNumber == 1:
        return math.pi/3        
    elif typeOfBead == "dendritic" and bondNumber == 2:
        return -math.pi/3         

#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------

# function to calculate the rotated vector based on a rotation matrix
def calculateRotatedVector(xValue, yValue, zValue, angle, length):
    
    # normalize xValue and yValue
    
    normFactor = 1 / math.sqrt((xValue**2) + (yValue**2) + (zValue**2))
    
    xValueNorm = normFactor * xValue
    yValueNorm = normFactor * yValue
    zValueNorm = normFactor * zValue
    
    # introduce new coordinates with 3D rotation matrix in xy plane
    newX = length * (xValueNorm * math.cos(angle) - yValueNorm * math.sin(angle))
    newY = length * (xValueNorm * math.sin(angle) + yValueNorm * math.cos(angle))
    newZ = length * zValueNorm
    
    return addRandomZ(newX, newY, newZ, length)


#-----------------------------------------------------------------------------

# function to apply the small pertubation in z direction to minimize overlapping probabilities
def addRandomZ(x, y, z, length):
     z += random.randint(-1000, 1000) / 8000
     
     normFactor = 1 / math.sqrt((x**2) + (y**2) + (z**2))
     
     newX = normFactor * length * x
     newY = normFactor * length * y
     newZ = normFactor * length * z
     
     return newX, newY, newZ


#-----------------------------------------------------------------------------

# calculate the angle between the different bond in the 3D space
def calculate3DAngle(partner2x, partner2y, partner2z,
                     partner1x, partner1y, partner1z,
                     beadx, beady, beadz):
    
    # define the first vector
    vector1x = partner1x - partner2x
    vector1y = partner1y - partner2y
    vector1z = partner1z - partner2z
    
    # define the second vector
    vector2x = beadx - partner1x
    vector2y = beady - partner1y
    vector2z = beadz - partner1z
    
    
    # calculate the values for the angle formula based on the dot product
    dotProduct = (vector1x * vector2x) + (vector1y * vector2y) + (vector1z * vector2z)
    
    absoluteValue1 = math.sqrt((vector1x**2) + (vector1y**2) + (vector1z**2))
    absoluteValue2 = math.sqrt((vector2x**2) + (vector2y**2) + (vector2z**2))
    
    tempValue = dotProduct / (absoluteValue1 * absoluteValue2)
    
    # to avoid rounding problems
    if tempValue > 1:
        tempValue = 1
    elif tempValue < -1:
        tempValue = -1                          
    
    angle3D = math.acos(tempValue)
       
    return angle3D    
     

#-----------------------------------------------------------------------------

# check for overlapping coordinates
def checkForOverlappingCoordinates(listToCheck):
    
    for subList in listToCheck:
        
        for comparingList in listToCheck:
            
            if ((subList[5] == comparingList[5])
                and (subList[6] == comparingList[6])
                and (subList[7] == comparingList[7])
                and (subList[1] != comparingList[1])):
                
                print('There are exactly overlapping beads')
                return True
                
    print('No overlapping beads detected.')
    return False


#-----------------------------------------------------------------------------

# create modified BondList
def createSortedBondList(bondList, idList):

    # create an empty List
    modifiedBondList = []
    
    # iterate through every id - therefore idList is needed
    for i in range(len(idList)):
        
        # create a new List for every bond
        modifiedBond = [i+1]
        
        for bond in bondList:
            # if the id is in the first place of the normal bond List, take the partner bond and append it to the new List
            if bond[0] == i+1:
                # take the partner bond
                modifiedBond.append(bond[1])
        
        # append this new bond with every bond partner of the id to the modiefied bond List
        modifiedBondList.append(modifiedBond)
        
    
    # create the actual list    
    sortedList = []
    
    # iterate through the new bond list, beginning from the last element
    for j in range(len(modifiedBondList), 0, -1):
        
        # every id is exactly two times in the list, so we need 2 dummy bonds
        bond1 = None
        bond2 = None
        
        # create an other empty List to store data in it
        emptyList = []
        
        # iterate through every bond in the modifiedBondList
        for bond in modifiedBondList:
            
            # if the id is found (what sould happen exactly two times
            if modifiedBondList[j-1][0] in bond:
                
                # fill both dummy bonds with the actual two bonds corresponding to the id
                if bond1 is None:
                    # bond1 is now the same reference as bond
                    # (***)
                    bond1 = bond
                    
                else:
                    bond2 = bond
                    
                    # iterate through the first of the bonds
                    for x in range(len(bond1)):
                        
                        # find the position, where id is in bond1
                        if bond1[x] == modifiedBondList[j-1][0]:
                            
                            # replace the id with the corresponding bond containing the id
                            # bond1 is a reference to the modifiedBondList entry assigned at (***)
                            # bond1[x] is the value -> therefore these entries also change here
                            bond1[x] = bond2
                            
                            # append this bond to the empty list
                            emptyList.append(bond1)
                        
                    # put the entry in the sorted list, so that the empty list can be empty again afterwards
                    sortedList = emptyList
                    
                    # set bond1 back to None and move to the next id
                    bond1 = None
                    break
                    
       
    return sortedList


#-----------------------------------------------------------------------------

def createSmilesString(sortedList, idList):
    # initialize smilesString
    smilesString = "N(C"
    
    
    # create String from the sortedList
    listString = str(sortedList)
    
    # remove first two and last two brackets from the string
    listString = listString[2:-2]
    # remove spaces
    listString = listString.replace(" ", "")
    
    # remove braces
    listString = listString.replace("[", "")
    listString = listString.replace("]", "")
    
    # generate a list of Strings from the single string
    stringList = listString.split(",")
    
    
    # iterate through the listString
    for beadid in stringList:
        
        # initialize beadType variable
        beadType = ""
        
        # first bead already initialized
        if beadid == "1":
            continue
        
        # iterate through the other sorted ids
        for i in range(len(idList)):
            
            if str(idList[i][1]) == beadid:
                beadType = idList[i][2]
                
                # there always has to be a C after a closing bracket
                if smilesString[-1:] == ")":
                    smilesString += "C"
                
                # append the atoms to the string depending on the kind of bead
                if beadType == "terminal":
                    smilesString += "CN)"
                elif beadType == "linear":
                    smilesString += "CNC"
                elif beadType == "dendritic":
                    smilesString += "CN(C"
                else:
                    print("Error in smilesString generation", beadType)
                    
                    
    # remove final bracket            
    smilesString = smilesString[:-1]
    
    return smilesString, stringList


#-----------------------------------------------------------------------------

def createSortedMappingList(stringList, idList):
    
    # count the number of Atoms
    numberOfAtoms = 1
    numberOfN = 1
    numberOfC = 2
    numberOfH = 4
    
    # create a sorted mapping List
    mappingList = []
    
    # iterate through the sorted stringList for the N Atoms
    for beadid in stringList:
        # append the necessary information to the list
        mappingList.append([numberOfAtoms, "N" + str(numberOfN), idList[int(beadid)-1][4]])
        numberOfAtoms += 1
        numberOfN += 1
        
    # the first C Atoms manually
    mappingList.append([numberOfAtoms, "C1", "S1"])
    numberOfAtoms += 1
    
    # create a list for unfinished beads (which means that they have another C Atom)
    unfinishedBeads = ["1"]
        
    # iteration process for C Atoms
    for beadid in stringList:
        # first bead already initilized
        if beadid == "1":
            continue
    
        #get the type of the bead
        typeOfBead = idList[int(beadid)-1][2]
        
        # different behaviour for different bead types
        if typeOfBead == "terminal":
            mappingList.append([numberOfAtoms, "C" + str(numberOfC), idList[int(beadid)-1][4]])
            numberOfAtoms += 1
            numberOfC += 1 
            
            # after a terminal bead, there always has to be another C atoms of a former dendritic bead
            # or the end is reached
            if unfinishedBeads:
                # take the last (dendritic) element out of the unfinished bead list
                beadid = unfinishedBeads.pop()
                mappingList.append([numberOfAtoms, "C" + str(numberOfC), idList[int(beadid)-1][4]])
                numberOfAtoms += 1
                numberOfC += 1
                
            
        elif typeOfBead == "linear":
            mappingList.append([numberOfAtoms, "C" + str(numberOfC), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms, "C" + str(numberOfC + 1), idList[int(beadid)-1][4]])
            numberOfAtoms += 2
            numberOfC += 2
            
        elif typeOfBead == "dendritic":
            mappingList.append([numberOfAtoms + 1, "C" + str(numberOfC), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 1, "C" + str(numberOfC + 1), idList[int(beadid)-1][4]])
            numberOfAtoms += 2
            numberOfC += 2
            # append the dendritic bead at the end of the unfinished bead list
            unfinishedBeads.append(beadid)
            
        else:
            print("Error in creation of sorted mappingList, C Atoms")
            
        
    # the first H Atoms manually
    mappingList.append([numberOfAtoms, "H1", "S1"])
    mappingList.append([numberOfAtoms + 1, "H2", "S1"])
    mappingList.append([numberOfAtoms + 2, "H3", "S1"])
    numberOfAtoms += 3   
    unfinishedBeads = ["1"]
        
    # same procedure with different numbers fo the H atoms
    for beadid in stringList:
        # first bead already initilized
        if beadid == "1":
            continue
    
        #get the type of the bead
        typeOfBead = idList[int(beadid)-1][2]
        
        # different behaviour for different bead types
        if typeOfBead == "terminal":
            mappingList.append([numberOfAtoms, "H" + str(numberOfH), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 1, "H" + str(numberOfH + 1), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 2, "H" + str(numberOfH + 2), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 3, "H" + str(numberOfH + 3), idList[int(beadid)-1][4]])
            numberOfAtoms += 4
            numberOfH += 4 
            
            # after a terminal bead, there always has to be another H atoms of a former dendritic bead
            # or the end is reached
            if unfinishedBeads:
                # take the last (dendritic) element out of the unfinished bead list
                beadid = unfinishedBeads.pop()
                mappingList.append([numberOfAtoms, "H" + str(numberOfH), idList[int(beadid)-1][4]])
                mappingList.append([numberOfAtoms + 1, "H" + str(numberOfH + 1), idList[int(beadid)-1][4]])
                numberOfAtoms += 2
                numberOfH += 2
                
            
        elif typeOfBead == "linear":
            mappingList.append([numberOfAtoms, "H" + str(numberOfH), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 1, "H" + str(numberOfH + 1), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 2, "H" + str(numberOfH + 2), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 3, "H" + str(numberOfH + 3), idList[int(beadid)-1][4]]) 
            mappingList.append([numberOfAtoms + 4, "H" + str(numberOfH + 4), idList[int(beadid)-1][4]])
            numberOfAtoms += 5
            numberOfH += 5
            
        elif typeOfBead == "dendritic":
            mappingList.append([numberOfAtoms, "H" + str(numberOfH), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 1, "H" + str(numberOfH + 1), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 2, "H" + str(numberOfH + 2), idList[int(beadid)-1][4]])
            mappingList.append([numberOfAtoms + 3, "H" + str(numberOfH + 3), idList[int(beadid)-1][4]])
            numberOfAtoms += 4
            numberOfH += 4
            # append the dendritic bead at the end of the unfinished bead list
            unfinishedBeads.append(beadid)
            
        else:
            print("Error in creation of sorted mappingList, H Atoms")
            
        
        
    
        
    return mappingList  



#-----------------------------------------------------------------------------
'''
With these two functions you can tune the values for bond length, equilibrium angle
and forece constants
They only affect the outputted ipt file and nothing inside the other programm
functions and features, therefore they only get called when creating the file
'''
# Write the correct values for each bond in the ipt file
def writeDistanceAndForceConstant(i, j, idList):
    
    # Check for type of beads (Dendritic, linear or terminal)
    bT1 = idList[i-1][2]
    bT2 = idList[j-1][2]
    
    typeList = [bT1, bT2]
    
    countTerminal = typeList.count("terminal")
    countLinear = typeList.count("linear")
    countDendritic = typeList.count("dendritic")
    
    resultString = ""
    
    
    # ps or sp
    if(countTerminal == 1 and countLinear == 1):
        resultString = '0.269'.rjust(8) + '51392'.rjust(8)   
    # pt or tp
    elif(countTerminal == 1 and countDendritic == 1):
        resultString = '0.325'.rjust(8) + '14968'.rjust(8)
    # ss
    elif(countLinear == 2):
        resultString = '0.367'.rjust(8) + '47244'.rjust(8)  
    #st or ts
    elif(countLinear == 1 and countDendritic == 1):
        resultString = '0.379'.rjust(8) + '32739'.rjust(8)   
    # tt
    elif(countDendritic == 2):
        resultString = '0.386'.rjust(8) + '16710'.rjust(8)   
    else:
        resultString = '1.000'.rjust(8) + '1000'.rjust(8)
    
    
    
    return resultString


#-----------------------------------------------------------------------------
# Write the correct values for each bond in the ipt file
def writeEquilibriumAngleAndForceConstant(i, j, k, idList):
        
    # Check for type of beads (Dendritic, linear or terminal)
    bT1 = idList[i-1][2]
    bT2 = idList[j-1][2]
    bT3 = idList[k-1][2]
    
    typeList = [bT1, bT2, bT3]
    
    countTerminal = typeList.count("terminal")
    countLinear = typeList.count("linear")
    countDendritic = typeList.count("dendritic")
    
    resultString = ""
    

    # pss or ssp
    if(countTerminal == 1 and countLinear == 2) and (bT2 == "linear"):
        resultString = '140.479'.rjust(8) + '524'.rjust(8)    
    # pst or tsp
    elif(countTerminal == 1 and countLinear == 1 and countDendritic == 1) and (bT2 == "linear"):
        resultString = '136.080'.rjust(8) + '306'.rjust(8)
    # ptp
    elif(countTerminal == 2 and bT2 == "dendritic"):
        resultString = '108.550'.rjust(8) + '211'.rjust(8)
    # pts or stp
    elif(countTerminal == 1 and countLinear == 1 and countDendritic == 1) and (bT2 == "dendritic"):
        resultString = '117.075'.rjust(8) + '169'.rjust(8)
    # ptt or ttp
    elif(countTerminal == 1 and countDendritic == 2) and (bT2 == "dendritic"):
        resultString = '116.768'.rjust(8) + '83'.rjust(8)
    # sss
    elif(countLinear == 3):
        resultString = '145.026'.rjust(8) + '650'.rjust(8)
    # sst or tss
    elif(countLinear == 2 and countDendritic == 1) and (bT2 == "linear"):
        resultString = '151.233'.rjust(8) + '729'.rjust(8)
    # sts
    elif(countLinear == 2  and bT2 == "dendritic"):
        resultString = '106.405'.rjust(8) + '196'.rjust(8)
    # stt or tts
    elif(countDendritic == 2 and countLinear == 1) and (bT2 == "dendritic"):
        resultString = '108.980'.rjust(8) + '152'.rjust(8)
    # tst
    elif(countDendritic == 2 and bT2 == "linear"):
        resultString = '162.837'.rjust(8) + '654'.rjust(8)
    # ttt
    elif(countDendritic == 3):
        resultString = '121.550'.rjust(8) + '168'.rjust(8)
    else:
        resultString = '100.000'.rjust(8) + '100'.rjust(8)
        
    
    return resultString      


#=============================================================================

''' PART 4: Start of the Main Programm ======================================================================================================
========================================================================================================================================='''


# bool variable to repeat process until a full polymerization without overlapping beads
needToRepeat = True


while needToRepeat == True:
    # Create a current State
    currentState = CurrentState()
    
    # add treeLevel variable
    treeLevel = 1
    
    # add a bool for termination
    isTerminating = False
    
    
    # Add the first bead manually
    currentState.addBead(Bead(LinearBead().name,
                              LinearBead().mass,
                              2,
                              "linear",
                              LinearBead().size,
                              currentState.currentBeadNumber,
                              treeLevel))
    
    
    # continue as long as the user defined polymer mass is reached
    while currentState.currentPolymerMass <= maximumPolymerMass:
          
        # As long as the point of starting to terminate is not reached yet
        if not isTerminating:
            isTerminating = iterateThroughTreeLevel(currentState,
                                                    terminatingProbability,
                                                    degreeOfBranching,
                                                    treeLevel,
                                                    maximumPolymerMass) 
            
            # stop programm if polimerization ended too early
            if currentState.currentNumberOfBranches == 0:
                print("Polymerisation ended early")
                break
            
        # if there is just enough mass left for terminating all the open bonds                    
        else:
            terminateOpenBonds(currentState, terminatingProbability, degreeOfBranching)
            
            # append and calculate additional information out of the polymer data
            currentState.appendBondLengths()
            currentState.appendResidues()
            currentState.createMappingList()
            currentState.calculateCoordinates()

            needToRepeat = checkForOverlappingCoordinates(currentState.idList)
    
    
        # increase the tree level variable to iterate through the next level
        treeLevel += 1
        
        

# for debugging purposes
# currentState.printPolymer()


''' create mapping file and smilesString '''
liste = createSortedBondList(currentState.bondList, currentState.idList)
smilesString, stringList = createSmilesString(liste, currentState.idList)
mappingList = createSortedMappingList(stringList, currentState.idList)




''' For debugging purposes '''
#print(currentState.idList)
#print("-----------------")
#print(currentState.bondList)






''' PART 5: Creation of output files ========================================================================================================
========================================================================================================================================='''



try:
    # itp file
    beadFile = open('file.itp', 'w')
    
    # File headers
    beadFile.write('; itp file \n\n')
    beadFile.write('; Martini 3 for PEI \n')
    beadFile.write('[moleculetype] \n')
    beadFile.write('; molname    nrexcl \n')
    beadFile.write('PEI    1 \n\n')
    beadFile.write('[atoms]\n')
    beadFile.write(';     id    Type   Resnr  Residu    Atom    Cgnr  charge\n')
    
    
    for idPair in currentState.idList:
        beadFile.write(str(idPair[1]).rjust(8)
                       + idPair[0].rjust(8)
                       + '1'.rjust(8))
        
        beadFile.write(idPair[0].rjust(8)
                       + idPair[4].rjust(8)
                       + str(idPair[1]).rjust(8)
                       + '0'.rjust(8)
                       + '\n')
    
    beadFile.write('\n[bonds] \n')
    beadFile.write(';      i       j   Funct  Length Force.c. \n')

    for bondPair in currentState.bondList:
        beadFile.write(str(bondPair[0]).rjust(8)
                       + str(bondPair[1]).rjust(8)
                       + '1'.rjust(8)
                       + writeDistanceAndForceConstant(bondPair[0], bondPair[1], currentState.idList)
                       + '\n')
        
    beadFile.write('\n[angles] \n')
    beadFile.write(';      i       j       k   Funct   Angle Force.c. \n')


    for bondPair in currentState.angleList:
        beadFile.write(str(bondPair[0]).rjust(8)
                       + str(bondPair[1]).rjust(8)
                       + str(bondPair[2]).rjust(8)
                       + '2'.rjust(8)
                       + writeEquilibriumAngleAndForceConstant(bondPair[0], bondPair[1], bondPair[2], currentState.idList)
                       + '\n')

    beadFile.close()
    
    
    
    
    # gro file
    coordinateFile = open('file.gro', 'w')
    
    # headers
    # coordinateFile.write('; gromacs file\n\n')
    coordinateFile.write('Mapped structure for PEI following MARTINI 3\n')
    coordinateFile.write(str(len(currentState.idList)) + '\n')
    # coordinateFile.write('; nrpoly    Residu    id    x    y    z\n')
    
    for subList in currentState.idList:
        coordinateFile.write(str(0).rjust(5)
                             + 'PEI'.rjust(5)
                             + str(subList[0]).rjust(5)
                             + str(subList[1]).rjust(5)
                             + str(format(subList[5], '.3f')).rjust(8)
                             + str(format(subList[6], '.3f')).rjust(8)
                             + str(format(subList[7], '.3f')).rjust(8)
                             + '\n')
    
    
    #coordinateFile.write('\n; box size x, y, z\n')
    coordinateFile.write('0.00000'.rjust(8) + '0.00000'.rjust(8) + '0.00000'.rjust(8) + '\n')
    
    coordinateFile.close()



    # mapping file
    mappingFile = open('file.mapBead', 'w')

    # headers
    mappingFile.write('[ to ]\n')
    mappingFile.write('martini\n')

    mappingFile.write('\n[ martini ]\n')
    for bead in currentState.idList:
        mappingFile.write(bead[4] + ' ')

    mappingFile.write('\n\n[ atoms ]\n')


    for subList in currentState.mappingList:
        mappingFile.write(str(subList[0]).ljust(8)
                             + subList[1].ljust(8)
                             + subList[2].ljust(8)
                             + '\n')

    mappingFile.close()
    
    
    
    # mapping file2
    mappingFile2 = open('file.map', 'w')

    # headers
    mappingFile2.write('[ to ]\n')
    mappingFile2.write('martini\n')

    mappingFile2.write('\n[ martini ]\n')
    for bead in currentState.idList:
        mappingFile2.write(bead[4] + ' ')

    mappingFile2.write('\n\n[ atoms ]\n')


    for subList in mappingList:
        mappingFile2.write(str(subList[0]).ljust(8)
                             + subList[1].ljust(8)
                             + subList[2].ljust(8)
                             + '\n')

    mappingFile2.close()
    
    
    
    # smilesString file
    smilesStringFile = open('file.smiles', 'w')

    smilesStringFile.write('SmilesString\n')
    smilesStringFile.write(smilesString)

    smilesStringFile.close()
    
    
    
    
    ''' Here starts the bonus part for the titration model
    if you don't need this just ignore the fileTitration.itp file
    
    Also tube the corresponding values here and not in the code before
    '''
    
    # create a mapping Dict for titration models
    mappingDict = {}
    
    shift = 0
    
    for i in range(len(currentState.idList)):
            
        mappingDict[i+1] = i+1 + shift
        
        # if (currentState.idList[i][0] == 'TN6d'):
            
        #     shift += 2
        
        shift += 2
    


    # create new idList and bondList
    titrationIdList = []
    titrationBondList = []
    
    
    
    # fill the new idList
    for idPair in currentState.idList:
        
        index = mappingDict[idPair[1]]
        
        
        if (idPair[0] != 'TN6d'):
            titrationIdList.append([index,
                                    'N2_10.2', # idPair[0],
                                    '1',
                                    'PEI',
                                    'N2_10.2', #idPair[0],
                                    index,
                                    '-1',
                                    '24'])
            
            titrationIdList.append([index + 1,
                                    'DB1',
                                    '1',
                                    'PEI',
                                    'DB1',
                                    index + 1,
                                    '0',
                                    '24'])
                                   
            titrationIdList.append([index + 2,
                                    'DB2',
                                    '1',
                                    'PEI',
                                    'DB2',
                                    index + 2,
                                    '1',
                                    '24'])
            
            
        else:   
            titrationIdList.append([index,
                                    'SN6d_10.6F',
                                    '1',
                                    'PEI',
                                    'SN6d_10.6F', #idPair[0],
                                    index,
                                    '-1',
                                    '24'])
             
            titrationIdList.append([index + 1,
                                    'DB3',
                                    '1',
                                    'PEI',
                                    'DB3',
                                    index + 1,
                                    '1',
                                    '24'])
                                   
            titrationIdList.append([index + 2,
                                    'DB4',
                                    '1',
                                    'PEI',
                                    'DB4',
                                    index + 2,
                                    '0',
                                    '24'])
            
          
    
    
    # fill the new bondList
    for i in range(len(currentState.bondList)):
        
        # insert the dummy bead bond before the other bonds because its index is smaller
        # check that there is still a next element in the list to compare
        if ((i+1) < len(currentState.bondList)):
            # check that that the bonds do not get inserted twice (case when it is 0 and not 1 or 2)
            # and that they get inserted for linear molecules (statement after the "or")
            # the "and" part is because the first bead initialized is this special linear bead to start the polymerization
            if ((currentState.bondList[i+1][0] - currentState.bondList[i][0] == 0)
                or (currentState.idList[currentState.bondList[i][0] - 1][2] == 'linear' and currentState.idList[currentState.bondList[i][0] - 1][1] != 1)):
                # DummyBead Bond for SN4a
                titrationBondList.append([mappingDict[currentState.bondList[i][0]],
                                          mappingDict[currentState.bondList[i][0]] + 2,
                                         '1',
                                         '0.00'.rjust(8) + '4000'.rjust(8)
                                         ])
        
        # change indizes of the regular bonds
        titrationBondList.append([mappingDict[currentState.bondList[i][0]],
                                  mappingDict[currentState.bondList[i][1]],
                                 '1',
                                 writeDistanceAndForceConstant(currentState.bondList[i][0], currentState.bondList[i][1], currentState.idList)
                                 ])
        
        # insert bonds with dummybeads for terminating beads because their indices were not in the original bond list
        # check if it is not the last element in the bond list
        if ((i+1) < len(currentState.bondList)):
            
            id_difference = currentState.bondList[i+1][0] - currentState.bondList[i][0]
            shift_index = 3
            
            while (id_difference > 1):
                titrationBondList.append([mappingDict[currentState.bondList[i][0]] + shift_index,
                                          mappingDict[currentState.bondList[i][0]] + shift_index + 1,
                                         '1',
                                         '0.00'.rjust(8) + '4000'.rjust(8)
                                         ])
                
                id_difference -= 1
                shift_index += 3
            
    
    
            
    # append the terminating bead bonds with their corresponding dummybeads
    for k in range(len(currentState.idList)):
        if (currentState.idList[k][1] <= currentState.bondList[len(currentState.bondList)-1][0]):
            continue
        
        titrationBondList.append([mappingDict[k+1],
                                 mappingDict[k+1] + 1,
                                 '1',
                                 '0.00'.rjust(8) + '4000'.rjust(8)
                                 ])

    # for debugging purposes
    # print(titrationBondList)
    
    
    
    # Titration itp file
    beadFile2 = open('fileTitration.itp', 'w')
    
    # File headers
    beadFile2.write('; itp file \n\n')
    beadFile2.write('; Martini 3 for PEI \n')
    beadFile2.write('[moleculetype] \n')
    beadFile2.write('; molname    nrexcl \n')
    beadFile2.write('PEI    1 \n\n')
    beadFile2.write('[atoms]\n')
    beadFile2.write(';     nr        type   resnr  residu        atom    cgnr  charge    mass\n')
    
    
    for bead in titrationIdList:
        beadFile2.write(str(bead[0]).rjust(8)
                       + bead[1].rjust(12)
                       + bead[2].rjust(8)
                       + bead[3].rjust(8)
                       + bead[4].rjust(12)
                       + str(bead[5]).rjust(8)
                       + bead[6].rjust(8)
                       + bead[7].rjust(8)
                       + '\n')
        
    
    beadFile2.write('\n[bonds] \n')
    beadFile2.write(';      i       j   Funct  Length Force.c. \n')

    
    for bond in titrationBondList:
        beadFile2.write(str(bond[0]).rjust(8)
                       + str(bond[1]).rjust(8)
                       + bond[2].rjust(8)
                       + bond[3]
                       + '\n')
        
    beadFile2.write('\n[angles] \n')
    beadFile2.write(';      i       j       k   Funct   Angle Force.c. \n')


    for bondPair in currentState.angleList:
        beadFile2.write(str(mappingDict[bondPair[0]]).rjust(8)
                       + str(mappingDict[bondPair[1]]).rjust(8)
                       + str(mappingDict[bondPair[2]]).rjust(8)
                       + '2'.rjust(8)
                       + writeEquilibriumAngleAndForceConstant(bondPair[0], bondPair[1], bondPair[2], currentState.idList)
                       + '\n')
    
        
    beadFile2.write('\n[exclusions] \n')    
    beadFile2.write('\n[bonds] \n')
    beadFile2.write(';      i       j   Funct  Length Force.c. \n')
    
    
    for bead in titrationIdList:
        if (bead[1] == 'SN6d_10.6F'):
            beadFile2.write(str(bead[0]).rjust(8)
                            + str(bead[0]+2).rjust(8)
                            + '1'.rjust(8)
                            + '0.187'.rjust(8)
                            + '10000'.rjust(8)
                            + '\n')
        elif (bead[1] == 'N2_10.2'):
            beadFile2.write(str(bead[0]).rjust(8)
                            + str(bead[0]+1).rjust(8)
                            + '1'.rjust(8)
                            + '0.200'.rjust(8)
                            + '10000'.rjust(8)
                            + '\n')
    

    beadFile2.close()
    
    
    
    # creation of a pH-dependent itp file
    if (pH == 7.4 or pH == 5.5):
        # itp file
        pHFile = open('filepH.itp', 'w')
        
        # File headers
        pHFile.write('; itp file \n\n')
        pHFile.write('; Martini 3 for PEI \n')
        pHFile.write('[moleculetype] \n')
        pHFile.write('; molname    nrexcl \n')
        pHFile.write('PEI    1 \n\n')
        pHFile.write('[atoms]\n')
        pHFile.write(';     id    Type   Resnr  Residu    Atom    Cgnr  charge\n')
        
        
        for idPair in currentState.idList:
            charge = 0
            name=''
            randomNumber = random.random()
            
            if (pH == 7.4):
                if (idPair[0] == 'TN6d'):
                    if (randomNumber <= 0.66):
                        name = 'TQ5p'
                        charge = 1
                    else:
                        name = 'TN6d'
                        charge = 0
                elif (idPair[0] == 'SN4'):
                    if (randomNumber <= 0.33):
                        name = 'SQ4p'
                        charge = 1
                    else:
                        name = 'SN4'
                        charge = 0
                elif (idPair[0] == 'SN3a'):
                    name = 'SN3a'
                    charge = 0
                    
            elif (pH == 5.5):
                if (idPair[0] == 'TN6d'):
                    name = 'TQ5p'
                    charge = 1
                elif (idPair[0] == 'SN4'):
                    if (randomNumber <= 0.66):
                        name = 'SQ4p'
                        charge = 1
                    else:
                        name = 'SN4'
                        charge = 0
                elif (idPair[0] == 'SN3a'):
                    if (randomNumber <= 0.33):
                        name = 'SQ3p'
                        charge = 1
                    else:
                        name = 'SN3a'
                        charge = 0
            
            
            
            
            pHFile.write(str(idPair[1]).rjust(8)
                           + name.rjust(8)
                           + '1'.rjust(8))
            
            pHFile.write(name.rjust(8)
                           + idPair[4].rjust(8)
                           + str(idPair[1]).rjust(8)
                           + str(charge).rjust(8)
                           + '\n')
        
        pHFile.write('\n[bonds] \n')
        pHFile.write(';      i       j   Funct  Length Force.c. \n')

        for bondPair in currentState.bondList:
            pHFile.write(str(bondPair[0]).rjust(8)
                           + str(bondPair[1]).rjust(8)
                           + '1'.rjust(8)
                           + writeDistanceAndForceConstant(bondPair[0], bondPair[1], currentState.idList)
                           + '\n')
            
        pHFile.write('\n[angles] \n')
        pHFile.write(';      i       j       k   Funct   Angle Force.c. \n')


        for bondPair in currentState.angleList:
            pHFile.write(str(bondPair[0]).rjust(8)
                           + str(bondPair[1]).rjust(8)
                           + str(bondPair[2]).rjust(8)
                           + '2'.rjust(8)
                           + writeEquilibriumAngleAndForceConstant(bondPair[0], bondPair[1], bondPair[2], currentState.idList)
                           + '\n')

        pHFile.close()
    



        print('Thank you for using the PolymerGenerator - We hope you enjoy this tool!')


except:
    print("Strange error that nobody can explain")