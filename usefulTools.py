# usefulTools.py
#
# Matt Churchfield
# National Renewable Energy Laboratory
# 8 May 2017
#
# This is a module that contains various useful tools.

 




import os
import numpy as np




# If the data is stored with each time sample in its own directory, this determines
# how many time directories there are, what they are called, and puts them in numerical 
#order.

def getOutputTimes(dir):

  data = os.listdir(dir)
  outputTimesI = []
  outputTimes = []
  nTimes = len(data)
  ii = 0
  for i in range(nTimes):
     if (data[i][0]).isnumeric():
        outputTimesI.append(data[i])
        ii = ii + 1

  nTimes = len(outputTimesI)

  outputTimesIndex = 0
  outputTimesSort = []
  for i in range(nTimes):
     outputTimesSort.append([i,float(outputTimesI[i])])

  outputTimesSort = sorted(outputTimesSort, key=lambda index: index[1])

  for i in range(nTimes):
     outputTimes.append(outputTimesI[outputTimesSort[i][0]])


  return nTimes, outputTimes





# Test to see if a string is contained in a list of strings.
def isStringInList(str,strList):
    
    inList = False
    
    for i in range(len(strList)):
        if (strList[i] == str):
            inList = True
            
    return inList





# Find the index and value of the point in a list nearest the specified point.
def findNearestIndex(array,value):
    array = np.asarray(array)
    idx = int((np.abs(array - value)).argmin())
    valNearest = array[idx]
    
    return idx, valNearest






# Find the nearest values in a list and their indices to a specified point.  Also return
# the weights if linearly interpolating between the bounding values.
def findBoundingIndices(array,value):
    idx,valNearest = findNearestIndex(array,value)
    
    isIncreasing = True if (array[1] > array[0]) else False

    if isIncreasing:
        if (valNearest <= value):
            idxLo = idx
            idxHi = idx+1
        else:
            idxLo = idx-1
            idxHi = idx
    else:
        if (valNearest <= value):
            idxLo = idx
            idxHi = idx-1
        else:
            idxLo = idx+1
            idxHi = idx

    indices = [idxLo,idxHi]
    
    weights = np.zeros((2,))
    weights[0] = (array[idxHi] - value) / (array[idxHi] - array[idxLo])
    weights[1] = (value - array[idxLo]) / (array[idxHi] - array[idxLo])
    

    return indices, weights, valNearest






# Convert from compass direction to Cartesian direction.
def CompassToCartesian(dir):
    dir += 180.0
    
    if (dir >= 360.0):
        dir -= 360.0

    dir = 90.0 - dir
    
    if (dir < 0.0):
        dir += 360.0

    return dir



def CompassToCartesianList(dir):
    for i in range(len(dir)):
        dir[i] = CompassToCartesian(dir[i])
        
    return dir






# Convert from Cartesian direction to compass direction.
def CartesianToCompass(dir):
    dir = 90.0 - dir
    
    if (dir < 0.0):
        dir += 360.0

    dir += 180.0
    
    if (dir >= 360.0):
        dir -= 360.0

    return dir



def CartesianToCompassList(dir):
    for i in range(len(dir)):
        dir[i] = CartesianToCompass(dir[i])
        
    return dir








