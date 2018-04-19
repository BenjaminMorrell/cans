import numpy as np 

import cans 


# initialise object
print(dir(cans.Object3D))

# Default init
obj = cans.Object3D()


# Init with points 
# knots
nKnots = 11
deg = 3
nCtrl = 11 - 3 - 1
knotU = np.array([0,0,0,0,0.1,0.3,0.7,1,1,1,1])
knotV = knotU 

controlP = np.zeros([3,nCtrl*nCtrl])

index = 0
for i in range(nCtrl):
  for j in range(nCtrl):
    controlP[0,index] = j
    controlP[1,index] = i
    controlP[2,index] = 0.0
    index += 1


# init
obj2 = cans.Object3D(deg,deg,knotU,knotV,controlP,nCtrl,nCtrl)


print(dir(obj2))

# ------------------------------------------------
# Get distance to point
point = np.ones([3,1])
point[0] = 3.5
point[1] = 3.5
point[2] = 1.0

dist = obj2.getDistanceFromPointToSurface(point, 25, 25)

print("Distance to point {} is {}".format(point,dist))

# ------------------------------------------------
# Get distance to batch of points 
nPoints = 25
points = np.ones([3,nPoints])

for i in range(nPoints):
  ratio = 1.0*i/nPoints
  points[0,i] = 2.0 + 2.0*ratio
  points[1,i] = 2.0 + 2.0*ratio
  points[2,i] = -2.0 + 4.0*ratio


# query the distance 
distances = np.zeros([1,nPoints])

distances=obj2.getBatchDistanceFromPointsToSurface(points,25,25)

print("Distances are: {}".format(distances))