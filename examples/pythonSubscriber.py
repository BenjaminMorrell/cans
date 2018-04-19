import rospy
import sys, os

import cans

from cans_msgs.msg import Object3D
import numpy as np

obj2 = cans.Object3D()

def readObject3DMessage(msg):
  
  # Print parameters
  print("ID: {}\nDegU: {}, DegV: {}\nnCtrl: ({}, {})".format(msg.ID,msg.degU,msg.degV,msg.nCtrlS,msg.nCtrlT))

  # Get data out
  knotU = msg.knotU
  knotV = msg.knotV

  controlPoints = np.zeros([3,msg.nCtrlS*msg.nCtrlT])
  controlPoints[0,:] = msg.ControlPoints_x
  controlPoints[1,:] = msg.ControlPoints_y
  controlPoints[2,:] = msg.ControlPoints_z

  print("Knots are:\nU\n{}\nV\n{}\nControl points:\n{}\n\n".format(knotU,knotV,controlPoints))

  # Create NURBS Object3D

  nKnots = np.size(knotU)
  deg = msg.degU
  nCtrl = msg.nCtrlS
  
  # init
  obj2.updateObject3D(deg,deg,knotU,knotV,controlPoints,nCtrl,nCtrl)

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

if __name__ == '__main__':

  # Start node
  rospy.init_node('object3DListener',anonymous=True)

  # Create Subscriber
  rospy.Subscriber("/object",Object3D,readObject3DMessage)

  rospy.spin()

  # r = rospy.Rate(rateHz) # 10hz
  # while not rospy.is_shutdown():
  #     # pub.publish("hello")
  #     msg = plan.getSetpointAtTime()
  #     setpoint_pub.publish(msg)
  #     # increment time
  #     plan.time += 1.0/rateHz # TODO WARNING - this is not going to accurately track time
  #     r.sleep()
      
      


  

  
