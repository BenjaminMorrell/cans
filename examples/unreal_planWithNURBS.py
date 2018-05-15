import rospy
import sys, os
from python_qt_binding.QtWidgets import QApplication
# from minsnap import *
import pickle
import tf

import cans
from cans_msgs.msg import Object3D
import numpy as np
import scipy.io

from geometry_msgs.msg import PoseStamped
from std_msgs.msg import Int32

import torq_gcs
from torq_gcs.plan.astro_plan import QRPolyTrajGUI
# import diffeo



# import struct
class Planner:

  def __init__(self):
    

    # Global goal variable
    self.goal = dict()
    # Global start variable
    self.start = dict()

    self.app = QApplication( sys.argv )

    self.global_dict = dict()

    # Time for publishing
    self.setpointCount = 0
    self.time = 0.0
    self.elapsedTime = 0.0
    self.tmax = 100.0
    self.averageVel = 0.01 # m/s

    # Times for replanning
    self.replanHz = 0.1
    self.timeOfReplan = 0.0
    self.startDeltaT = 0.5 # Time ahead of current time to use as the start location
    self.firstPlan = True

    self.blockUpdate = False

    # Initialise the map
    self.nurbs = []#cans.Object3D()

    # Init planner GCS object
    self.planner = torq_gcs.plan.astro_plan.QRPolyTrajGUI(self.global_dict,defer=True,curv_func=False)

    self.global_dict['fsp_out_map'] = None

    self.planner.inflate_buffer = 0.0
    self.planner.quad_buffer = 0.1
    # self.planner.nurbs_weight = 1e-27 # UNREAL
    self.planner.nurbs_weight = 1e-22

  def initialisePlanner(self):
    # Waypoints
    waypoints = dict()
    for key in self.start.keys():
      waypoints[key] = np.zeros([1,2])
      waypoints[key][0,0] = self.start[key][0]
      waypoints[key][0,1] = self.goal[key][0]

    # Set time to complete the trajectory
    self.computeTrajTime()
    self.planner.seed_times = np.array([self.tmax])
    
    self.global_dict['disc_out_waypoints'] = waypoints

    # Initial guess to load waypoints and initialise planner
    self.planner.on_disc_updated_signal()

    # Set mutation strength
    # self.planner.qr_polytraj.mutation_strength = 1e-26 # UNREAL
    self.planner.qr_polytraj.mutation_strength = 1e-23

  def computeTrajTime(self):
    if not'x' in self.start.keys():
      print("Need start and end first")
      return 
    
    # Compute the total time based on the Euclidean distance from start to end and the desired average velocity
    dist =  0.0
    for key in ['x','y','z']:
      dist += (self.goal[key][0] - self.start[key][0])**2

    dist = np.sqrt(dist)

    self.tmax = dist/self.averageVel

    print("\n\t\tDistance is {}, Traj time is {} s".format(dist,self.tmax))

  def updateStart(self,position):
    # position is a dict with keys x, y, z
    self.planner.on_waypoint_control_callback(position, 0, "main")

  def updateGoal(self,position):
    # position is a dict with keys x, y, z
    self.planner.on_waypoint_control_callback(position, np.shape(self.planner.qr_polytraj.waypoints['x'])[1]-1,"main")

  def updateWaypoints(self,start, goal):

    self.computeTrajTime()
    self.updateStart(start)
    self.updateGoal(goal)

  def updateNURBSObstacle(self):
    
    # Remove nurbs constraints if there are any already
    self.planner.qr_polytraj.remove_nurbs_constraint()

    # Add NURBS constraint 
    # Custom weight
    # self.planner.load_nurbs_obstacle(self.nurbs)
    
    # Fixed weight
    for i in range(0,len(self.nurbs)):
      self.planner.load_nurbs_obstacle(self.nurbs[i],sum_func=True,custom_weighting=True)
      # self.planner.nurbs_weight = 1e-5#12
      # self.planner.on_nurbs_weight_update_button_clicked()
      # TODO(BM) find a more efficient way to update the NURBS for planning 

  def resetGoalinClass(self):
    # Take goal from goal in the planner - if moved in the gui
    for key in self.goal.keys():
      self.goal[key] = self.planner.qr_polytraj.waypoints[key][:,-1] # Last waypoint

  def planTrajectory(self):
    # Runs the optimisation and updates the trajectory 
    # self.planner.on_run_astro_button_click()

    self.planner.qr_polytraj.exit_on_feasible = True
    self.planner.qr_polytraj.optimise(mutate_serial=5)
    self.planner.qr_polytraj.get_trajectory()
    self.planner.update_path_markers()
    self.saveTrajectory()
    
  def saveTrajectory(self,filename="test_traj.traj"):

    qrp_out = self.planner.qr_polytraj

    with open(filename, 'wb') as f:
      # Remove ESDF constraint if there are any already
      qrp_out.remove_esdf_constraint()

      pickle.dump(qrp_out, f, 2 )

    scipy.io.savemat('/home/bjm/Dropbox/PhD_Code/Results/Traj/traj_opt.mat', qrp_out.state_combined)

  def readNURBSMessage(self,msg):
    
    if self.blockUpdate:
      return
    # Print parameters
    print("ID: {}\nDegU: {}, DegV: {}\nnCtrl: ({}, {})".format(msg.ID,msg.degU,msg.degV,msg.nCtrlS,msg.nCtrlT))

    # Get data out
    controlPoints = np.zeros([3,msg.nCtrlS*msg.nCtrlT])
    controlPoints[0,:] = msg.ControlPoints_x
    controlPoints[1,:] = msg.ControlPoints_y
    controlPoints[2,:] = msg.ControlPoints_z

    knotU = np.array(msg.knotU)
    knotV = np.array(msg.knotV)

    print("Knots are:\nU\n{}\nV\n{}\nControl points:\n{}\n\n".format(knotU,knotV,controlPoints))
    
    # Update the object
    if msg.ID >= len(self.nurbs):
      print("Adding objects to map with received ID: {}, and existing map size: {}".format(msg.ID,len(self.nurbs)))
      for i in range(0,msg.ID - len(self.nurbs)+1):
        self.nurbs.append(cans.Object3D())
    
    self.nurbs[msg.ID].updateObject3D(msg.degU,msg.degV,knotU,knotV,controlPoints,msg.nCtrlS,msg.nCtrlT)

    # import pdb; pdb.set_trace()

    self.updateNURBSObstacle()

    # # Run planner
    # if self.time > 1/self.replanHz or self.firstPlan: # If the time since the last replan is more than the desired period
    #   self.resetStartFromTraj()
    #   print("\n\nTime to replan ({}): Running ASTRO\n\n".format(self.time))
    #   self.time = 0.0 # Starting at the start of the new trajectory
    #   self.updateNURBSObstacle()
    #   self.planTrajectory()
    #   print("\n\n\t\t COMPLETED TRAJECTORY PLAN \n\n")
      
    #   # Reset times:
    #   # self.timeOfReplan = self.time
      

    #   self.firstPlan = False
    #   self.saveTrajectory()
      

  def getSetpointAtTime(self):

    # Get the state at the current time being tracked
    output = self.planner.on_animate_eval_callback(self.time)
    # output format is: (t_max, x, y, z, q[0], q[1], q[2], q[3])

    # Set tmax to match 
    self.tmax = output[0]

    if self.time > self.tmax:
      self.time = self.tmax
      print("Time set to tmax = {}".format(self.tmax))

    # Fill message
    msg = PoseStamped()

    # Header
    msg.header.seq = self.setpointCount
    self.setpointCount = self.setpointCount+1 # increment count
    msg.header.stamp = rospy.get_rostime()
    msg.header.frame_id = "local_origin"

    # Position
    msg.pose.position.x = output[1]
    msg.pose.position.y = output[2]
    msg.pose.position.z = output[3]

    # Orientation
    msg.pose.orientation.w = output[4]
    msg.pose.orientation.x = output[5]
    msg.pose.orientation.y = output[6]
    msg.pose.orientation.z = output[7]
    
    # Publish message
    # print("Computed setpoint is:\n{}".format(msg))

    return msg

  

  def resetStartFromTraj(self):
    # Resets the start location and the planned time

    # Time to start replan
    # startTime = self.time + self.startDeltaT
    startTime = self.elapsedTime + self.startDeltaT
    
    print("Getting start at time: {}".format(startTime))
    # State for the start
    self.start = self.planner.get_state_at_time(startTime)
    
    self.updateStart(self.start)

    # Get goal from GUI modifications
    self.resetGoalinClass()

    # Reset traj time
    self.computeTrajTime() # updates self.tmax
    
    # Reset planned trajectory time
    self.planner.qr_polytraj.update_times([0],self.tmax,defer=True)
    # self.planner.qr_polytraj.update_times([0],self.tmax-startTime,defer=True)

    
    # self.tmax -= startTime
    

    print("\n\nRESET: New duration is {}\nStart Location is: {}".format(self.tmax,self.start))

  def goalCallback(self, msg):

    if self.blockUpdate:
      return
    # Reads a call message from Unreal and resets the goal, then replans
    self.goal['x'][0] = msg.pose.position.x
    self.goal['y'][0] = msg.pose.position.y
    self.goal['z'][0] = msg.pose.position.z

    print("\n\n\t\t READ new goal \n\n".format(self.goal))

    # Update goal
    self.updateGoal(self.goal)

    # Reset traj time
    self.computeTrajTime() # updates self.tmax
    
    # Reset planned trajectory time
    self.planner.qr_polytraj.update_times([0],self.tmax,defer=True)
    # self.planner.qr_polytraj.update_times([0],self.tmax-startTime,defer=True)

    # # Run planner
    # if True: #self.time > 1/self.replanHz or self.firstPlan: # If the time since the last replan is more than the desired period
    #   self.setupAndRunTrajectory()
      # self.resetStartFromTraj()
      # print("\n\nTime to replan ({}): Running ASTRO\n\n".format(self.time))
      # self.time = 0.0 # Starting at the start of the new trajectory
      # self.updateEsdfObstacle()
      # self.planTrajectory()
      # print("\n\n\t\t COMPLETED TRAJECTORY PLAN \n\n")


  def updateStartFromTF(self,trans,rot):
    # Callback to update start from tf
    
    # Only updates the position (not the acceleration)
    self.start['x'][0] = trans[0]
    self.start['y'][0] = trans[1]
    self.start['z'][0] = trans[2]

    self.updateStart(self.start)

  def updateElapsedTime(self, msg):

    self.elapsedTime = float(msg.data)/10.0**6

    print("Elapsed time updated to {}".format(self.elapsedTime))


  def setupAndRunTrajectory(self):

    self.blockUpdate = True
    
    self.resetStartFromTraj()
    print("\n\nTime to replan ({}): Running ASTRO\n\n".format(self.time))
    self.time = 0.0 # Starting at the start of the new trajectory

    if not self.firstPlan:
      self.updateNURBSObstacle()
    
    self.planTrajectory()
    print("\n\n\t\t COMPLETED TRAJECTORY PLAN \n\n")
    
    print("\n\n\t\t SENDING TRAJECTORY... \n\n")
    self.planner.on_send_trajectory_button_click()
    

    self.firstPlan = False
    # self.saveTrajectory()

    self.blockUpdate = False

if __name__ == '__main__':

  # Start node
  rospy.init_node('astro',anonymous=True)

  # Init class
  plan = Planner()

  # Set up start and goal 
  # For simple test case
  plan.start['x'] = [0.8]
  plan.start['y'] = [-2.0]
  plan.start['z'] = [0.6]
  plan.start['yaw'] = [0.0]
  plan.goal['x'] = [0.8]
  plan.goal['y'] = [1.0]
  plan.goal['z'] = [0.6]
  plan.goal['yaw'] = [0.0]

  # For 67P test case
  # plan.start['x'] = [-12.0]
  # plan.start['y'] = [0.0]
  # plan.start['z'] = [0.0]
  # plan.start['yaw'] = [0.0]
  # plan.goal['x'] = [11.5]
  # plan.goal['y'] = [-2.0]
  # plan.goal['z'] = [0.0]
  # plan.goal['yaw'] = [0.0]

  plan.initialisePlanner()

  # plan.delayTime = 10.0

  # Example use case:
  plan.updateWaypoints(plan.start,plan.goal)

  # Create Subscriber for object3D objects
  rospy.Subscriber("/object",Object3D,plan.readNURBSMessage)

  # Create Subscriber for goal
  rospy.Subscriber("/goal_unreal",PoseStamped,plan.goalCallback)

  # Create Subscriber for elapsed time 
  rospy.Subscriber("/elapsed_time_unreal",Int32,plan.updateElapsedTime)

  # pub = rospy.Publisher("topic",String,queue_size=1)
  # setpoint_pub = rospy.Publisher("setpoint",PoseStamped,queue_size=1)

  # Flag to subscribe to the tf from unreal to update the start position for planning
  bUseTFToUpdateStart = False

  if bUseTFToUpdateStart:
    # TF listener
    tf_listener = tf.TransformListener()

  rateHz = 5.0
  timer = 0.0

  r = rospy.Rate(rateHz) # 10hz
  while not rospy.is_shutdown():
      # pub.publish("hello")
      # msg = plan.getSetpointAtTime()
      # setpoint_pub.publish(msg)
      # increment time
      # plan.time += 1.0/rateHz # TODO WARNING - this is not going to accurately track time
      
      if timer > 1.0/ plan.replanHz:
        # plan.firstPlan = False
        timer = 0.0
        plan.setupAndRunTrajectory()
        

      plan.time += 1.0/rateHz
      timer += 1.0/rateHz 


      if bUseTFToUpdateStart:
        try:
          (trans, rot) = tf_listener.lookupTransform('/world', '/cam_optical', rospy.Time(0)) # time argument just gets the latest
        except (tf.LookupException, tf.ConnectivityException, tf.ExtrapolationException):
          pass
    
      #   plan.updateStartFromTF(trans,rot)

      
        
      r.sleep()
      
      


  

  
