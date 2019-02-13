#!/usr/bin/env python
import rospy
import math
from geometry_msgs.msg import Point32
from geometry_msgs.msg import PointStamped
from nav_msgs.msg import Odometry
from sensor_msgs.msg import NavSatFix
from custom_msgs.msg import positionEstimate
from custom_msgs.msg import orientationEstimate
#import numpy as np

class CollectWaypoints():

	def __init__(self):
                #Need to parse data csv files
		
                lf = 1.3314298821150170049065764033003	
                lr = 1.5185105279655697341212317041936

                #### Subscribe to relevant topics
                # gps
                rospy.Subscriber("/vehicle/perfect_gps/fix",NavSatFix,self.gps_cb) #convert to utm
                # speed
                rospy.Subscriber("/vehicle/twist",TwistStamped,self.twist_cb)
                # orientation
                rospy.Subscriber("/vehicle/perfect_gps/utm",Vector3Stamped,self.ori_cb)
                rospy.Subscriber("/vehicle/imu/data_raw",Imu,self.imu_cb)
                #cmd inputs
                rospy.Subscriber("/vehicle/steering_cmd",SteeringCmd,self.steer_cb) #steering wheel angle
                rospy.Subscriber("/vehicle/req_accel",Float64,self.accel_cb)
                rate = ropsy.Rate(50)
                while not rospy.is_shutdown():

                    Fk = np.array([[0, 0, -V*np.sin(Psi), np.cos(Psi), 0],
                                  [0, 0, V*np.cos(Psi), np.sin(Psi), 0],
                                  [0, 0, 0, lr/(lf+lr)*np.tan(Delf), 0],
                                  [0, 0, 0, 0, 0],
                                  [0, 0, 0, 0, 1]])
                    Gk = np.array([[0,0],
                                   [0,0],
                                   [V/(lf+lr)*Delf,0],
                                   [0,1],
                                   [0,0])
                    Hk = np.array([1,0,0,0,0],
                                  [0,1,0,0,0],
                                  [0,0,1,0,0],
                                  [0,0,0,1,0],
                                  [0,0,1,0,0],
                                  [0,0,0,0,1]])


			
	def gpsCallback(self,data):

		position = Point32()
		position.x = data.latitude
		position.y = data.longitude
		print position, "\n"

		if len(self.points)<1:
			self.points.append([position.x, position.y])
			#negative indices access arrays from the right instead of the left
		else:
			self.points.append([position.x, position.y])
			#print self.points[len(self.points)-1][0]
			#print self.points
			#write to a file line by line
			f = open(self.filename, 'w')
			for state in self.points:
				f.write(str(state[0]) + ','+str(state[1]) + '\n')
			f.close()

if __name__=='__main__':
	rospy.init_node('collect_waypoints')
	collect_waypoints = CollectWaypoints()
	rospy.spin()
	

