#!/usr/bin/env python
# THIS SCRIPT WILL BE USED TO DETERMINE Q, R NOISE COVARIANCE MATRICES
#  and to parse bag files
import csv
import rosbag
import rosbag_pandas
#import rospy
import math
from geometry_msgs.msg import Point32
from geometry_msgs.msg import PointStamped
from nav_msgs.msg import Odometry
from sensor_msgs.msg import NavSatFix
from custom_msgs.msg import positionEstimate
from custom_msgs.msg import orientationEstimate
import numpy as np
import pdb

import rosbag
#### bag file handle
bag = rosbag.Bag('./Data/2018-12-16-17-23-03.bag')
StartTime = bag.get_start_time()
#### Initialize Arrays
#measurements
rtkx = [] 
rtky = []
gpslat = [] 
gpslon = []
gpsyaw = [] #Estimate from lat/lon
carvel = []
imuyaw = [] #Estimate from yawd
imuyawd = []
# inputs
uswheel = [] #read steering angle
uaccel = [] #read accel


####  read_messages
for topic, msg, t in bag.read_messages(topics=['/gps/rtkfix']):
    rtkx.append([t.to_sec()-StartTime,msg.pose.pose.position.x])
    rtky.append([t.to_sec()-StartTime,msg.pose.pose.position.y])
rtkx = np.asarray(rtkx)
rtky = np.asarray(rtky)

for topic, msg, t in bag.read_messages(topics=['/vehicle/gps/fix']):
    gpslat.append([t.to_sec()-StartTime,msg.latitude])
    gpslon.append([t.to_sec()-StartTime,msg.longitude])
gpslat = np.asarray(gpslat)
gpslon = np.asarray(gpslon)

#estimate gpsyaw from lat/lon

for topic, msg, t in bag.read_messages(topics=['/vehicle/twist']):
    carvel.append([t.to_sec()-StartTime,msg.twist.linear.x])
carvel = np.asarray(carvel)

for topic, msg, t in bag.read_messages(topics=['/vehicle/imu/withcovariance']):
    imuyawd.append([t.to_sec()-StartTime,msg.angular_velocity.z])
imuyawd = np.asarray(imuyawd)

#estimate imuyaw from imuyawd. use initial yaw from utm coords

for topic, msg, t in bag.read_messages(topics=['/vehicle/steering_report']):
    uswheel.append([t.to_sec()-StartTime,msg.steering_wheel_angle])
uswheel = np.asarray(uswheel)

for topic, msg, t in bag.read_messages(topics=['/vehicle/filtered_accel']):
    uaccel.append([t.to_sec()-StartTime,msg.data])
uaccel = np.asarray(uaccel)

bag.close()

lf = 1.3314298821150170049065764033003	
lr = 1.5185105279655697341212317041936
base = [30.6286298219,-96.4821284934]


#Fk = np.array([[0, 0, -V*np.sin(Psi), np.cos(Psi), 0],
#          [0, 0, V*np.cos(Psi), np.sin(Psi), 0],
#          [0, 0, 0, lr/(lf+lr)*np.tan(Delf), 0],
#          [0, 0, 0, 0, 0],
#          [0, 0, 0, 0, 1]])
#Gk = np.array([[0,0],
#           [0,0],
#           [V/(lf+lr)*Delf,0],
#           [0,1],
#           [0,0]])
#Hk = np.array([[1,0,0,0,0],
#          [0,1,0,0,0],
#          [0,0,1,0,0],
#          [0,0,0,1,0],
#          [0,0,1,0,0],
#          [0,0,0,0,1]])

#initial covariance



			
	

