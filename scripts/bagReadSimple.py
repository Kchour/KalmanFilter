#!/usr/bin/env python
# THIS SCRIPT WILL BE USED TO DETERMINE Q, R NOISE COVARIANCE MATRICES
#  and to parse bag files
import csv
import rosbag
import math
import numpy as np
import pdb
import rosbag

#### bag file handle
bag = rosbag.Bag('2018-12-21-12-09-33.bag')
StartTime = bag.get_start_time()

#### Initialize Arrays
rtkx = [] 
rtky = []

####  read_messages
print "Reading messages, please wait"
for topic, msg, t in bag.read_messages(topics=['/gps/rtkfix']):
    rtkx.append([t.to_sec()-StartTime,msg.pose.pose.position.x])
    rtky.append([t.to_sec()-StartTime,msg.pose.pose.position.y])
rtkx = np.asarray(rtkx)
rtky = np.asarray(rtky)
bag.close()

pdb.set_trace()
#### Write to file
print "Writing to file, please wait"
filename = "coords.dat"
fopen = open(filename,"w")
for i in range(len(rtkx)):
    fopen.write(str(rtkx[i][1])+", "+str(rtky[i][1])+"\n")
fopen.close()
print "DONE"


####Alternatively
merge = np.vstack((rtkx[:,1],rtky[:,1]))
np.savetxt("test.dat",merge.T,delimiter=',')






			
	

