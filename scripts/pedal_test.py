#!/usr/bin/env python
import rospy
import rosbag
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import pdb

#### INITIALIZE AND DEFINE THINGS
filename = "2018-12-21-12-09-33.bag"
topic1 = "/vehicle/throttle_report"
topic2 = "/vehicle/filtered_accel"
accelVect = [[0,0]]
throttleVect = [[0,0]]

#### READ BAG FILE AND USE READ_MESSAGES METHOD TO MSGS
#### APPEND THEM INTO RESPECTIVE ARRAYS DECLARED BEFORE
bag = rosbag.Bag(filename)
StartTime = bag.get_start_time()
EndTime = bag.get_end_time()
Duration = EndTime - StartTime
for topic, msg, t in bag.read_messages(topics=[topic1,topic2]):
    if topic == topic1:
        throttleVect.append([t.to_sec()-StartTime,msg.pedal_input])
    if topic == topic2:
        accelVect.append([t.to_sec()-StartTime,msg.data])


#### APPEND END TIME TO BOTH VECTS AND CONVERT
throttleVect.append([Duration,0])
accelVect.append([Duration,0])
throttleVect = np.asarray(throttleVect)
accelVect = np.asarray(accelVect)


#pdb.set_trace() # DEBUGGING
#### CREATE INTER1D OBJECTS AND TIME VECTOR
thF = interp1d(throttleVect[:,0],throttleVect[:,1])
acF = interp1d(accelVect[:,0],accelVect[:,1])
tVect = np.linspace(0,Duration,100)

#### Plot
# normalized
temp1 = thF(tVect)/max(thF(tVect))
temp2 = acF(tVect)/max(acF(tVect))

plt.figure()
plt.plot(tVect,temp1)
plt.plot(tVect,temp2)
plt.plot(tVect,abs(temp2/temp1))
plt.legend(['throttle','accel','th:accel'])
plt.show()

    

