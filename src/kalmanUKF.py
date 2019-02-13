#!/usr/bin/env python
# THIS SCRIPT WILL BE USED TO DETERMINE Q, R NOISE COVARIANCE MATRICES
#  AND TO PARSE BAG FILES
# Currently fuses: built-in gps, mti-gps, wheel odom(vel and yawrate), built-in imu(lin acc and ang vel), mti-imu(yaw-angle)
# TODO: USE RAW MTI IMU DATA NOT FILTERED
# NOTE CONVERTING TO UTM ALWAYS RESULTS IN DISTORTION!!!
# NOTE A GAIN OF 1 INDICATES GREATER WEIGHTING OF MEASUREMENTS, WHILE 0 INDICATES MORE WEIGHTING ON MODEL
import rospy
import rosbag
import math
from kalman_filter.msg import States
from std_msgs.msg import Float64
from geometry_msgs.msg import Point32
from geometry_msgs.msg import PointStamped
from geometry_msgs.msg import TwistStamped
from nav_msgs.msg import Odometry
from sensor_msgs.msg import NavSatFix
from sensor_msgs.msg import JointState
from sensor_msgs.msg import Imu
from custom_msgs.msg import positionEstimate
from custom_msgs.msg import orientationEstimate
from scipy.interpolate import interp1d
import scipy as sp
import numpy as np
import pdb
import tf
#import utm
import geodesy.utm as utm
import sympy as sym
import sys
from scipy.signal import savgol_filter

#TODO: Create custom messages for the states

##### IMPORT TROUBLESOME TOPICS FIRST
##IMPORT STEERING_REPORT BECAUSE THE BAG MESSAGE IS OLDER THAN INSTALLEDTO
#filename = "./2018-12-21-12-09-33.bag" 
#cmdSteering = [[0,0]]
#bag = rosbag.Bag(filename)
#startTime = bag.get_start_time()
#endTime = bag.get_end_time()
#duration = endTime - startTime
#for topic, msg, t in bag.read_messages(topics=['/vehicle/steering_report']):
#    cmdSteering.append([t.to_sec()-startTime,msg.steering_wheel_angle_cmd])
##pdb.set_trace()
#cmdSteering.append([duration,0])
#cmdSteering = np.asarray(cmdSteering)
#bag.close()
#
##### INTERPOLATION OBJECT
## TO CALL DO CMDSTEERINGF(SOME X)
## DONT FORGET TO DIVIDE BY STEERING RATIO := 14.8
#cmdSteeringF = interp1d(cmdSteering[:,0],cmdSteering[:,1]/14.8)

global wvsteer, wvaccd, wvspeed, wvyawd, jointstate_cbTime, vel_wheels_prev, cb_flag
jointstate_cbTime = [0,0]
vel_wheels_prev = 0
global vimuAVC, vimuLAC
global mtiOP
global vehUTMx, vehUTMy
global mtiUTMx, mtiUTMy
global vehAcc
global vehgpsvel
global piksiUTMx, piksiUTMy
piksiUTMx = 0
piksiUTMy = 0
mtiUTMx = 0
mtiUTMy = 0

cb_flag = [0,0,0,0,0,0,0]    
init_flag = 0

t = 0
tprev = 0
cmdAccprev = 0
wvsteer_prev=0
rospy.init_node("kalman_ukf_filter_testnode")
#### FUNCTION DEFINTIONS
def timer_cb(event):
    global cb_flag
    for i in range(0,len(cb_flag)):
        if i != 4:                                              ####TESTING PURPOSES REMOVE LATER
            cb_flag[i] -= 1
            cb_flag[i] = clamp(cb_flag[i],6,0)

def kalman(Fk,G,H,Zk,Qk,R,Pk,Xk,Uk):
    #### GATHER NUMERIC LIST FOR SUBSTITUTION (ALLVARS)
    # Prediction
    #pdb.set_trace()
    Xk = Fk.dot(Xk) + G.dot(Uk)
    Pk = Fk.dot(Pk).dot(np.transpose(Fk))+Qk
    # Update
    S = Hnum.dot(Pk).dot(Hnum.T) + R        # Measurement Covariance
    if np.linalg.cond(S) > 1/sys.float_info.epsilon:
            pdb.set_trace()
    #Gk = (Pk.dot(np.transpose(H))).dot(np.linalg.inv(H.dot(Pk).dot(np.transpose(H))+R))
    Gk = (Pk.dot(np.transpose(H))).dot(np.linalg.inv(S))
    Xk = Xk + Gk.dot((Zk - H.dot(Xk)))
    Pk = (np.eye(np.shape(Pk)[0])-Gk.dot(H)).dot(Pk)
    return Xk,Pk
    

def jointstate_cb(msg):
    global wvsteer,wvaccd, wvspeed, wvyawd, jointstate_cbTime, vel_wheels_prev,cb_flag
    #JointState
    R = 0.2413#wheel radius
    W = 1.58242#m track width
    L = 2.84988#m wheel base
    #### OBTAIN MESSAGES
    jointstate_cbTime[0] = jointstate_cbTime[1]
    jointstate_cbTime[1] = msg.header.stamp.secs + msg.header.stamp.nsecs*(10**-9)
    steer_flp = msg.position[4] #steer_fl_phi
    steer_frp = msg.position[5] #steer_fr_phi
    #### VEL_WHEELS CALCULATIONS
    vel_wheels = float(sum([msg.velocity[2],msg.velocity[3]]))/2.0  #average vehicle speed rear tires
    dt = jointstate_cbTime[1] - jointstate_cbTime[0]
    vel_wheels_dot = (vel_wheels - vel_wheels_prev)/dt
    vel_wheels_prev = vel_wheels
    #### GET INNER AND OUTER STEERING ANGLES
    if steer_flp >= steer_frp:
        phi_i = steer_flp
        phi_o = steer_frp
    else:
        phi_i = steer_frp
        phi_o = steer_flp
    #### STEERING ANGLE AND VEH YAW RATE CALCULATIONS
    tanphi = 0.5*(np.tan(phi_i)/(1+np.tan(phi_i)*(W/(2*L))) + np.tan(phi_o)/(1-np.tan(phi_o)*(W/(2*L)))) 
    veh_speed = vel_wheels*R
    #alpha = veh_speed/L*tanphi - thetaD    #ang accel: current - previous
    thetaD = veh_speed/L*tanphi
    #theta += (thetaD/50) + 0.5*alpha*(1/2500)  #globalangle
    #### STORE RESULTS AS GLOBAL VALUES
    offset = 0.
    wvspeed = veh_speed+offset
    wvyawd = thetaD
    wvaccd = vel_wheels_dot*R #xdotdot = R*phidotdot
    wvsteer = np.arctan(tanphi)
    #### SET CALLBACK FLAG TO 1
    cb_flag[0] += 1
    cb_flag[0] = clamp(cb_flag[0],5,0)

def vehimu_cb(msg):
    global vimuAVC, vimuLAC, vimuPD, vimuVD, cb_flag
    #Imu, only has ang vel,cov, lin accel,cov
    #row major x,y,z covariances
    vimuPD = msg.angular_velocity.z #psi dot
    vimuVD = msg.linear_acceleration.x #vdot

    temp = msg.angular_velocity_covariance
    vimuAVC = np.reshape(temp,(3,3))    #angular velocity covariance rpy

    temp = msg.linear_acceleration_covariance
    vimuLAC = np.reshape(temp,(3,3))    #linear acceleration covariance xyz
    cb_flag[1] += 1
    cb_flag[1] = clamp(cb_flag[1],5,0)

def mtiori_cb(msg):
    global mtiOP, cb_flag
    #orientationEstimate
    mtiOP=np.radians(msg.yaw) #Orientation Psi
    cb_flag[2] += 1
    cb_flag[2] = clamp(cb_flag[2],5,0)

def vehgps_cb(msg):
    global vehUTMx, vehUTMy, cb_flag
    #NavSatFix
    #temp = utm.from_latlon(msg.latitude, msg.longitude)    #from utm module
    temp = utm.fromLatLong(msg.latitude, msg.longitude)      #from geodesy.utm
    vehUTMx = temp.easting#temp[0]
    vehUTMy = temp.northing#temp[1]
    cb_flag[3] += 1
    cb_flag[3] = clamp(cb_flag[3],5,0)

def mtigps_cb(msg):
    global mtiUTMx, mtiUTMy, cb_flag
    #positionEstimate
    #temp = utm.from_latlon(msg.latitude, msg.longitude)    #using utm module
    temp = utm.fromLatLong(msg.latitude, msg.longitude)     #using geodesy.utm
    mtiUTMx = temp.easting#temp[0]
    mtiUTMy = temp.northing#temp[1]
    cb_flag[4] += 1
    cb_flag[4] = clamp(cb_flag[4],5,0)

def vehacc_cb(msg):
    global vehAcc, cb_flag
    #float64
    #### THIS IS VEHICLE STATE NOT COMMAND!!!
    vehAcc = msg.data			
    cb_flag[5] += 1
    cb_flag[5] = clamp(cb_flag[5],5,0)

def vehgpsvel_cb(msg):
    global vehgpsVel
    #TwistStamped
    vehgpsVel = np.sqrt(msg.twist.linear.x**2 + msg.twist.linear.y**2)
    cb_flag[6] += 1
    cb_flag[6] = clamp(cb_flag[6],6,0)


#### NOT FUSING THIS ONE BECAUSE IT'S CHEATING
def piksigps_cb(msg):
    global piksiUTMx, piksiUTMy
    temp = utm.fromLatLong(msg.latitude,msg.longitude)
    piksiUTMx = temp.easting
    piksiUTMy = temp.northing

def clamp(a,b,c):
    if a >=b:
        return b
    elif a <= c:
        return c
    else:
        return a

def motion_model_f(x,dt,u,lr,lf):
    #Kinematic Vehicle model
    #x0 = Xpos                  u[0] = steering_cmd
    #x1 = Ypos                  u[1] = accel_cmd
    #x2 = Psi
    #x3 = Psidot
    #x4 = v body frame
    #x5 = vdot body frame

    bet = np.arctan(lr/(lf+lr)*np.tan(u[0]))
    x[5] = u[1]
    x[4] = u[1]*dt+x[4]
    x[3] = x[4]/lr*np.sin(bet)
    x[2] = x[4]/lr*np.sin(bet)*dt+x[2]
    x[1] = x[4]*np.sin(x[2]+bet)*dt+x[1]
    x[0] = x[4]*np.cos(x[2]+bet)*dt+x[0]

    return x

def wrapToPi(x):
    x = np.fmod(x+np.pi,2*np.pi)
    if x<0:
        x+= 2*np.pi
    return x-np.pi


def measure_model_h(x):
    # 10 different quantities
    #vehGps X
    #vehGps Y
    #mtiGps X
    #mtiGps Y
    #vehGps SpeedBody
    #vehWheel SpeedBody
    #vehImu AccelBody
    #mtiImuFiltered Yaw
    #vehWheel YawRate
    #vehImu YawRate
    H=np.array([[1,0,0,0,0,0],
             [0,1,0,0,0,0],
             [1,0,0,0,0,0],
             [0,1,0,0,0,0],
             [0,0,0,0,1,0], #vehgpsVel
             [0,0,0,0,1,0],
             [0,0,0,0,0,1],
             [0,0,1,0,0,0],
             [0,0,0,1,0,0],
             [0,0,0,1,0,0]])
    return H.dot(x)

def gen_weights(n,kappa,alpha,bet):
    lambda_ = alpha**2*(n+kappa)-n
    Wc = np.full(2*n+1, 1./(2.*(n+lambda_)))
    Wm = np.copy(Wc)
    Wc[0] = lambda_ /(n+lambda_)+(1. - alpha**2.+bet)
    Wm[0] = lambda_/(n+lambda_)
    return Wc,Wm,lambda_

####COMPUTE SIGMA POINTS
" Inputs: n=system order, P:= initial state covariance matrix, X:=mean"
def gen_sigma(n,P,x,lambda_):
    sigmas = np.zeros((2*n+1,n))
    U = sp.linalg.cholesky((n+lambda_)*P)
    sigmas[0] = x
    for k in range (n):
        sigmas[k+1] = x + U[k]
        sigmas[n+k+1] = x - U[k]
    return sigmas




#### ODOMETRY
rospy.Subscriber("test/vehicle/joint_states",JointState,jointstate_cb)
rospy.Subscriber("test/vehicle/imu/withcovariance",Imu, vehimu_cb)
rospy.Subscriber("test/mti/filter/orientation",orientationEstimate, mtiori_cb)
#### GPS
rospy.Subscriber("test/vehicle/gps/fix/withcovariance",NavSatFix, vehgps_cb)
rospy.Subscriber("test/mti/filter/position", positionEstimate, mtigps_cb)      #TESTING PURPOSES
rospy.Subscriber("test/gps/fix",NavSatFix,piksigps_cb)
#### GPS VEL
rospy.Subscriber("test/vehicle/gps/vel",TwistStamped,vehgpsvel_cb)
#### ADDITIONAL FILTERED FEEDBACK INFO
rospy.Subscriber("test/vehicle/filtered_accel",Float64,vehacc_cb)
#### WATCHDOG TIMER
rospy.Timer(rospy.Duration(2), timer_cb)

#### PUBLISHER
pubinfo = rospy.Publisher("test/kalman/states",States,queue_size=1000)
pubtest = rospy.Publisher("test/kalman/test",States,queue_size=1000)
SI = States()
SI.name = ['X','Y','Psi','Psidot','V','Vdot']
numStates = 6
SI.states = numStates*[None]

#### TESTING PURPOSES
Stest = States()
Stest.states=21*[0]

lf = 1.3314298821150170049065764033003	
lr = 1.5185105279655697341212317041936

# Let's run the prediction with a larger dt, while maintaining a different loop rate
des_rate = 100.0
#rate = rospy.Rate(des_rate)
rate = rospy.Rate(des_rate)

count = 0           # ADAPTIVE KALMAN FILTER COUNTER
eps = 0             # ADAPTIVE KALMAN FILTER EPS
kfwin = 100         # SAVITZKY GOLAY, KF WIND
svgwin = 25         # SVG WIN
svgorder = 3
kfvect = []         # kf vector of size win
svresults=np.array([kfwin*[0]])

n = 6           #dimension of state space
nz = 10         #dimension of measurement space
alp = 0.6
bet = 2.0
kap = 3.0 #6 is the order...change 3 to 7?
Xk = np.zeros(n)     #mean of states?
Px = np.full((n,n),0.0)      #covariance associated with x, the mean states
Pz = np.full((nz,nz),0.0)    #covariance associated with z, the mean measurements
Pxz = np.full((n,nz),0.0)    #cross covariance associated with x and z
Z = np.full((2*n+1,10),0.0 )     #result of passing Sigmas through function h
Y = np.full((2*n+1,n),0.0)  #result of passing Sigmas through function F

timeInit = rospy.get_time()
print "READY TO BEGIN, START ROSBAG"
while not rospy.is_shutdown():
    if (sum(cb_flag )>= 5) and init_flag==1:

        #### CALCULATE UK VECTOR
        Ukdot = np.array([[(wvsteer-wvsteer_prev)*des_rate], [(wvaccd-cmdAccprev)*des_rate]])
        cmdAccprev = wvaccd
        wvsteer_prev = wvsteer
        #Uk = np.array([[wvsteer],[wvaccd]])   #steering tire angle, desired veh acceleration
        Uk = np.array([[wvsteer],[vimuPD]])   #steering tire angle, desired veh acceleration
        "PREDICTION PHASE"
        #### GET WEIGHTS AND SIGMAS
        #x is like the states' mean. Input params: n=6, kap:=3-n, alp:=[0...1], bet:=2.
        #gen_weights returns weights Wc, Wm and parameter lambda_. gen_sigma returns a matrix
        #of sigma points. Input params: n=6, Pk:=state cov matrix"
        Wc,Wm,lambda_=gen_weights(n,kap,alp,bet)
        sigmas = np.copy(gen_sigma(n,Pk,Xk,lambda_))
        #### PASS SIGMAS THRU NONLINEAR FUNCTION F
        for k in range(2*n+1):
            Y[k,:] = np.copy(motion_model_f(sigmas[k,:],1./des_rate,Uk,lr,lf))
        #####GET MEAN x AND COVARIANCE USING Y
        x = np.dot(Wm,Y)   #Scale each row of Y by an element in Wm
        #scale each matrix
        for k in range(2*n+1):
            Px += Wc[k] * np.outer(Y[k]-x,Y[k]-x)
        Px+= Qk
        " UPDATE PHASE"
        #### Pass new sigma point, Y, through function h, which converts them to points in mesurement space. I        have 10 measuredments
        for k in range(2*n+1):
            Z[k] = measure_model_h(Y[k])

        ### GET MEAN z AND COVARIANCE USING Z
        z = np.dot(Wm,Z)    #Scales each row of Z by an element in Wm
        for k in range(2*n+1):
            Pz += Wc[k] * np.outer(Z[k]-z,Z[k]-z)   #Scales each matrix by an element in Wc
        Pz += R 
        
        #### GET MEASUREMENT VALUES FROM SENSORS
        Zm = np.array([vehUTMx-initPosX+X0[0],
                       vehUTMy-initPosY+X0[1],
                       mtiUTMx-initPosX+X0[0],
                       mtiUTMy-initPosY+X0[1],
                       vehgpsVel,
                       wvspeed,
                       vimuVD,
                       mtiOP,
                       wvyawd,
                       vimuPD])
        
        #### compute innovation or residue
        y = np.copy(Zm - z)
        
        #### Compute cross covariance of state and measurements
        for k in range(2*n+1):
            Pxz+= Wc[k]*np.outer(Y[k]-x,Z[k]-z)
        
        ####Kalman Gain
        K = Pxz.dot(np.linalg.inv(Pz))
        
        #### UPDATE
        Xk = x + K.dot(y)
        Xk[2] = wrapToPi(Xk[2])
        Pk = Px - K.dot(Pz).dot(K.T)

        #print "PK: ", Pk,"Pk norm: ",np.linalg.norm(Pk)
        #### Try Savitzky golay filter
        #### create a flatten list out of lists of lists
        temp = [item for item in Xk]
        if len(kfvect) < kfwin:
            kfvect.append(temp)  
        else:
            kfvect.pop(0)
            kfvect.append(temp)
        if (len(kfvect) > svgorder+2):
            len_ = len(kfvect)
            if np.fmod(len_,2)==0:
                len_ -= 1
            temp = np.asarray(kfvect)
            svresults = savgol_filter(temp,len_,svgorder,axis=0)
        for i in range(0,numStates):
            SI.states[i] = Xk[i]    
        pubinfo.publish(SI)



        Stest.states[0] = Zm[0]  #vehGps
        Stest.states[1] = Zm[1]  #vehGps
        Stest.states[2] = Zm[2]  #MtiGps
        Stest.states[3] = Zm[3]  #MtiGps
        Stest.states[4] = piksiUTMx-initPosX  #Piksi
        Stest.states[5] = piksiUTMy-initPosY  #Piksi
        #pdb.set_trace()
        len1_ = np.ceil(len(svresults)/3.0)
        len2_ = np.ceil(len(svresults)/10.0)
        len3_ = np.ceil(len(svresults)/5.0)
        Stest.states[6] = svresults[int(-len1_)][0]  #svgol filter posX
        Stest.states[7] = svresults[int(-len1_)][1]  #svgol filter posY
        Stest.states[8] = svresults[-1][3]          #svgol filter yawrate
        Stest.states[9] = svresults[int(-len3_)][4]          #svgol filter velocity
        Stest.states[10] = np.degrees(Xk[2])         # estimated yaw state
        pubtest.publish(Stest)
    else:
        if init_flag == 0 and all(i >= 1 for i in cb_flag):
            #Set initial state
            initPosX = np.mean([vehUTMx,mtiUTMx])
            initPosY = np.mean([vehUTMy,mtiUTMy])
            #X0 = np.array([[7.3794], [-11.45], [mtiOP], [np.mean([wvyawd,vimuPD])], [np.mean([vehgpsVel,wvspeed])],[vimuVD]])
            X0 = np.array([0, 0, mtiOP, np.mean([wvyawd,vimuPD]), np.mean([vehgpsVel,wvspeed]),vimuVD])
            Xk = np.copy(X0)
            for i in range(0,numStates):
                SI.states[i] = Xk[i] 
            cmdAccprev = wvaccd
            ####initial covariance matrices
            #initial covariance
            #Pk = 0.5*np.identity(6) #Initial State Covariance
            Pk = np.full((n,n),1)+50.0*np.identity(6)
            #Qk = 0.05*np.identity(6) #Process Noise covariance..we can scale this later by residual 500
            dt = 1./des_rate
            Qk = np.diag([0.5*10.0*dt**2,0.5*10.0*dt**2,0.5*1.0*dt**2,1.0*dt,10.0*dt,10.0])*1e-12
            # xgps, ygps, xmgps, ymgps, vehgpsvel, vwo, vimudot, psimti, psiwodot, psiimudot 
            #R = 2.0*np.diag([5.0,5.0,5.0,5.0,1.0,3.0,2.0,0.8,0.5,0.5])    #Measurement Noise Covariance 5
            #R = 10*np.identity(10)
            # Calculate acceleration from vwo?
            R = np.diag([5.0,5.0,2.0,2.0,5.0,25.00,0.05,0.25,1.00,0.5])    #Measurement Noise Covariance
            #R = np.diag([1,1,0.05,0.05,1,0.045,0.20,0.002,0.045,0.01])    #Measurement Noise Covariance
            #TODO: if cb_flag is zero then set that measurement to something high!
            init_flag = 1
    rate.sleep()

	

