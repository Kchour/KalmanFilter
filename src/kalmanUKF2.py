#!/usr/bin/env python
# THIS SCRIPT WILL BE USED TO DETERMINE Q, R NOISE COVARIANCE MATRICES
#  AND TO PARSE BAG FILES
# Currently fuses: built-in gps, mti-gps, wheel odom(vel and yawrate), built-in imu(lin acc and ang vel), mti-imu(yaw-angle)
# TODO: USE RAW MTI IMU DATA NOT FILTERED
# NOTE CONVERTING TO UTM ALWAYS RESULTS IN DISTORTION!!!
# NOTE A GAIN OF 1 INDICATES GREATER WEIGHTING OF MEASUREMENTS, WHILE 0 INDICATES MORE WEIGHTING ON MODEL

# FIX MODEL YAW TO BE BETWEEN -180 and 180!
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
jointstate_cbTime = [0.0,0.0]
vel_wheels_prev = 0.0
vel_wheels_dot_prev = 0.0
global vimuAVC, vimuLAC
global mtiOP
global vehUTMx, vehUTMy
global mtiUTMx, mtiUTMy
global vehAcc
global vehgpsvel
global piksiUTMx, piksiUTMy
piksiUTMx = 0.0
piksiUTMy = 0.0
mtiUTMx = 0.0
mtiUTMy = 0.0

#1=VGPSx, 2=VGPSy, 3=VGPSVEL, 4=VWheelVEL, 5=VIMUAccel, 6=MTIYaw, 7=VWheelYawD, 8=VImuYAWD,9=MTIGPSx, 10=MTIGPSy
# Test without using MTI gps for now
cb_flag = [0,0,0,0,0,0,0,0,0,0]   #Set to 1 and disable callback (Lines 277-287) if you wish to remove that sensor 
                            # AND TIMER TOO

cb_flag_prev = list(cb_flag)
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
        cb_flag[i] -= 1
        cb_flag[i] = clamp_reset(cb_flag[i],1000,1)

def kalman_predict(n,kap,alp,bet,Pk,Xk,des_rate,Uk,lr,lf,Y,Qk):

    #### compute weights using scalar parameters as inputs
    #### outputs: covariance and mean weights, and lambda_ (scalar parameter)
    Wc,Wm,lambda_=gen_weights(n,kap,alp,bet)

    #### compute sigma points using order, state covariance, state mean, and lambda_
    sigmas = np.copy(gen_sigma(n,Pk,Xk,lambda_))

    #### PASS SIGMAS THRU NONLINEAR FUNCTION F to obtain sigmas in function space
    for k in range(2*n+1):
        Y[k,:] = np.copy(motion_model_f(sigmas[k,:],1./des_rate,Uk,lr,lf))
    
    #### GET MEAN x AND COVARIANCE USING Y
    #### To obtain x, we scale each row of Y by each scalar element in Wm
    #### To obtain Px, we scale each matrix by scalar element in Wc
    #### Also need to wrapToPi due to computing angle difference in yres
    #### yres := function space residual
    x = np.dot(Wm,Y)   
    Px = np.full((n,n),0.0)
    yres = Y-x
    for k in range(len(yres)):
        yres[k,2] = wrapToPi(yres[k,2])
    for k in range(2*n+1):
        if k == 0:
            Px = Wc[k] * np.outer(yres[k],yres[k])
        else:
            Px += Wc[k] * np.outer(yres[k],yres[k])
    #### Add process noise covariance Qk
    Px+= Qk
    return Px,x,Wc,Wm,Y

def kalman_update(H,Zm,R,Pk,x,Y,Wm,Wc):
   
    #### Convert sigmas from function space Y into measurement space Z
    #### and then find the estimated measurement mean z 
    #### To find z, we scales each row of Z by an element in Wm, the mean weights
    Z = H.dot(Y.T).T
    z = np.dot(Wm,Z)    
    
    #### COMPUTE ALL RESIDUALS
    #### 1:= sensor measurement residual
    #### 2:= estimated measurement residual
    #### 3:= function space residual
    y = (Zm.reshape(z.shape) - z)
    zres = Z - z
    yres = Y - x
    
    #### NORMALIZE INNOVATION ANGLE RESIDUALS IF THEY EXIST
    #### BY SEARCHING THROUGH H:=(STATE TO MEASUREMENT SPACE MAPPING) and saving the row index
    Hind = []
    for k in range(len(H)):
        #### RHS corresponds to yaw measurements
        if H[k][2]==1:
            y[k] = wrapToPi(y[k])
            Hind.append(k)
            #if np.degrees(Zm[k]) > 175.0:
            #    pdb.set_trace()
            Stest.states[11] = np.degrees(y[k])
    
    #### NORMALIZE ALL ANGLE RESIDUALS
    for k in range(len(zres)):
        if len(Hind)>0:
            for h in Hind:
                zres[k,h] = wrapToPi(zres[k,h])
        yres[k,2] = wrapToPi(yres[k,2])


    #### COMPUTE COVARIANCE OF MEASUREMENTS WITH MEAS. SPACE AND MEAN
    #### We scales each matrix by an element in Wc, Covariance weighting
    for k in range(2*n+1):
        if k==0:
            Pz = Wc[k] * np.outer(zres[k],zres[k])
        else:
            Pz += Wc[k] * np.outer(zres[k],zres[k])   
    Pz += R 

    #### Compute cross covariance of state and measurements
    for k in range(2*n+1):
        if k == 0:
            Pxz = Wc[k]*np.outer(yres[k],zres[k])
        else:
            Pxz+= Wc[k]*np.outer(yres[k],zres[k])
    
    ####COMPUTE Kalman Gain
    K = Pxz.dot(np.linalg.inv(Pz))
    
    #### UPDATE STATE MEAN AND COVARIANCE
    #### attempt to wrap angle here
    Xk = x + K.dot(y)
    Xk[2] = wrapToPi(Xk[2])
    Pk = Px - K.dot(Pz).dot(K.T)
    return Xk,Pk

def jointstate_cb(msg):
    global wvsteer,wvaccd, wvspeed, wvyawd, jointstate_cbTime, vel_wheels_prev, vel_wheels_dot_prev,cb_flag
    #JointState
    R = 0.2413#wheel radius
    W = 1.58242#m track width
    L = 2.84988#m wheel base
    #### OBTAIN MESSAGES
    jointstate_cbTime[0] = jointstate_cbTime[1]
    jointstate_cbTime[1] = msg.header.stamp.secs + msg.header.stamp.nsecs*np.power(10.0,-9.0)
    steer_flp = msg.position[4] #steer_fl_phi
    steer_frp = msg.position[5] #steer_fr_phi
    #### VEL_WHEELS CALCULATIONS
    vel_wheels = float(sum([msg.velocity[2],msg.velocity[3]]))/2.0  #average vehicle speed rear tires
    dt = jointstate_cbTime[1] - jointstate_cbTime[0]
    vel_wheels_dot = 0.45*(vel_wheels - vel_wheels_prev)/dt + 0.55*vel_wheels_dot_prev
    vel_wheels_prev = vel_wheels
    vel_wheels_dot_prev = vel_wheels_dot
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
    offset = 1.
    wvspeed = veh_speed + offset #TO offset bias
    wvyawd = thetaD
    wvaccd = vel_wheels_dot*R #xdotdot = R*phidotdot
    wvsteer = np.arctan(tanphi)

    #print "Steer angle: ",wvsteer,phi_i,phi_o,"Accel desired: ", wvaccd
    #### SET CALLBACK FLAG TO 1
    cb_flag[3] += 1
    cb_flag[3] = clamp_reset(cb_flag[3],1000,1)
    cb_flag[6] += 1
    cb_flag[6] = clamp_reset(cb_flag[6],1000,1)

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
    cb_flag[4] += 1
    cb_flag[4] = clamp_reset(cb_flag[4],1000,1)
    cb_flag[7] += 1
    cb_flag[7] = clamp_reset(cb_flag[7],1000,1)
def mtiori_cb(msg):
    global mtiOP, cb_flag
    #orientationEstimate
    mtiOP=np.radians(msg.yaw) #Orientation Psi
    cb_flag[5] += 1
    cb_flag[5] = clamp_reset(cb_flag[5],1000,1)

def vehgps_cb(msg):
    global vehUTMx, vehUTMy, cb_flag
    #NavSatFix
    #temp = utm.from_latlon(msg.latitude, msg.longitude)    #from utm module
    temp = utm.fromLatLong(msg.latitude, msg.longitude)      #from geodesy.utm
    vehUTMx = temp.easting#temp[0]
    vehUTMy = temp.northing#temp[1]
    cb_flag[0] += 1
    cb_flag[0] = clamp_reset(cb_flag[0],1000,1)

    cb_flag[1] += 1
    cb_flag[1] = clamp_reset(cb_flag[1],1000,1)
def mtigps_cb(msg):
    global mtiUTMx, mtiUTMy, cb_flag
    #positionEstimate
    #temp = utm.from_latlon(msg.latitude, msg.longitude)    #using utm module
    temp = utm.fromLatLong(msg.latitude, msg.longitude)     #using geodesy.utm
    mtiUTMx = temp.easting#temp[0]
    mtiUTMy = temp.northing#temp[1]
    cb_flag[8] += 1
    cb_flag[8] = clamp_reset(cb_flag[8],1000,1)

    cb_flag[9] += 1
    cb_flag[9] = clamp_reset(cb_flag[9],1000,1)
def vehacc_cb(msg):
    global vehAcc, cb_flag
    #float64
    #### THIS IS VEHICLE STATE NOT COMMAND!!!
    vehAcc = msg.data			
    #cb_flag[5] += 1
    #cb_flag[5] = clamp_reset(cb_flag[5],1000,1)

def vehgpsvel_cb(msg):
    global vehgpsVel
    #TwistStamped
    vehgpsVel = np.sqrt(msg.twist.linear.x**2 + msg.twist.linear.y**2)
    cb_flag[2] += 1
    cb_flag[2] = clamp_reset(cb_flag[2],1000,1)


#### NOT FUSING THIS ONE BECAUSE IT'S CHEATING
def piksigps_cb(msg):
    global piksiUTMx, piksiUTMy
    temp = utm.fromLatLong(msg.latitude,msg.longitude)
    piksiUTMx = temp.easting
    piksiUTMy = temp.northing

def clamp_reset(a,b,c):
    if a >=b:
        return c
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
    x[5] = u[1]                             #vdot
    x[4] = x[5]*dt+x[4]                     #v
    x[3] = x[4]/lr*np.sin(bet)              #yawdot
    x[2] = x[3]*dt+x[2]  #yaw: Need to wrap #don't need to do this?
    x[1] = x[4]*np.sin(x[2]+bet)*dt+x[1]    #y
    x[0] = x[4]*np.cos(x[2]+bet)*dt+x[0]    #x

    return x

def wrapToPi(x):
    x = np.fmod(x+np.pi,2*np.pi)
    if x<0:
        x+= 2*np.pi
    return x-np.pi

def angleDiff(x):
    x = np.fmod(x,2*np.pi)
    if x>=np.pi:
        x-= 2*np.pi
    elif x<=-np.pi:
        x+= 2*np.pi
    return x

def measure_model_h(x):
    # 10 different quantities
    #vehGps X
    #vehGps Y
    #vehGps SpeedBody
    #vehWheel SpeedBody
    #vehImu AccelBody
    #mtiImuFiltered Yaw
    #vehWheel YawRate
    #vehImu YawRate
    H=np.array([[1,0,0,0,0,0],
             [0,1,0,0,0,0],
             [0,0,0,0,1,0], #vehgpsVel
             [0,0,0,0,1,0],
             [0,0,0,0,0,1], ####change last place back to 1
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
#rospy.Timer(rospy.Duration(2), timer_cb)

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

des_rate = 30.0
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
kap = 3 #6 is the order...change 3 to 7?
Xk = np.zeros(n)     #mean of states?
Px = np.full((n,n),0.0)      #covariance associated with x, the mean states
Pz = np.full((nz,nz),0.0)    #covariance associated with z, the mean measurements
Pxz = np.full((n,nz),0.0)    #cross covariance associated with x and z
Y = np.full((2*n+1,n),0.0)  #result of passing Sigmas through function F

timeInit = rospy.get_time()
print "READY TO BEGIN, START ROSBAG"
while not rospy.is_shutdown():
    if (sum(cb_flag )>= 5) and init_flag==1:

        #### CALCULATE UK VECTOR
        Ukdot = np.array([(wvsteer-wvsteer_prev)*des_rate, (wvaccd-cmdAccprev)*des_rate])
        cmdAccprev = wvaccd
        wvsteer_prev = wvsteer
        #Uk = np.array([wvsteer,wvaccd])   #steering tire angle, desired veh acceleration
        Uk = np.array([wvsteer,vimuVD])   #steering tire angle, desired veh acceleration
        "PREDICTION PHASE" #Returns the mean and covariance
        Px,x,Wc,Wm,Y=kalman_predict(n,kap,alp,bet,Pk,Xk,des_rate,Uk,lr,lf,Y,Qk)
        " UPDATE PHASE"
        #### Pass new sigma point, Y, through function h, which converts them to points in mesurement space. I        have nz measuredments
        
        #### GET MEASUREMENT VALUES FROM SENSORS
        #1=VGPSx, 2=VGPSy, 3=VGPSVEL, 4=VWheelVEL, 5=VIMUAccel, 6=MTIYaw, 7=VWheelYawD, 8=VImuYAWD,9=MTIGPSx, 10=MTIGPSy
        Zm = np.array([vehUTMx-initPosX+X0[0],
                       vehUTMy-initPosY+X0[1],
                       vehgpsVel,
                       wvspeed,
                       vimuVD,
                       mtiOP,
                       wvyawd,
                       vimuPD,
                       mtiUTMx-initPosX+X0[0],
                       mtiUTMy-initPosY+X0[1]])
    
        # Test without using MTI gps for now
        #1=VGPSx, 2=VGPSy, 3=VGPSVEL, 4=VWheelVEL, 5=VIMUAccel, 6=MTIYaw, 7=VWheelYawD, 8=VImuYAWD,9=MTIGPSx, 10=MTIGPSy
        Hnum=np.array([[1,0,0,0,0,0],
                         [0,1,0,0,0,0],
                         [0,0,0,0,1,0], #vehgpsVel
                         [0,0,0,0,1,0],
                         [0,0,0,0,0,1], ####change last place back to 1
                         [0,0,1,0,0,0],
                         [0,0,0,1,0,0],
                         [0,0,0,1,0,0],
                         [1,0,0,0,0,0],
                         [0,1,0,0,0,0]])


        #### UPDATE ONLY ON CB CHANGE
        cb_diff = list(np.asarray(cb_flag) - np.asarray(cb_flag_prev))
        cb_flag_prev = list(cb_flag)
        if sum(cb_diff)>0:
            Hnum_pass = []
            R_pass = []
            Z_pass = []
            for k in range(len(cb_diff)):
                if cb_diff[k] != 0:
                    Hnum_pass.append(Hnum[k])
                    R_pass.append(Rk[k][k])
                    Z_pass.append(Zm[k])
            if len(Hnum_pass) > 0:
                Hnum_pass = np.asarray(Hnum_pass)
                #print Hnum_pass, cb_diff,"\n"
                R_pass = np.diag(R_pass)
                Z_pass = np.asarray(Z_pass).reshape(len(Z_pass),1)
                Xk,Pk=kalman_update(Hnum_pass,Z_pass,R_pass,Px,x,Y,Wm,Wc)
        
        #### Try Savitzky golay filter
        #### create a flatten list out of lists of lists
        # TRY WITH x instead of Xk
        temp = [item for item in x]
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
            #pdb.set_trace()
            svresults = savgol_filter(temp,len_,svgorder,axis=0)
        #### TRY WITH x INSTEAD OF Xk
        for i in range(0,numStates):
            SI.states[i] = x[i]    
        pubinfo.publish(SI)
        #print x[2]*180.0/np.pi,mtiOP*180.0/np.pi

        Stest.states[0] = Zm[0]  #vehGps
        Stest.states[1] = Zm[1]  #vehGps
        Stest.states[2] = Zm[8]  #MtiGps            ##DELETED
        Stest.states[3] = Zm[9]  #MtiGps            ##DELETED
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
        Stest.states[10] = np.degrees(x[2])         # estimated yaw state
        Stest.states[12] = Xktest[2]
        Stest.states[13] = Xktest[3]
        Stest.states[14] = Xktest[4]
        Stest.states[15] = Xktest[5]
        Stest.states[16] = Zm[3]    #wvspeed
        Stest.states[17] = Zm[4]    #vehimuspeed
        Stest.states[18] = Zm[2]    #gps speed
        Stest.states[19] = Uk[0]    #Steer Command
        Stest.states[20] = Uk[1]    #Accel Command
        pubtest.publish(Stest)
    else:
        if init_flag == 0 and all(i >= 1 for i in cb_flag):
            #Set initial state
            initPosX = np.mean([vehUTMx])
            initPosY = np.mean([vehUTMy])
            #X0 = np.array([[7.3794], [-11.45], [mtiOP], [np.mean([wvyawd,vimuPD])], [np.mean([vehgpsVel,wvspeed])],[vimuVD]])
            X0 = np.array([0, 0, mtiOP, np.mean([wvyawd,vimuPD]), np.mean([vehgpsVel]),vimuVD])
            Xk = np.copy(X0)
            Xktest = np.copy(X0)
            for i in range(0,numStates):
                SI.states[i] = Xk[i] 
            cmdAccprev = wvaccd
            ####initial covariance matrices
            #initial covariance
            #Pk = 0.5*np.identity(6) #Initial State Covariance
            Pk = np.full((n,n),0)+50.0*np.identity(n)
            #Qk = 100.0*np.identity(n) #Process Noise covariance..we can scale this later by residual 500
            dt = 1./des_rate
            Qk = np.diag([0.5*10.0*dt**2,0.5*10.0*dt**2,0.5*1.0*dt**2,1.0*dt,10.0*dt,10.0])
            #Qk = np.diag([6,6,0.1,0.95,1.0,1.0])

            #1=VGPSx, 2=VGPSy, 3=VGPSVEL, 4=VWheelVEL, 5=VIMUAccel, 6=MTIYaw, 7=VWheelYawD, 8=VImuYAWD,9=MTIGPSx, 10=MTIGPSy
            Rk = np.diag([5.0,5.0,1.0,25.0,1.0,0.5,0.85,0.5,5.0,5.0])    #Measurement Noise Covariance 5
            #R = np.identity(nz)
            #R = np.diag([5.0,5.0,0.25,0.25,5.0,1.00,1.00,0.25,1.00,0.5])    #Measurement Noise Covariance
            #R = np.diag([1,1,0.05,0.05,1,0.045,0.20,0.002,0.045,0.01])    #Measurement Noise Covariance
            #TODO: if cb_flag is zero then set that measurement to something high!
            init_flag = 1
    rate.sleep()

	

