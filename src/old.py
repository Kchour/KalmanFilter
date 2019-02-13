#!/usr/bin/env python
# THIS SCRIPT WILL BE USED TO DETERMINE Q, R NOISE COVARIANCE MATRICES
#  AND TO PARSE BAG FILES
import rospy
import rosbag
import math
from geometry_msgs.msg import Point32
from geometry_msgs.msg import PointStamped
from nav_msgs.msg import Odometry
from sensor_msgs.msg import NavSatFix
from sensor_msgs.msg import JointState
from custom_msgs.msg import positionEstimate
from custom_msgs.msg import orientationEstimate
from scipy.interpolate import interp1d
import numpy as np
import pdb
import tf
import utm
import sympy as sym




#### IMPORT TROUBLESOME TOPICS FIRST
#IMPORT STEERING_REPORT BECAUSE THE BAG MESSAGE IS OLDER THAN INSTALLEDTO
filename = "./2018-12-21-12-09-33.bag" 
cmdSteering = [[0,0]]
bag = rosbag.Bag(filename)
startTime = bag.get_start_time()
endTime = bag.get_end_time()
duration = endTime - startTime
for topic, msg, t in bag.read_messages(topics=['/vehicle/steering_report']):
    cmdSteering.append([t.to_sec()-startTime,msg.steering_wheel_angle_cmd])
#pdb.set_trace()
cmdSteering.append([duration,0])
cmdSteering = np.asarray(cmdSteering)
bag.close()

#### INTERPOLATION OBJECT
# TO CALL DO CMDSTEERINGF(SOME X)
# DONT FORGET TO DIVIDE BY STEERING RATIO := 14.8
cmdSteeringF = interp1d(cmdSteering[:,0],cmdSteering[:,1]/14.8)

global wvaccd, wvspeed, wvyawd, jointstate_cbTime, vel_wheels_prev
jointstate_cbTime = [0,0]
vel_wheels_prev = 0
global vimuAVC, vimuLAC
global mtiOP
global vehUTMx, vehUTMy
global mtiUTMx, mtiUTMy
global vehAcc
cb_flag = [0,0,0,0,0,0]
init_flag = 0

t = 0
tprev = 0
cmdAccprev = 0
rospy.init_node("kalman_filter_testnode")
#### FUNCTION DEFINTIONS
def timer_cb(event):
    global cb_flag
    for i in len(cb_flag):
        cb_flag[i] -= 1
        cb_flag[i] = clamp(cb_flag[i],5,0)

def kalman(Fk,G,H,Zk,Qk,R,Pk,Xk,Uk):
    #### GATHER NUMERIC LIST FOR SUBSTITUTION (ALLVARS)
    # Prediction
    Xk = Fk*Xk + G*Uk
    Pk = Fk*Pk*np.transpose(Fk)
    # Update
    Gk = Pk*np.transpose(H)*(H*Pk*np.transpose(H)+R)^-1
    Xk = Xk + Gk*(Zk - H*Xk)
    Pk = (eye(np.shape(Pk)[0])-Gk*H)*Pk
    return Xk,Pk
    

def jointstate_cb(msg):
    global wvaccd, wvspeed, wvyawd, jointstate_cbTime, vel_wheels_prev,cb_flag
    #JointState
    R = 0.2413#wheel radius
    W = 1.58242#m track width
    L = 2.84988#m wheel base
    #### OBTAIN MESSAGES
    jointstate_cbTime[0] = jointstate_cbTime[1]
    jointstate_cbTime[1] = msg.header.stamp.secs + msg.header.stamp.nsecs*10^-9
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
    wvspeed = veh_speed
    wvyawd = thetaD
    wvaccd = vel_wheels_dot*R #xdotdot = R*phidotdot
    #### SET CALLBACK FLAG TO 1
    cb_flag[0] += 1
    cb_flag[0] = clamp(cb_flag[0],5,0)

def vehimu_cb(msg):
    global vimuAVC, vimuLAC, vimuPD, vimuVD, cb_flag
    #Imu, only has ang vel,cov, lin accel,cov
    #row major x,y,z covariances
    vimuPD = msg.angular_velocity.z #psi dot
    vimuVD = msg.linear_acceleration.x #vdot

    temp = msg.angular_acceleration_covariance
    vimuAVC = np.reshape(temp,(3,3))    #angular velocity covariance rpy

    temp = msg.linear_acceleration_covariance
    vimuLAC = np.reshape(temp,(3,3))    #linear acceleration covariance xyz
    cb_flag[1] += 1
    cb_flag[1] = clamp(cb_flag[1],5,0)

def mtiori_cb(msg):
    global mtiOP
    #orientationEstimate
    mtiOP=np.radians(msg.yaw) #Orientation Psi
    cb_flag[2] += 1
    cb_flag[2] = clamp(cb_flag[2],5,0)

def vehgps_cb(msg):
    global vehUTMx, vehUTMy, cb_flag
    #NavSatFix
    temp = utm.from_latlon(msg.latitude, msg.longitude)
    vehUTMx = temp[0]
    vehUTMy = temp[1]
    cb_flag[3] += 1
    cb_flag[3] = clamp(cb_flag[3],5,0)

def mtigps_cb(msg):
    global mtiUTMx, mtiUTMy, cb_flag
    #positionEstimate
    temp = utm.from_latlon(msg.latitude, msg.longitude)
    mtiUTMx = temp[0]
    mtiUTMy = temp[1]
    cb_flag[4] += 1
    cb_flag[4] = clamp(cb_flag[4],5,0)

def cmdacc_cb(msg):
    global vehAcc, cb_flag
    #float64
    #### THIS IS VEHICLE STATE NOT COMMAND!!!
    vehAcc = msg.data			
    cb_flag[5] += 1
    cb_flag[5] = clamp(cb_flag[5],5,0)
    
def clamp(a,b,c):
    if a >b:
        return b
    elif a < c:
        return c

#### ODOMETRY
rospy.Subscriber("test/vehicle/joint_states",JointState,jointstate_cb)
rospy.Subscriber("test/vehicle/imu/withcovariance",Imu, vehimu_cb)
rospy.Subscriber("test/mti/filter/orientation",orientationEstimate, mtiori_cb)
#### GPS
rospy.Subscriber("test/vehicle/gps/fix/withcovariance",NavSatFix, vehgps_cb)
rospy.Subscriber("test/mti/filter/position", positionEstimate, mtigps_cb)
#### CONTROL INPUTS
rospy.Subscriber("test/vehicle/filtered_accel",Float64,cmdacc_cb)
rospy.Timer(rospy.Duration(1), timer_cb)

lf = 1.3314298821150170049065764033003	
lr = 1.5185105279655697341212317041936
#base = [30.6286298219,-96.4821284934] 

#### LINEARIZE EQUATIONS USING JACOBIAN
myvarsF = sym.sympify('x(t), y(t), psi(t), diff(psi(t),t), v(t), diff(v(t),t)')
myvarsG = sym.sympify('delta_f(t), a(t)')
bet = sym.sympify('atan(l_r/(l_f+l_r)*tan(delta_f(t)))')
f1 = sym.sympify('v(t)*cos(psi(t)+bet)')   # xdot
f1 = f1.subs({'bet':bet})
f2 = sym.sympify('v(t)*sin(psi(t)+bet)')   # ydot
f2 = f2.subs({'bet':bet})
f3 = sym.sympify('v(t)/l_r*sin(bet)')  # psidot
f3 = f3.subs({'bet':bet})           #substitute bet expression
f4 = sym.sympify('a(t)')            #vdot 
f5 = sym.diff(f3,sym.symbols('t'))  #psidotdot
f5 = f5.subs({'bet':bet})           #sub bet expression
J = sym.Matrix([f1,f2,f3,f4,f5])    #define vector-valued function
F = J.jacobian(myvarsF)              # Perform jacobiam wrt to myvars
G = J.jacobian(myvarsG)

#### symbolic variables to be replaced with numerics
allvars = ['l_f','l_r','delta_f(t)','diff(delta_f(t),t)','v(t)','diff(v(t),t)','psi(t)','diff(psi(t),t)','a(t)']   #9 vars to be replaced

timeInit = rospy.Time.now()
rate = rospy.Rate(50)
while not rospy.is_shutdown():
    if sum(cb_flag >= 20):
        #### SHIFT U FOR UPPER CASE
        #### CALCULATE DERIVATIVES OF INPUT AND INPUT ITSELF
        tprev = t
        t = rospy.Time.now() - timeInit
        Ukdot = [(cmdSteeringF(t)-cmdSteeringF(tprev))/50, (wvaccd-cmdAccprev)/50] 
        cmdAccprev = wvaccd
        Uk = [cmdSteeringF(t),wvaccd]   #steering tire angle, desired veh acceleration
        #numvars = [lf,lr,Uk[0],Ukdot[0],wvspeed,vehAcc,mtiOP,wvyawd,wvaccd]#use est states, not measured
        numvars = [lf,lr,Uk[0],Ukdot[0],x[4],x[5],x[2],x[3],Uk[1]]#use estimated states, not measured
        #### SUB VARS 1 BY 1 IN F
        for i in allvars:
            temp = sym.sympify(allvars[i])
            Fnum = F.subs({temp:numvars[i]})
            Gnum = G.subs({temp:numvars[i]})
        Fnum = np.array(Fnum).astype(np.float64)
        Gnum = np.array(Gnum).astype(np.float64)
        Hnum = np.array([[1,0,0,0,0,0],
                         [0,1,0,0,0,0],
                         [1,0,0,0,0,0],
                         [0,1,0,0,0,0],
                         [0,0,0,0,1,0],
                         [0,0,0,0,0,1],
                         [0,0,1,0,0,0],
                         [0,0,0,1,0,0],
                         [0,0,0,1,0,0]])
        Zk = np.array([[vehUTMx-initPosX],
                       [vehUTMy-initPosY],
                       [mtiUTMx-initPosX],
                       [mtiUTMy-initPosY],
                       [wvspeed],
                       [vimuVD],
                       [mtiOP],
                       [wvyawd],
                       [vimuPD]])
        Xk,Pk=kalman(Fnum,Gnum,Hnum,Zk,Qk,R,Pk,Xk,Uk)
        init_flag = 1
    else:
        if init_flag == 0 and all(i > 0 for i in cb_flag):
            #Set initial state
            initPosX = np.mean(vehUTMx,mtiUTMx)
            initPosY = np.mean(vehUTMy,mtiUTMy)
            X0 = np.array([[0], [0], [mtiOP], [wvyawd], [wvspeed]])
            Xk = np.copy(X0)
            cmdAccprev = wvaccd
            ####initial covariance matrices
            #initial covariance
            Pk = np.identity(5) #5 states
            Qk = np.identity(5)
            # xgps, ygps, xmgps, ymgps, vwo, vimudot, psimti, psiwodot, psiimudot 
            R = np.diag([5.0,5.0,0.25,0.25,0.45,1.00,0.02,0.45,0.5])
    rate.sleep()

	

