#!/usr/bin/env python
# THIS SCRIPT IMPLEMENTS AN EXTENDED KALMAN FILTER FOR THE LINCOLN VEHICLE
# WE USE A KINEMATIC MODEL BASED ON THE PAPER BY KONG
# UPDATE ONLY ON NEW DATA (CALLBACK) NOT IN WHILE LOOP
# PREDICTION RUNS IN THE WHILE LOOP
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

global wvsteer, wvaccd, wvspeed, wvyawd, jointstate_cbTime, vel_wheels, vel_wheels_prev, cb_flag
jointstate_cbTime = [0,0]
vel_wheels_prev = 0.0
vel_wheels = 0.0
global vimuAVC, vimuLAC
global mtiOP
global vehUTMx, vehUTMy
global mtiUTMx, mtiUTMy
global vehAcc
global vehgpsvel
global piksiUTMx, piksiUTMy
piksiUTMx = 0
piksiUTMy = 0

cb_flag = [0,0,0,0,0,0,0]
cb_flag_prev = list(cb_flag)
init_flag = 0

t = 0
tprev = 0
cmdAccprev = 0
wvsteer_prev=0
rospy.init_node("kalman_filter_testnode")
#### FUNCTION DEFINTIONS
def timer_cb(event):
    global cb_flag
    for i in range(0,len(cb_flag)):
        cb_flag[i] -= 1
        cb_flag[i] = clamp(cb_flag[i],6,0)

def kalman_predict(Fk,G,Pk,Qk,Xk,Uk):
    #Prediction
    #Xk = Fk.dot(Xk) + G.dot(Uk)
    Xk = Fk.dot(Xk)
    Pk = Fk.dot(Pk).dot(np.transpose(Fk))+Qk
    return Xk, Pk


def kalman_update(H,Zk,R,Pk,Xk):
    # Update
    S = H.dot(Pk).dot(H.T) + R        # Measurement Covariance
    if np.linalg.cond(S) > 1/sys.float_info.epsilon:
        print "MEASUREMENT COVARIANCE HAS BECOME SINGULAR"    
        pdb.set_trace()
    #Gk = (Pk.dot(np.transpose(H))).dot(np.linalg.inv(H.dot(Pk).dot(np.transpose(H))+R))
    Gk = (Pk.dot(np.transpose(H))).dot(np.linalg.inv(S))

    try:
        Xk = np.copy(Xk + Gk.dot((Zk - H.dot(Xk))))
        Pk = (np.eye(np.shape(Pk)[0])-Gk.dot(H)).dot(Pk)
    except:
        e = sys.exc_info()[0]
        write_to_page( "<p>Error: %s</p>" % e )
        pdb.set_trace()
    return Xk,Pk
    

def jointstate_cb(msg):
    global wvsteer,wvaccd, wvspeed, wvyawd, jointstate_cbTime, vel_wheels, vel_wheels_prev,cb_flag
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
    vel_wheels_prev = vel_wheels
    vel_wheels = 0.85*float(sum([msg.velocity[2],msg.velocity[3]]))/2.0 + 0.15*vel_wheels_prev  #average vehicle speed rear tires
    dt = jointstate_cbTime[1] - jointstate_cbTime[0]
    vel_wheels_dot = (vel_wheels - vel_wheels_prev)/dt 
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
    wvsteer = np.arctan(tanphi)
    #### SET CALLBACK FLAG TO 1
    cb_flag[0] += 1
    cb_flag[0] = clamp_reset(cb_flag[0],1000,1)

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
    cb_flag[1] = clamp_reset(cb_flag[1],1000,0)
    
def mtiori_cb(msg):
    global mtiOP, cb_flag

    #orientationEstimate
    mtiOP=np.radians(msg.yaw) #Orientation Psi
    cb_flag[2] += 1
    cb_flag[2] = clamp_reset(cb_flag[2],1000,0)
    #### After setting initial position
    

def vehgps_cb(msg):
    global vehUTMx, vehUTMy, cb_flag

    #NavSatFix
    #temp = utm.from_latlon(msg.latitude, msg.longitude)    #from utm module
    temp = utm.fromLatLong(msg.latitude, msg.longitude)      #from geodesy.utm
    vehUTMx = temp.easting#temp[0]
    vehUTMy = temp.northing#temp[1]
    cb_flag[3] += 1
    cb_flag[3] = clamp_reset(cb_flag[3],1000,0)
   
def mtigps_cb(msg):
    global mtiUTMx, mtiUTMy, cb_flag

    #positionEstimate
    #temp = utm.from_latlon(msg.latitude, msg.longitude)    #using utm module
    temp = utm.fromLatLong(msg.latitude, msg.longitude)     #using geodesy.utm
    mtiUTMx = temp.easting#temp[0]
    mtiUTMy = temp.northing#temp[1]
    cb_flag[4] += 1
    cb_flag[4] = clamp_reset(cb_flag[4],1000,0)

           
def vehacc_cb(msg):
    global vehAcc, cb_flag
    #float64
    #### THIS IS VEHICLE STATE NOT COMMAND!!!
    vehAcc = msg.data			
    cb_flag[5] += 1
    cb_flag[5] = clamp_reset(cb_flag[5],1000,0)

def vehgpsvel_cb(msg):
    global vehgpsVel, cb_flag

    #TwistStamped
    vehgpsVel = np.sqrt(msg.twist.linear.x**2 + msg.twist.linear.y**2)
    cb_flag[6] += 1
    cb_flag[6] = clamp_reset(cb_flag[6],1000,0)
   

#### NOT FUSING THIS ONE BECAUSE IT'S CHEATING
def piksigps_cb(msg):
    global piksiUTMx, piksiUTMy
    temp = utm.fromLatLong(msg.latitude,msg.longitude)
    piksiUTMx = temp.easting
    piksiUTMy = temp.northing

def clamp_reset(a,b,c):
    if a >=b:
        return c
    elif a < c:
        return c
    else:
        return a

#### ODOMETRY
rospy.Subscriber("test/vehicle/joint_states",JointState,jointstate_cb)
rospy.Subscriber("test/vehicle/imu/withcovariance",Imu, vehimu_cb)
rospy.Subscriber("test/mti/filter/orientation",orientationEstimate, mtiori_cb)
#### GPS
rospy.Subscriber("test/vehicle/gps/fix/withcovariance",NavSatFix, vehgps_cb)
rospy.Subscriber("test/mti/filter/position", positionEstimate, mtigps_cb)
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
Stest.states=10*[0]

lf = 1.3314298821150170049065764033003	
lr = 1.5185105279655697341212317041936
#base = [30.6286298219,-96.4821284934] 

#### LINEARIZE EQUATIONS USING JACOBIAN
sym.init_printing(use_latex=True)
myvarsF = sym.sympify('x(t), y(t), psi(t), diff(psi(t),t), v(t), diff(v(t),t)')
myvarsG = sym.sympify('delta_f(t), a(t)')
bet = sym.sympify('atan(l_r/(l_f+l_r)*tan(delta_f(t)))')
f1 = sym.sympify('v(t)*cos(psi(t)+bet)')   # xdot
f1 = f1.subs({'bet':bet})
f2 = sym.sympify('v(t)*sin(psi(t)+bet)')   # ydot
f2 = f2.subs({'bet':bet})
f3 = sym.sympify('v(t)/l_r*sin(bet)')  # psidot
f3 = f3.subs({'bet':bet})           #substitute bet expression
f4 = sym.diff(f3,sym.symbols('t'))  #psidotdot
f4 = f4.subs({'bet':bet})           #sub bet expression
f5 = sym.sympify('a(t)')            #vdot 
f6 = f5.diff('t')                   #vdotdot
J = sym.Matrix([f1,f2,f3,f4,f5,f6]) #define vector-valued function
F = J.jacobian(myvarsF)              # Perform jacobiam wrt to myvars
G = J.jacobian(myvarsG)

#### All "deriviatves" must be bare-bone variables

#### symbolic variables to be replaced with numerics
oldvars = ['l_f','l_r','diff(delta_f(t),t)','delta_f(t)','diff(v(t),t)','v(t)','diff(psi(t),t)','psi(t)','diff(a(t),t)','a(t)']   #10 vars to be replaced

newvars = ['l_f','l_r','dd','d','dv','v','dp','p','da','a']   #10 vars to be replaced
for i in range(0,len(oldvars)):
    temp = sym.sympify(oldvars[i])
    F = F.subs({temp:newvars[i]})
    G = G.subs({temp:newvars[i]})

#### NEED TO ADD I AND CONSIDER LOOP RATE. DISCRETE TIME REQUIREMENTS
des_rate = 500.0
rate = rospy.Rate(des_rate)
F = F/des_rate+sym.eye(F.shape[0])
G = G/des_rate

### All "deriviatves" must be bare-bone variables
Fobj = sym.lambdify((newvars),F,"numpy")
Gobj = sym.lambdify((newvars),G,"numpy")
count = 0           # ADAPTIVE KALMAN FILTER COUNTER
eps = 0             # ADAPTIVE KALMAN FILTER EPS
kfwin = 100         # SAVITZKY GOLAY, KF WIND
svgwin = 25         # SVG WIN
svgorder = 3
kfvect = []         # kf vector of size win
svresults=np.array([kfwin*[0]])
timeInit = rospy.get_time()

#### TESTING PURPOSES
Stest = States()
Stest.states=21*[0]

print "READY TO BEGIN, START ROSBAG"
while not rospy.is_shutdown():
    #TODO: Switch the else and if statement parts so the initializer is on top
    #TODO: PUT CALLBACK flags in here (while loop) and perform update for those with nonzero counter
    #       differences 
    if (sum(cb_flag )>= 5) and init_flag==1:
        #### SHIFT U FOR UPPER CASE
        #### CALCULATE DERIVATIVES OF INPUT AND INPUT ITSELF
        #t = rospy.get_time() - timeInit
        #Ukdot = np.array([[(cmdSteeringF(t)-cmdSteeringF(tprev))*50], [(wvaccd-cmdAccprev)*50]])
        #Ukdot = np.array([[(wvsteer-wvsteer_prev)*des_rate], [(wvaccd-cmdAccprev)*des_rate]])
        Ukdot = np.array([[(wvsteer-wvsteer_prev)*des_rate], [(vimuVD-cmdAccprev)*des_rate]])
        #tprev = t
        cmdAccprev = vimuVD
        wvsteer_prev = wvsteer
        #Uk = np.array([[cmdSteeringF(t)],[wvaccd]])   #steering tire angle, desired veh acceleration
        #Uk = np.array([[wvsteer],[wvaccd]])   #steering tire angle, desired veh acceleration
        Uk = np.array([[wvsteer],[vimuVD]])   #steering tire angle, desired veh acceleration
        numvars = [lf,lr,Ukdot[0][0],Uk[0][0],Xk[5][0],Xk[4][0],Xk[3][0],Xk[2][0],Ukdot[1][0],Uk[1][0]]#use est states, not measured
        #pdb.set_trace()

        "The * operator unpacks the list and passes them as positional arguments to whatever accepts them"

        Fnum = Fobj(*numvars)
        Gnum = Gobj(*numvars)

        Fnum = np.array(Fnum).astype(np.float64)
        Gnum = np.array(Gnum).astype(np.float64)
        Xk,Pk = kalman_predict(Fnum,Gnum,Pk,Qk,Xk,Uk)

        # xgps, ygps, xmgps, ymgps, vehgpsvel, vwo, vimudot, psimti, psiwodot, psiimudot 
        Hnum = np.array([[1,0,0,0,0,0],
                         [0,1,0,0,0,0],
                         [1,0,0,0,0,0],
                         [0,1,0,0,0,0],
                         [0,0,0,0,1,0], #vehgpsVel
                         [0,0,0,0,1,0],
                         [0,0,0,0,0,1],
                         [0,0,1,0,0,0],
                         [0,0,0,1,0,0],
                         [0,0,0,1,0,0]])
        Zk = np.array([[vehUTMx-initPosX+X0[0][0]],
                       [vehUTMy-initPosY+X0[1][0]],
                       [mtiUTMx-initPosX+X0[0][0]],
                       [mtiUTMy-initPosY+X0[1][0]],
                       [vehgpsVel],
                       [wvspeed],
                       [vimuVD],
                       [mtiOP],
                       [wvyawd],
                       [vimuPD]])

        #Xk,Pk=kalman(Fnum,Gnum,Hnum,Zk,Qk,R,Pk,Xk,Uk)
        ####Try Adapative Process Noise technique
        #y = Zk - Hnum.dot(Xk)                        # Compute Residual
        #S = Hnum.dot(Pk).dot(Hnum.T) + R        # Measurement Covariance
        ##if np.linalg.cond(S) < 1/sys.float_info.epsilon:
        #eps = np.dot(y.T,np.linalg.inv(S)).dot(y)       #Bar-Shalom
        #if eps > 36 and count < 6:
        #    Qk *= 10
        #    count += 1
        #elif eps <= 36 and count > 0:
        #    Qk /= 10
        #    count -=1
        #print "PK: ", Pk,"Pk norm: ",np.linalg.norm(Pk), "eps: ",eps[0][0]
        #### Try Savitzky golay filter
        #### create a flatten list out of lists of lists
        temp = [item[0] for item in Xk]
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
       
        #### UPDATE ONLY ON CB CHANGE
        cb_diff = list(np.asarray(cb_flag) - np.asarray(cb_flag_prev))
        cb_flag_prev = list(cb_flag)
        #print cb_diff
        if sum(cb_diff)>0:
            Hnum_pass = []
            R_pass = []
            Z_pass = []
            for k in range(len(cb_diff)):
                if cb_diff[k] != 0:
                    Hnum_pass.append(Hnum[k])
                    R_pass.append(R[k][k])
                    Z_pass.append(Zk[k][0])
            if len(Hnum_pass) > 0:
                Hnum_pass = np.asarray(Hnum_pass)
                R_pass = np.diag(R_pass)
                Z_pass = np.asarray(Z_pass).reshape(len(Z_pass),1)
                Xk,Pk=kalman_update(Hnum_pass,Z_pass,R_pass,Pk,Xk)
        
        for i in range(0,numStates):
            SI.states[i] = Xk[i][0]     
        pubinfo.publish(SI)

        Stest.states[0] = vehUTMx-initPosX+X0[0][0]  #vehGps
        Stest.states[1] = vehUTMy-initPosY+X0[1][0]  #vehGps
        Stest.states[2] = mtiUTMx-initPosX+X0[0][0] #MtiGps
        Stest.states[3] = mtiUTMy-initPosY+X0[1][0]  #MtiGps
        Stest.states[4] = piksiUTMx-initPosX+X0[0][0]  #Piksi
        Stest.states[5] = piksiUTMy-initPosY+X0[1][0]  #Piksi
        #pdb.set_trace()
        len1_ = np.ceil(len(svresults)/3.0)
        len2_ = np.ceil(len(svresults)/10.0)
        len3_ = np.ceil(len(svresults)/5.0)
        Stest.states[6] = svresults[int(-len1_)][0]  #svgol filter posX
        Stest.states[7] = svresults[int(-len1_)][1]  #svgol filter posY
        Stest.states[8] = svresults[-1][3]          #svgol filter yawrate
        Stest.states[9] = svresults[int(-len3_)][4]          #svgol filter velocity
        Stest.states[10] = Xk[2]*180.0/np.pi
        Stest.states[11] = vimuPD
        Stest.states[12] = wvyawd
        Stest.states[18] = vehgpsVel
        Stest.states[16] = wvspeed

        pubtest.publish(Stest)
    else:
        if init_flag == 0 and all(i > 1 for i in cb_flag):
            #Set initial state
            initPosX = np.mean([vehUTMx,mtiUTMx])
            initPosY = np.mean([vehUTMy,mtiUTMy])
            #X0 = np.array([[7.3794], [-11.45], [mtiOP], [np.mean([wvyawd,vimuPD])], [np.mean([vehgpsVel,wvspeed])],[vimuVD]])
            X0 = np.array([[0], [0], [mtiOP], [np.mean([wvyawd,vimuPD])], [np.mean([vehgpsVel,wvspeed])],[vimuVD]])
            Xk = np.copy(X0)
            for i in range(0,numStates):
                SI.states[i] = Xk[i][0] 
            cmdAccprev = wvaccd
            ####initial covariance matrices
            #initial covariance
            Pk = 10000.0*np.identity(6) #Initial State Covariance
            #Qk = 20.0*np.identity(6) #Process Noise covariance..we can scale this later by residual 500
            dt = 1./des_rate
            Qk = np.diag([0.5*8.8*dt**2,0.5*8.8*dt**2,0.5*0.5*dt**2,0.5*dt,8.8*dt,8.8])
            #Gtest = np.array([[0.5/(des_rate)**2],[0.5/(des_rate)**2],[0.5/(des_rate)**2],[1./des_rate],[1./des_rate],[1]])
            #Qk = Gtest.dot(Gtest.T)
            # xgps, ygps, xmgps, ymgps, vehgpsvel, vwo, vimudot, psimti, psiwodot, psiimudot 
            #R = 2.0*np.diag([5.0,5.0,5.0,5.0,1.0,3.0,2.0,0.8,0.5,0.5])    #Measurement Noise Covariance 5
            #R = np.identity(10)
            R = np.diag([5.0,5.0,0.05,0.05,0.5,5.0,0.05,0.025,0.0005,0.0008])    #Measurement Noise Covariance
            #R = np.diag([1,1,0.05,0.05,1,0.045,0.20,0.002,0.045,0.01])    #Measurement Noise Covariance
            #TODO: if cb_flag is zero then set that measurement to something high!
            init_flag = 1
    rate.sleep()

	

