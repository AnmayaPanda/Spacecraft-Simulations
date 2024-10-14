'''1D PROJECTILE DYNAMICS WITH NEWTONIAN GRAVITY'''
###Import all the modules
import numpy as np ###numeric python
import matplotlib.pyplot as plt ###matlab style plotting
import scipy.integrate as sci ##integration toolbox

plt.close('all')

##DEFINE SOME CONSTANT PARAMETERS
G = 6.6742*10**-11 ##Gravitational Constant (SI Units)

##PLANET
##EARTH
Rplanet = 6371000 ##meters
mplanet = 5.972e24 #kg

##KERBIN
#Rplanet = 600000 #meters
#mplanet = 5.2915158*10**22 #kg

##ROCKET
mass = 640.0/1000.0 ##kg

##Gravitational Accelerattion Model
def gravity(z):
    global Rplanet,mplanet

    r = np.sqrt(z**2) ##Norm from centre of planet to spacecraft
    if r<Rplanet:
        accel = 0.0
    else:
        accel = G*mplanet/(Rplanet**3)*r
    return accel

###Equations of Motion
###F=m*a=m*zddot
## z is the altitude of the surface
## this is in meter
## zdot is the velocity
## zddot is the acceleration
###Second Order Differential Equation
def Derivatives(state,t):
    ###Globals
    global mass
    #state vector 
    z = state[0]
    velz = state[1]

    #Compute zdot - Kinematic Relationship
    zdot = velz

    ###Compute the total Forces
    ###Gravity
    gravityF = -gravity(z)*mass

    ###Aerodynamics
    aeroF = 0.0

    ###Thrust
    thrustF = 0.0

    Forces = gravityF+aeroF+thrustF

    #Compute Acceleration
    zddot = Forces/mass

    #Compute the statedot
    statedot = np.asarray([zdot,zddot])

    return statedot

#############EVERYTHING BELOW HERE IS THE MAIN SCRIPT###################

###Test Surface Gravity
print('Surface Gravity (m/s^2)',gravity(Rplanet))

###Initial Conditions
z0 = Rplanet ##m
velz0 = 33*331.0 ##m/s
stateinitial = np.asarray([z0,velz0])

##Time Window
tout = np.linspace(0,35,1000)

##Numerical Integration Call
stateout = sci.odeint(Derivatives,stateinitial,tout)

##Rename Variables
zout = stateout[:,0]
altitude = zout - Rplanet
velzout = stateout[:,1]

##Plot

###ALTITUDE
plt.plot(tout,altitude)
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')
plt.grid()
plt.show()

###VELOCITY
plt.figure
plt.plot(tout,velzout)
plt.xlabel('Time (sec)')
plt.ylabel('Normal Speed (m/s)')
plt.grid()
plt.show()