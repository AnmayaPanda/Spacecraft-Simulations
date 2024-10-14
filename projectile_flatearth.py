'''SIMULATING A PRJECTILE ON A FLAT EARTH'''
###Import all the modules
import numpy as np ###numeric python
import matplotlib.pyplot as plt ###matlab style plotting
import scipy.integrate as sci ##integration toolbox

##DEFINE SOME CONSTANT PARAMETERS

##ROCKET
mass = 640.0/1000.0 ##kg

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
    gravity = -9.81*mass

    ###Aerodynamics
    aero = 0.0

    ###Thrust
    thrust = 0.0

    Forces = gravity+aero+thrust

    #Compute Acceleration
    zddot = Forces/mass

    #Compute the statedot
    statedot = np.asarray([zdot,zddot])

    return statedot

#############EVERYTHING BELOW HERE IS THE MAIN SCRIPT###################

###Initial Conditions
z0 = 0.0 ##m
velz0 = 164.0 ##m/s
stateinitial = np.asarray([z0,velz0])

##Time Window
tout = np.linspace(0,35,1000)

##Numerical Integration Call
stateout = sci.odeint(Derivatives,stateinitial,tout)

##Rename Variables
zout = stateout[:,0]
velzout = stateout[:,1]

##Plot

###ALTITUDE
plt.plot(tout,zout)
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