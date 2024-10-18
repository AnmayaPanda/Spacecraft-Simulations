'''NUMERICAL INTEGRATION OF 2D ORBIT'''
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
def gravity(x,z):
    global Rplanet,mplanet

    r = np.sqrt(x**2 + z**2) ##Norm from centre of planet to spacecraft
    if r<Rplanet:
        accelx = 0.0
        accelz = 0.0
    else:
        accelx = G*mplanet/(r**3)*x
        accelz = G*mplanet/(r**3)*z
    return np.asarray([accelx,accelz])

###Equations of Motion
###F=m*a=m*zddot
## z is the altitude from the center of the planet along the north pole
## x is the altitude from the center of the planet along the equator through Africa
## this is in meter
## zdot is the velocity along z
## zddot is the acceleration along z
## zdot is the velocity along x
## xddot is the acceleration along x
###Second Order Differential Equation
def Derivatives(state,t):
    ###Globals
    global mass
    #state vector 
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]

    #Compute xdot and zdot - Kinematic Relationship
    xdot = velx
    zdot = velz

    ###Compute the total Forces
    ###Gravity
    gravityF = -gravity(x,z)*mass

    ###Aerodynamics
    aeroF = np.asarray([0.0,0.0])

    ###Thrust
    thrustF = np.asarray([0.0,0.0])

    Forces = gravityF+aeroF+thrustF

    #Compute Acceleration
    ddot = Forces/mass

    #Compute the statedot
    statedot = np.asarray([xdot,zdot,ddot[0],ddot[1]])

    return statedot

#############EVERYTHING BELOW HERE IS THE MAIN SCRIPT###################

###Test Surface Gravity
print('Surface Gravity (m/s^2)',gravity(0,Rplanet))

###Initial Conditions
x0 = Rplanet+600000 ##m
#To start at surface of earth x0 = Rplanet with some velz0 or vertical sppeed, say 100.0 m/s, it sort of crashes
z0 = 0.0
r0 = np.sqrt(x0**2+z0**2)
velx0 = 0.0
velz0 = np.sqrt(G*mplanet/r0)*1.1 ##m/s  multiplication factor to make it in elliptic orbit
stateinitial = np.asarray([x0,z0,velx0,velz0])

##Time Window
period = 2*np.pi/np.sqrt(G*mplanet)*r0**(3.0/2.0)*1.5 ##multiplication factor to move it in elliptic orbit
tout = np.linspace(0,period,1000)

##Numerical Integration Call
stateout = sci.odeint(Derivatives,stateinitial,tout)

##Rename Variables
xout = stateout[:,0]
zout = stateout[:,1]
altitude = np.sqrt(xout**2+zout*2) - Rplanet
velxout = stateout[:,2]
velzout = stateout[:,3]
velout = np.sqrt(velxout**2+velzout**2)

##Plot

###ALTITUDE
plt.plot(tout,altitude)
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')
plt.grid()
plt.show()

###VELOCITY
plt.figure()
plt.plot(tout,velout)
plt.xlabel('Time (sec)')
plt.ylabel('Total Speed (m/s)')
plt.grid()
plt.show()

##2D ORBIT
plt.figure()
plt.plot(xout,zout,'r-',label = 'Orbit')
plt.plot(xout[0],zout[0],'g*')
theta = np.linspace(0,2*np.pi,1000)
xplanet = Rplanet*np.sin(theta)
yplanet = Rplanet*np.cos(theta)
plt.plot(xplanet,yplanet,'b-',label = 'Planet')
plt.grid()
plt.legend()
plt.show()