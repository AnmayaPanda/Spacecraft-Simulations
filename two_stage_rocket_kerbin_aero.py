'''KERBIN AERODYNAMICS MODEL IN A TWO STAGE ROCKET'''
###Import all the modules
import numpy as np ###numeric python
import matplotlib.pyplot as plt ###matlab style plotting
import scipy.integrate as sci ##integration toolbox

plt.close('all')

##DEFINE SOME CONSTANT PARAMETERS
G = 6.6742*10**-11 ##Gravitational Constant (SI Units)

##PLANET
##EARTH
#Rplanet = 6371000 ##meters
#mplanet = 5.972e24 #kg
#name = 'Earth'

##KERBIN
Rplanet = 600000 #meters
mplanet = 5.2915158*10**22 #kg
name = 'Kerbin'

##PARAMETERS OF ROCKET
###Initial Conditions - For Single Stage Rocket
x0 = Rplanet
z0 = 0.0
r0 = 200000 + Rplanet
velx0 = 0.0
velz0 = 0.0
period = 12000 #2*np.pi/np.sqrt(G*mplanet)*r0**(3.0/2.0)*1.5 ##multiplication factor to move it in elliptic orbit
weighttons = 5.3
mass0 = weighttons*2000/2.2 ##kg ---> Initial mass
max_thrust = 167970.0 ##Newtons
Isp1 = 250.0 #seconds
Isp2 = 400.0 #seconds
tMECO = 38.0 #seconds ------> Main Engine CutOff time
tsep1 = 2.0 #seconds -------> Length of time to remove first stage
mass1tons = 0.2
mass1 = mass1tons*2000/2.2 ##kg
t2start = 264.0 #seconds ----> Second stage burn time
t2end = t2start + 20.0 #seconds --> we keep it 20 sec as we still have mass to burn, keeping 17.5 crashes it as apogee comes closer due to smaller orbit
D = 0.85 ## Diameter in meters
# S = np.pi*(D/2.0)**2 #Cross sectional area of the rocket in sq. meters
CD = 0.1 ##Drag coefficient

##Create a Aerodynamics class
class Aerodynamics():
    def __init__(self,name):
        self.name = name
        if name == 'Kerbin':
        ##Import the Aero model for Interpolation
            data = np.loadtxt('kerbin_aerodynamics.txt',delimiter=',')
            #print(data)
            self.altitude = data[:,0]
            #print(self.altitude)
            self.density = data[:,3]
            print(self.density)
            self.rhos = self.density[0]
            self.beta = 0.0
        elif name == 'Earth':
        ##We use the Earth Aero model
            self.beta = 0.1354/1000.0 #density constant
            self.rhos = 1.225 #kg/m^3

    def getDensity(self,altitude):
        if self.name == 'Kerbin':
            ##interpolate
            rho = np.interp(altitude,self.altitude,self.density)
        elif self.name == 'Earth':
            rho = self.rhos*np.exp(-self.beta*altitude)
        return rho
        
    
##Create a aeroModel variable which is an instance of the class Aerodynamics and putting it up here to make it global
aeroModel = Aerodynamics(name)


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
    return np.asarray([accelx,accelz]),r
'''
## Considering open loop
def thrust(t):
    global max_thrust,Isp
    if t<5:
        thrustF = max_thrust
    else:
        thrustF = 0.0
    ##Angle of thrust
    theta = 10*np.pi/180.0
    thrustx = thrustF*np.cos(theta)
    thrustz = thrustF*np.sin(theta)
    #mdot
    #ve = Isp*9.81 #Exit velocity in m/s
    #mdot = -thrustF/ve 
    return np.asarray([thrustx,thrustz])
'''
## It is no longer thrust as we also considering change in mass, hence, propulsion

def propulsion(t): 
    global max_thrust,Isp1,Isp2,tMECO
    ## Timing for thrusters
    if t < tMECO:
    ## We are firing the main thruster
        theta = 10*np.pi/180.0 
        thrustF = max_thrust
        ve = Isp1*9.81 #Exit velocity in m/s
        mdot = -thrustF/ve
    if t > tMECO and t < (tMECO + tsep1): ## tsep1 is Time taken to separate the first stage
        theta = 0.0
        thrustF = 0.0
        ##masslost = mass1
        mdot = -mass1/tsep1
    if t > (tMECO + tsep1):
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0
    if t > t2start and t < t2end:
        theta = 90.0*np.pi/180.0
        thrustF = max_thrust
        ve = Isp2*9.81 #Exit velocity in m/s
        mdot = -thrustF/ve

    if t > t2end:
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0
    
    ##Angle of thrust
    thrustx = thrustF*np.cos(theta)
    thrustz = thrustF*np.sin(theta)
    return np.asarray([thrustx,thrustz]), mdot

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
    global aeroModel, Rplanet
    #state vector 
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    mass = state[4]

    #Compute xdot and zdot - Kinematic Relationship
    xdot = velx
    zdot = velz

    ###Compute the total Forces
    ###Gravity
    accel,r = gravity(x,z) ## r = distance from the center
    gravityF = -accel*mass

    ###Aerodynamics
    altitude = r - Rplanet ## Altitude above the surface
    rho = aeroModel.getDensity(altitude) ## Air Density
    V = np.sqrt(velx**2+velz**2) ## Total velocity
    qinf = (np.pi/8.0)*rho*(D**2)*abs(V) ##Dynamic Pressure
    aeroF = -qinf*CD*np.asarray([velx,velz])

    ###Thrust
    thrustF,mdot = propulsion(t)

    Forces = gravityF+aeroF+thrustF

    #Compute Acceleration
    if mass > 0:
        ddot = Forces/mass
    else:
        ddot = 0.0
        mdot = 0.0

    #Compute the statedot
    statedot = np.asarray([xdot,zdot,ddot[0],ddot[1],mdot])

    return statedot

#############EVERYTHING BELOW HERE IS THE MAIN SCRIPT###################

###Test Surface Gravity
print('Surface Gravity (m/s^2)',gravity(0,Rplanet))

###Plot the Air Density as a function of AGL (Above Ground Level)
test_altitude = np.linspace(0,100000,100)
test_rho = aeroModel.getDensity(test_altitude)
plt.figure()
plt.subplot(2,3,1)
plt.plot(test_altitude,test_rho,'b-')
#plt.title('Figure 1 - Air Density as a Function of Above Ground Level')
plt.xlabel('Altitude (m)')
plt.ylabel('Air Density (kg/m^3)')
plt.grid()
#plt.show()

###Initial Conditions - For Orbit
'''
x0 = Rplanet+600000 ##m
#To start at surface of earth x0 = Rplanet with some velz0 or vertical sppeed, say 100.0 m/s, it sort of crashes
z0 = 0.0
r0 = np.sqrt(x0**2+z0**2)
velx0 = 0.0
velz0 = np.sqrt(G*mplanet/r0)*1.1 ##m/s  multiplication factor to make it in elliptic orbit
stateinitial = np.asarray([x0,z0,velx0,velz0])
period = 2*np.pi/np.sqrt(G*mplanet)*r0**(3.0/2.0)*1.5 ##multiplication factor to move it in elliptic orbit
'''
#Populate Initial Condition Vector
stateinitial = np.asarray([x0,z0,velx0,velz0,mass0])
##Time Window
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
massout = stateout[:,4]

##Plot

###ALTITUDE
plt.subplot(2,3,2)
plt.plot(tout,altitude)
#plt.title('Figure 2 - Altitude as a Function of Time')
plt.xlabel('Time (sec)')
plt.ylabel('Altitude (m)')
plt.grid()
#plt.show()

###VELOCITY
plt.subplot(2,3,3)
plt.plot(tout,velout)
#plt.title('Figure 3 - Velocity as a Function of Time')
plt.xlabel('Time (sec)')
plt.ylabel('Total Speed (m/s)')
plt.grid()
#plt.show()

##MASS
plt.subplot(2,3,4)
plt.plot(tout,massout)
#plt.title('Figure 4 - Mass as a Function of Time')
plt.xlabel('Time (sec)')
plt.ylabel('Mass (kg)')
plt.grid()
#plt.show()

##2D ORBIT
plt.subplot(2,3,5)
plt.plot(xout,zout,'r-',label = 'Orbit')
plt.plot(xout[0],zout[0],'g*')
theta = np.linspace(0,2*np.pi,1000)
xplanet = Rplanet*np.sin(theta)
yplanet = Rplanet*np.cos(theta)
plt.plot(xplanet,yplanet,'b-',label = 'Planet')
#plt.title('Figure 5 - 2D Orbit')
plt.grid()
plt.legend()
plt.show()
