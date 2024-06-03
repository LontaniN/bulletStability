# bulletStability
Codes for the analysis of stability of bullets and projectiles through the linearized theory that can be found in the "Modern Exterior Ballistics: The Launch and Flight Dynamics of Symmetric Projectiles" by Robert L. McCoy.

The bullet loaded is a .308 Sierra International Bullet tabulated in McCoy's book, with all the relevant aerodynamic coefficients.

## MAIN
Run this code to compute the epycyclic trajectory of the bullet. It loads the data of the chosen bullet from the DATA folder and checks if the aerodynamic coefficients loaded are arrays, which means that they depend on Mach number. If they are it computes the trajectory taking into account the linear deceleration due to drag and possibily also the changing roll rate, if **"flags.RollDamp"** is set to true. Otherwise it computes the simplified trajectory that asumes constant speed and constat roll rate. The dependence of the aerodynamic coefficients from the angle of attack is not taken into account as the theory behind the code is the "linearized pitching and yawing motion of projectiles" from **"Modern Exterior Ballistics: The Launch and Flight Dynamics of Symmetric Projectiles"** by Robert L. McCoy.

### The parameters and initial conditions that need to be set for computations are:
* sMax_metres, maximum downrange in metres
* twistRateInch, twist rate of spin in the rifle barrel in inches per turn
* deltaMax, maximum yaw at muzzle (usually measured with cameras), can be modified at will to change the magnitude of the perturbation
* Enviroment data: air temperature, air density and gravitational acceleration constant
### PLOTS
Containts the plots code, including the animation one.

