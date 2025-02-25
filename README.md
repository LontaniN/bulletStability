# bulletStability
Codes for the analysis of stability of bullets and projectiles through the linearized theory that can be found in the "Modern Exterior Ballistics: The Launch and Flight Dynamics of Symmetric Projectiles" by Robert L. McCoy.

The bullet loaded is a .308 Sierra International Bullet tabulated in McCoy's book, with all the relevant aerodynamic coefficients.
The file "Report.pdf" contains a comprehensive explanation of the methods applied and of the results obtained. It was written for the course of Fundamental of Ballistics (2024) at Politecnico di Milano.
## MAIN
Run this code to compute the epicyclic trajectory of the bullet. It loads the data of the chosen bullet from the DATA folder and checks if the aerodynamic coefficients loaded are arrays, which means that they depend on Mach number. If they are it computes the trajectory taking into account the linear deceleration due to drag and possibily also the changing roll rate, if **"flags.RollDamp"** is set to true. Otherwise it computes the simplified trajectory that assumes constant speed and constat roll rate. The dependence of the aerodynamic coefficients from the angle of attack is not taken into account as the theory behind the code is the "linearized pitching and yawing motion of projectiles" from **"Modern Exterior Ballistics: The Launch and Flight Dynamics of Symmetric Projectiles"** by Robert L. McCoy.

### The parameters and initial conditions that need to be set for computations are:
* sMax_metres, maximum downrange in metres
* twistRateInch, twist rate of spin in the rifle barrel in inches per turn
* deltaMax, maximum yaw at muzzle (usually measured with cameras), can be modified at will to change the magnitude of the perturbation
* Enviroment data: air temperature, air density and gravitational acceleration constant
### PLOTS
Containts the plots code, including the animation one.


![immagine](https://github.com/LontaniN/bulletStability/assets/93401408/c2c905f9-65ce-48ab-bcf3-e813e102bb97)

![immagine](https://github.com/LontaniN/bulletStability/assets/93401408/644f9364-7cf3-4709-a968-465529a1206d)


