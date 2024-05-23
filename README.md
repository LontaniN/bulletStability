# bulletStability
Codes for the analysis of stability of bullets and projectiles through the linearized theory that can be found in the "Modern Exterior Ballistics: The Launch and Flight Dynamics of Symmetric Projectiles" by Robert L. McCoy.

The bullet loaded is a .308 Sierra International Bullet tabulated in McCoy's book, with all the relevant aerodynamic coefficients.

## MAIN
Run this code to compute the basic version of the bullet stability theory. This means that the horizontal speed of the bullet and the roll rate are constants during the flight.

## MAIN_pseudosim
Run this code to compute a pseudo-simulation of the motion of the bullet druing flight. The code changes the speed of the bullet accordingly to the same linearized theory and also computes the aerodynamic coefficients at each iteration depending only on Mach number.
The changing roll rate during flight can be activated via the flag "flags.RollDamp".

### PLOTS
Containts the plots code

