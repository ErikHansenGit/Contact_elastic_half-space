The following MATLAB codes compute the pressure-deformation equilibrium when a flat rigid plane is loaded against an ideal elastic-plastic profile. The profile deformation is computed with the elastic half-space model using the Boundary Element Method (BEM). Subsurface plastic deformation is neglected, only plastic flow on the surface is taken into account to limit the maximum contact pressure. Outside of the contact zone, a constant hydrodynamic pressure can be prescribed. This pressure also serves as a lower limit for the pressure within the contact zone.
 
 The scripts consist of an exemplary case setup, the contact solver functions and plot commands. By executing the script, the exemplary contact problem is solved by calling the embedded solver and the results are visualized in plots.
 
 The codes denoted with akchurin compute the pressure profile based on an imposed rigid body seperation h_s.
 
 The codes denoted with polonsky compute the pressure profile based on an imposed normal load force W_aim.
 
 
 The algorithm is based on the conjugate gradient method and fast Fourier transformations. 
 
 The codes denoted with linear use a linear convolution to determine the elastic deformation due to the pressure profile. Use this e.g. for single concentrated contacts.
 
 The codes denoted with circular use a circular convolution to determine the elastic deformation due to the pressure profile. This corresponds to periodic boundary conditions. Use this e.g. for periodic roughness profiles.

This code  is free to use by anyone. Still, citing this repository is greatly appreciated if it is of use to your work.

Erik Hansen, 07.09.2020
