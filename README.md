 Erik Hansen, 26.08.2020
 
 The following codes compute the pressure equilibrium when a flat rigid plane is loaded against an ideal elastic-plastic profile. The profile deformation is computed with the elastic half-space model. Subsurface plastic deformation is neglected, only plastic flow on the surface is taken into account. Outside of the contact zone, a constant hydrodynamic pressure can be prescribed. This pressure also serves as a lower limit for the pressure in the contact zone.
 The codes denoted with akchurin compute the pressure profile based on an imposed rigid body seperation h_s.
 The codes denoted with polonsky compute the pressure profile based on an imposed normal load force W_aim.
 The algorithm is based on the conjugate gardient method and fast Fourier transformations. 
 The codes denoted with linear use a linear convolution to determine the elastic deformation due to the pressure profile. Use this e.g. for single concentrated contacts.
 The codes denoted with circular use a circular convolution to determine the elastic deformation due to the pressure profile. This corresponds to periodic boundary conditions. Use this e.g. for periodic roughness profiles.
