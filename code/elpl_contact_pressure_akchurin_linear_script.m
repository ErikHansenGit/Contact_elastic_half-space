close all; clc; clear all;

%% Settings
% Algorithm settings:
h_s     = 4e-6;                                                             % [m] seperating distance between upper and lower reference planes
err_tol = 1e-10;                                                            % [-] relative error tolerance
h_ref   = 1e-6;                                                             % [m] reference length of relative error
it_max  = 1e2;                                                              % [-] maximum number of iterations

% Fluid settings:
p_hd = 1e5;                                                                 % [Pa] fluid pressure in non-contact zone

% Solid settings:
sld.nu_low      = 0.3;                                                      % [-] Poisson ratio of lower body
sld.E_low       = 210e9;                                                    % [Pa] Young's modulus of lower body
sld.nu_up       = 0.208;                                                    % [-] Poisson ratio of upper body
sld.E_up        = 81e9;                                                     % [Pa] Young's modulus of upper body
sld.E_dash      = 2/((1 - sld.nu_low^2)/sld.E_low + (1 - sld.nu_up^2)/sld.E_up); % [Pa] reduced modulus of elasticity
sld.H           = 5e8;                                                      % [Pa] hardness of the softer material as a maximum limit of the contact pressure

% Geometry settings: 
geo.Nx1         = 2^7;                                                      % [-] number of discretization points in the x1-direction
geo.Nx2         = 2^7;                                                      % [-] number of discretization points in the x2-direction
x1_Wb           = -0.25e-3;                                                 % [m] x1-coordinate of cell at the West boundary
x1_Eb           = 0.25e-3;                                                  % [m] x1-coordinate of cell at the East boundary
x2_Sb           = -0.25e-3;                                                 % [m] x2-coordinate of cell at the South boundary
x2_Nb           = 0.25e-3;                                                  % [m] x2-coordinate of cell at the North boundary
geo.Rx1         = 1.27e-2;                                                  % [m] ball radius

% Construct geometry:
geo.dx1         = (x1_Eb - x1_Wb)/(geo.Nx1 - 1);                            % [m] Spacing in x1-direction
geo.dx2         = (x2_Nb - x2_Sb)/(geo.Nx2 - 1);                            % [m] Spacing in x2-direction
geo.x1          = linspace(x1_Wb,x1_Eb,geo.Nx1);                            % [m] x1-coordinates with uniform discretization
geo.x2          = linspace(x2_Sb,x2_Nb,geo.Nx2);                            % [m] x2-coordinates with uniform discretization
[x1_matr, x2_matr] = ndgrid(geo.x1,geo.x2);
clear x1_Wb, clear x1_Eb; clear x2_Sb; clear x2_Nb;
% Ball-on-disc tribometer:
h_prof   = (x1_matr.^2)/(2*geo.Rx1) + (x2_matr.^2)/(2*geo.Rx1);             % [m] gap height variation induced by profile
h_prof   = h_prof - max(max(h_prof));
geo.prof = -h_prof;                                                         % [m] profile
clear x1_matr; clear x2_matr; clear h_prof;

% Add sinosoidal roughness:
a           = 1e-7;                                                         % [m] roughness amplitude
lambda_x1   = (geo.x1(geo.Nx1) - geo.x1(1) + geo.dx1)/4;                    % [m] roughness wavelengh in x1-direction
lambda_x2   = (geo.x2(geo.Nx2) - geo.x2(1) + geo.dx2)/7;                    % [m] roughness wavelengh in x2-direction
[x1, x2]    = ndgrid(geo.x1,geo.x2);                                        % Convert vectors to matrizes
rough_prof  = a*sin(x1/lambda_x1*2*pi).*cos(x2/lambda_x2*2*pi);             % [m] roughness profile
clear x1; clear x2; clear lambda_x1; clear lambda_x2; clear a;
% Superimpose roughness on profile:
geo.prof    = geo.prof + rough_prof;        
clear rough_prof;

%% Computation
% Compute 2-D fast Fourier transform of the elastic half-space Kernel
[fft2_Kernel]  = construct_linear_Kernel(geo.Nx1,geo.Nx2,geo.dx1,geo.dx2,sld.E_dash);

% Compute pressure field:
p_con_ini = zeros(geo.Nx1,geo.Nx2);                                         % [Pa] first guess of contact pressure field
[p_con,g,err] = elpl_contact_pressure_akchurin_linear(p_hd,sld.H,p_con_ini,geo.prof,h_s,geo.Nx1,geo.Nx2,fft2_Kernel,err_tol,it_max,h_ref);
clear p_con_ini;
% Compute resulting normal load:
W = geo.dx1*geo.dx2*sum(p_con(:));

% Compute deformations and deformed profiles:
[h_el] = compute_h_el(p_con,geo.Nx1,geo.Nx2,fft2_Kernel);                   % [m] elastic deformation
clear fft2_Kernel;
z_el = geo.prof - h_el;                                                     % [m] elastically deformed profile
A_pl = find(p_con>=sld.H);                                                  % [-] indices of plastically deformed points
h_pl = zeros(geo.Nx1,geo.Nx2);
h_pl(A_pl) = z_el(A_pl) - h_s;                                              % [m] plastic deformation                          
clear A_pl;
z_elpl = z_el - h_pl;                                                       % [m] elastically and plastically deformed profile

%% Plot
% Plot seetings:
% Line plots:
KIT_colorlist={[0,150,130]/255,[162 34 35]/255,[70 100 170]/255,[252 229 0]/255,[140 182 60]/256,[223 155 27]/255,[167 130 46]/255,[163 16 124]/255,[35 161 224]/255};
widthlines          = 1.2;
% Surface plots:
colmap = parula;

% Plots:
figure('Units','centimeters','Position',[01 13 7 7])
surf(geo.x1,geo.x2,geo.prof')
xlabel('x_1 [m]')
ylabel('x_2 [m]')
zlabel('z [m]')
title('Undeformed profile')
shading interp
material dull
colormap(colmap)
camlight

figure('Units','centimeters','Position',[09 13 7 7])
surf(geo.x1,geo.x2,h_el')
xlabel('x_1 [m]')
ylabel('x_2 [m]')
zlabel('h_{el} [m]')
title('Elastic deformation')
shading interp
material dull
colormap(colmap)
camlight

figure('Units','centimeters','Position',[17 13 7 7])
surf(geo.x1,geo.x2,h_pl')
xlabel('x_1 [m]')
ylabel('x_2 [m]')
zlabel('h_{pl} [m]')
title('Plastic deformation')
shading interp
material dull
colormap(colmap)
camlight

figure('Units','centimeters','Position',[25 13 7 7])
surf(geo.x1,geo.x2,z_elpl')
xlabel('x_1 [m]')
ylabel('x_2 [m]')
zlabel('z_{elpl} [m]')
title('Deformed profile')
shading interp
material dull
colormap(colmap)
camlight

figure('Units','centimeters','Position',[01 00 7 7])
semilogy(err,'LineWidth',widthlines,'Color',KIT_colorlist{1})
xlabel('i_{it} [m]')
ylabel('err [-]')
title('Relative error')

figure('Units','centimeters','Position',[11 00 7 7])
surf(geo.x1,geo.x2,p_con')
xlabel('x_1 [m]')
ylabel('x_2 [m]')
zlabel('p_{con} [m]')
title('Contact pressure')
shading interp
material dull
colormap(colmap)
camlight

figure('Units','centimeters','Position',[21 00 7 7])
surf(geo.x1,geo.x2,g')
xlabel('x_1 [m]')
ylabel('x_2 [m]')
zlabel('g [m]')
title({'Residual'; '{\itAim in hydrodynamic region: g<=0}'; '{\itAim in elastic region: g=0}'; '{\itAim in plastic region: g>=0}'})
shading interp
material dull
colormap(colmap)
camlight

clear KIT_colorlist; clear colmap; clear widthlines;

function [p_con,g,err] = elpl_contact_pressure_akchurin_linear(p_min,H,p_con_ini,z,h_s,Nx1,Nx2,fft2_Kernel,err_tol,it_max,h_ref)
% Erik Hansen, 26.08.2020
% Calculates the contact pressure occuring when a rigid smooth surface is
% loaded against an elastic half-space. The minimum and maximum contact pressures can be set to limiting values. Subsurface plastic flow within the half-space is not modelled.
% It is assumed that plastic flow only occures on the the half-space's surface if the maximum contact pressure is reached.
% The profile deformation is computed with the elastic half-space model using the Boundary Element Method (BEM).
% The code is very similar to and consists largely of the 
% algorithm described Akchurin et al., 2015:
% "On a model for the prediction of the friction coefficient in mixed lubrication based on a load-sharing concept with measured surface roughness"
% However, some changes were introduced in the treatment of non-contact and
% plastic contact points. Inspiration was taken from Vollebregt, 2014:
% "A new solver for the elastic normal contact problem using conjugate gradients, deflation, and an FFT-based preconditioner"
% The Matlab implementation of the linear convolution in the fourier space
% was largely inspired by Sainsot and Lubrecht, 2011:
% "Efficient solution of the dry contact of rough surfaces: a comparison of fast Fourier transform and multigrid methods"
% Input:
% p_min             [Pa]    pressure in non contact zones as a lower limit for the contact pressure
% H                 [Pa]    hardness of the softer material as a maximum limit of the contact pressure
% p_con_ini         [Pa]    initial guess of contact pressure field
% z                 [m]     undeformed profile
% h_s               [m]     seperating distance between upper and lower reference planes
% Nx1               [-]     number of points in x1-direction
% Nx2               [-]     number pf points in x2-direction
% fft2_Kernel       [m/Pa]  2-D fast Fourier transform of the Kernel
% err_tol           [-]     relative error tolerance
% it_max            [-]     maximum number of iterations
% h_ref             [m]     reference length of relative error
% Output:
% p_con             [Pa]    contact pressure field
% g                 [m]     residual of the gap height distribution
% err               [-]     relative error
% -------------------------------------------------------------------------
p_con = p_con_ini;                              % [Pa] pressure field
[u] = compute_h_el(p_con,Nx1,Nx2,fft2_Kernel);  % [m] elastic deformation 
g = -u + z - h_s;                               % [m] residual of the gap height distribution

% Find indices of points that are in the non-contact, elastic and plastic
% domain due to the pressure distribution. Non-contact and plastic points 
% are also evaluated whether they are correctly or not in surface contact
% due to the gap height distribution
A_el    = find(p_con>p_min  &p_con<H);
A_nc_cr = find(p_con<=p_min &g<=0);
A_nc_wr = find(p_con<=p_min &g>0);
A_pl_cr = find(p_con>=H     &g>=0);
A_pl_wr = find(p_con>=H     &g<0);
% Within these points, the pressure distribution needs to be adjusted:
A_free = union(union(A_el,A_nc_wr),A_pl_wr);
 
G = sum(g(A_free).*g(A_free));  % [m^2] norm of the residual
G_old = 1;                      % [m^2] previous norm of the residual
delta = 0;                      % [-] flag whether to use conjugate gradient or steepest descend
err = zeros(it_max,1);          % [-] relative error
i_it=0;                         % [-] iteration counter        
t = zeros(Nx1,Nx2);             % [Pa] search direction
while i_it == 0 || (err(i_it)>err_tol && i_it<it_max)
    i_it = i_it + 1;    
    
    % Find search direction:
    t(A_free)   = g(A_free) + delta*(G/G_old)*t(A_free);
    t(A_nc_cr)  = 0;
    t(A_pl_cr)  = 0;
    clear delta;
    
    % Determine step length tau:
    [t_conv]    = compute_h_el(t,Nx1,Nx2,fft2_Kernel);     
    r           = -t_conv;
    clear t_conv;
    tau = (sum(g(A_free).*t(A_free)))/(sum(r(A_free).*t(A_free)));
    clear r;

    % Update pressure:
    p_con(A_free)       = p_con(A_free) - tau*t(A_free);
    p_con(p_con<p_min)  = p_min;
    p_con(p_con>H)      = H;
    
    % Compute elsatic deformation to find new residual of the gap height distribution
    [u] = compute_h_el(p_con,Nx1,Nx2,fft2_Kernel); 
    g = -u + z - h_s;
    clear u;
    
    % Find indices of points that are in the non-contact, elastic and plastic domain:
    A_el    = find(p_con>p_min  &p_con<H);
    A_nc_cr = find(p_con<=p_min &g<=0);
    A_nc_wr = find(p_con<=p_min &g>0);
    A_pl_cr = find(p_con>=H     &g>=0);
    A_pl_wr = find(p_con>=H     &g<0);
    % Within these points, the pressure distribution needs to be adjusted:
    A_free = union(union(A_el,A_nc_wr),A_pl_wr);
    
    % Determine whether to use conjugate gradient or steepest descend in the next ieration
    if isempty(A_nc_wr) && isempty(A_pl_wr)
    	delta = 1;
    else
        delta = 0;
    end

    % Save G for the next iteration:
    G_old = G;
    % Compute norm of the residual:
    G = sum(g(A_free).*g(A_free));
    % Compute relative error
    err(i_it) = sqrt(G)/h_ref;
end
% Resize err:
err = err(1:i_it,1);
end


function [h_el] = compute_h_el(p,Nx1,Nx2,fft2_Kernel)
% Computes elastic displacement due to pressure field with linear
% convolution in Fourier space
% Input:
% p                 [Pa]        pressure field
% Nx1               [-]         number of discretized points in x1-direction
% Nx2               [-]       	number of discretized points in x2-direction
% fft2_Kernel       [?]         2-D fast Fourier transform of the Kernel
% Output:
% h_el              [m]         elastic deformation
% -------------------------------------------------------------------------
% Extend pressure field for linear convolution:
p_ext               = zeros(2*Nx1,2*Nx2);
p_ext(1:Nx1,1:Nx2)  = p(:,:);
% Compute convolution in Fourier space for better pefromance:
h_el_ext            = real(ifft2(fft2_Kernel.*fft2(p_ext)));
% Extract linear convolution:
h_el                = h_el_ext(1:Nx1,1:Nx2);
end

function [fft2_Kernel]  = construct_linear_Kernel(Nx1,Nx2,dx1,dx2,E_dash)
% Calculates the Kernel function for a linear convolution in the
% influence area of size Nx1*Nx2
% The Kernel is constructed for an imposed normal load pressure as explained by 
% Johnson, K. L., 2004. Contact mechanics. Cambridge: Cambridge University Press
% in equation (3.25) on P.54
% and the definition of 
% sld.E_dash = 2/((1 - sld.nu_low^2)/sld.E_low + (1 - sld.nu_up^2)/sld.E_up).
% The Kernel center is at the edge of the domain
% The Kernel domain is of size 2*Nx1*2*Nx2
% Input:
% Nx1                 [-]       Number of discretized points in x1-direction
% Nx2                 [-]       Number of discretized points in x2-direction
% dx1                 [m]       Length of discretized cell in x1-direction
% dx2                 [m]       Length of discretized cell in x2-direction
% E_dash              [Pa]      reduced modulus of elasticity
% Output:
% fft2_Kernel         [?]       2-D fast Fourier transform of the Kernel
% -------------------------------------------------------------------------
Nx1_K = 2*Nx1;      % [-]       Number of discretized Kernel points in x1-direction
Nx2_K = 2*Nx2;      % [-]       Number of discretized Kernel points in x2-direction
dx1_mod = dx1/2;
dx2_mod = dx2/2;
% Determine distances:
i           = 1:Nx1_K;
j           = 1:Nx2_K;
i_cond      = (i <= floor(Nx1_K/2));
x1          = -(((floor(Nx1_K/2) + 1) - (i - (ceil(Nx1_K/2) + 1))) - 1)*dx1;
x1(i_cond)  = (i(i_cond) - 1)*dx1;
j_cond      = (j <= floor(Nx2_K/2));
x2          = -(((floor(Nx2_K/2) + 1) - (j - (ceil(Nx2_K/2) + 1))) - 1)*dx2;
x2(j_cond)  = (j(j_cond) - 1)*dx2;
clear i; clear j; clear i_cond; clear j_cond;
clear Nx1; clear Nx2; clear dx1; clear dx2;
[x1, x2]    = ndgrid(x1,x2); % Convert vectors to matrizes
% Construct Kernel:
term_1 = (x1 + dx1_mod).*log(ext_sqrt(x2 + dx2_mod, x1 + dx1_mod)./ext_sqrt(x2 - dx2_mod, x1 + dx1_mod));
term_2 = (x2 + dx2_mod).*log(ext_sqrt(x1 + dx1_mod, x2 + dx2_mod)./ext_sqrt(x1 - dx1_mod, x2 + dx2_mod));
term_3 = (x1 - dx1_mod).*log(ext_sqrt(x2 - dx2_mod, x1 - dx1_mod)./ext_sqrt(x2 + dx2_mod, x1 - dx1_mod));
term_4 = (x2 - dx2_mod).*log(ext_sqrt(x1 - dx1_mod, x2 - dx2_mod)./ext_sqrt(x1 + dx1_mod, x2 - dx2_mod));
clear x1;  clear x2; clear dx1_mod;  clear dx2_mod; 
Kernel        = 2/(pi*E_dash)*(term_1 + term_2 + term_3 + term_4);
clear term_1; clear term_2; clear term_3; clear term_4; 
fft2_Kernel   = fft2(Kernel);
function return_value = ext_sqrt(p,q)
% Auxiliary function for Kernel construction
return_value = p + sqrt(p.^2 + q.^2);
end
end