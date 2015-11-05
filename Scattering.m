%% General description of the program
%
% Numerical Algorithm for the resolution of the scattering problem governed
% by the following Helmholtz equation 
%
% $$ \Delta u^s + k^2 u^s  = 0 \quad in\quad R\backslash \overline\Omega$$
% 
% with the boundary condition:
%
% $$ \partial_n u^s + Z u^s = -\partial_n u^{i} - Z u^{i}\quad on\quad \Gamma$$
% 
% Complemented with the radiation condition at $\infty$
%
% $$\lim_{|r|\rightarrow\infty}r^{\frac{1}{2}} (\partial_r u^s - ik u^s) =0,$$



%% parameters of the obstacle scattering problem
myScattering_problem.TYPE     =['Dirichlet ';'Neumann   ';'Impedance ';'GImpedance'];% padded
myScattering_problem.n        = 20;   % \pi/n discretization points on the boundary $\Gamma$
myScattering_problem.eta      =.2;    % coupling parameter
myScattering_problem.k        =.4;    % k: wave number of the scattering problem 
myScattering_problem.d        =[1; 0];
myScattering_problem.kk       = myScattering_problem.k^2;
myScattering_problem.type     = myScattering_problem.TYPE(1,:);
                               
%% instantiation of the scatterer's shape unite Circle
%myShape=Circle(1,myScattering_problem.n);      %  circle of radius 1 composed by 2*n discret points
%% Plot the graph-options of myShape
%graph(myShape);     % plot the shape as parametric two-dimension function


%% instantiation of the scatterer's shape Ellipse
%myShape=Ellipse(1,.2,myScattering_problem.n);      %  circle of radius 1 composed by 2*n discret points
%% Plot the graph-options of myShape Ellipse
%graph(myShape);     % plot the shape as parametric two-dimension function


%% instantiation of the scatterer's shape Kite C-shaped
%myShape=Lune(1,myScattering_problem.n);      %  Kite C-shaped with rotation -90°
%% Plot the graph-options of myShape Kite C-shaped
%graph(myShape);     % plot the shape as parametric two-dimension function


%% instantiation of the scatterer's shape Kite  
myShape=Kite(1,myScattering_problem.n);      %  Kite shaped with rotation -90°
%% Plot the graph-options of myShape Kite
%graph(myShape);     % plot the shape as parametric two-dimension function


%% Properties (read-only) of myShape are:

%myShape


%% Instantiation class "single-layer" potential integral operator
%tic
%disp 'Single-layer operator class construction'
%SL=Single_Layer(myShape,myScattering_problem)
%toc


%% Instantiation class "double-layer" potential integral operator
%tic 
%disp 'Double-layer operator class construction'
%DL=Double_Layer(myShape,myScattering_problem)
%toc


%% Instantiation of the right hand side vector $g$ for the linear system $A\Phi=g$




%% Solve system $A\Phi=g$ with the iterative solver GMRES
% Build up the discret linear system $A\Phi=g$ that corresponds ...
% to the combined acoustic double- and single-layer potential
%
% $$u^s := DL -i\eta SL$$ 
% 
tic
Phi=Solve(myShape,myScattering_problem);
toc 
disp ' ** Succesfull end **'


%% Authors - Copyright : 
% Mohamed Kamel RIAHI : riahi@njit.edu
% Joint work with Yassine BOUBENDIR and Catalin TURC
%
% *Copyright*: New Jersey Institute of Technology
%
% Last modified Nov 17 2014