%% General description of the program
%
% Numerical Algorithm for the resolution of the scattering problem governed
% by the following Helmholtz equation
%
% $$ \Delta u^s + k^2 u^s  = 0 \quad in\quad R\backslash \overline\Omega$$
%
% with the boundary condition:
%
% $$ \partial_n u^s - Z u^s = -\partial_n u^{i} - Z u^{i}\quad on\quad \Gamma$$
%
% Complemented with the radiation condition at $\infty$
%
% $$\lim_{|r|\rightarrow\infty}r^{\frac{1}{2}} (\partial_r u^s - ik u^s) =0,$$
%clf();clc();
!rm tabular*
%!rm ../../NOTES/note_1/tabular*.tex
clc;clear;
for tpe=[1,2,3]%,4]
    %disp '-------- ~ Scatt problem type ~--------'
    for wn=[1,5]
        for n=[8,16,32,64,128]%,256,512,1024]
            for regul=0:1
                %disp '-------- ~ Reg ~--------'
                
                %tpe=1;regul=1;n=8
                %% parameters of the obstacle scattering problem
                %digits(16);
                myScattering_problem.TYPE     =['Dirichlet ';'Neumann   ';'Impedance ';'GImpedance'];% padded
                myScattering_problem.d        =[1; 0];
                myScattering_problem.k        = wn;    % k: wave number of the scattering problem
                myScattering_problem.ik       = 1i*wn;    % k: wave number of the scattering problem
                myScattering_problem.eta      = wn;    % coupling parameter
                myScattering_problem.zed      = wn;  % Impedance multiplicatif coefficient
                myScattering_problem.T        = 2*pi/myScattering_problem.k;
                myScattering_problem.n        = n;  % \pi/n discretization points on $\Gamma$
                myScattering_problem.regu     = regul;
                myScattering_problem.type     = myScattering_problem.TYPE(tpe,:);
                
                fprintf('\n %s Scatt-Pb,k=%d, n = %d ,Reg=%d. '...
                    ,myScattering_problem.type ...
                    ,myScattering_problem.k ...
                    ,myScattering_problem.n    ...
                    ,myScattering_problem.regu);
                %% instantiation of the scatterer's shape Kite
                myShape=Kite(.65,myScattering_problem.n);      %  Kite shaped with rotation -90°
                %myShape=Circle(1,myScattering_problem.n);      %  Kite shaped with rotation -90°
                %% Plot the graph-options of myShape Kite
                %graph(myShape);     % plot the shape as parametric two-dimension function
                
                
                %% Solve system $A\Phi=g$ with the iterative solver e.g. GMRES
                % Build up the discret linear system $A\Phi=g$ that corresponds ...
                % to the combined acoustic double- and single-layer potential
                %
                % $$u^s := DL - i\eta SL$$
                %
                t   = discretize(myScattering_problem.n);
                %g   = Plane_wave(myShape,myScattering_problem);plot3(myShape.x,myShape.y,g);
                %U   = Exact2D(myScattering_problem,g,100);
                
                %tic
                lin_sys  = Problem(myShape,myScattering_problem);
                
                Sol      = Solve  (myShape,lin_sys,myScattering_problem);
                
                %Sol.density(1)
                fprintf(' far field (R=%0d) -->%f,%f\n',regul,real(Sol.far_field(1)),imag(Sol.far_field(1)));
                
                TAB(regul+1)=Sol.far_field(1);
               
            end;
            CREATE_TABULAR(myScattering_problem,TAB);
            TAB=[];
        end;
    end;
end;
%-----------------------------------
!cp tabular*.tex ../../NOTES/note_1/


%toc
%-----------------------------------
%figure(1);hold on;
%plot3(myShape.x,myShape.y,real(Sol.density));
%savefig([myShape.title 'density.fig']);
%hold off;
%-----------------------------------
%disp '-------------------------';
%disp '  ** Succesfull end **   ';
%disp '-------------------------';


%% Authors - Copyright :
% Mohamed Kamel RIAHI : riahi@njit.edu
% Joint work with Yassine BOUBENDIR and Catalin TURC
%
% *Copyright*: New Jersey Institute of Technology
%
% Last modified Nov 17 2014