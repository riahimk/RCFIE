classdef Shape <  handle
    properties (SetAccess=private,GetAccess=public)
        title;
        radius;
        gon;
         n;              % number of discretizations points pi/n
         x;  y;          % cordinate of any point in the boundary X=(x,y)
        x_t; y_t;        % current celerity
        x_tt;y_tt;       % current curvature
        x_ttt;y_ttt;
        arclength;       % Arclength of the shape
        ARCL;
        normal;          % outward normal on the boundary perponducular to the tnagent
        tangent;         % tangent vector perpenducular to the normal 
        distance;        % matrix M_{i,j} := ||X(t_i)-X(t_j)|| -> X= [x;y];
        factor ;  %x_t_dot_x_tau_ov_arcl % Matrix for neumann pre - calculus 
    end
    methods (Access=public)
        function s=Shape(title,n,t,x,x_t,x_tt,x_ttt,y,y_t,y_tt,y_ttt) % class constructor
            s.n=n;
            s.x    = x;         
            s.y    = y;
            s.x_t  = x_t;       s.y_t  =y_t;  
            s.x_tt = x_tt;      s.y_tt =y_tt;
            s.x_ttt= x_ttt;     s.y_ttt=y_ttt;
            s.arclength = sqrt(s.x_t.^2 + s.y_t.^2);
            s.ARCL       = s.arclength'*s.arclength;
            for i=1:(2*n)
                s.tangent(:,i) = [s.x_t(i) ; s.y_t(i)]; 
                s.normal (:,i) = [s.y_t(i) ;-s.x_t(i)];
                for j=1:(2*n)
                    s.distance(i,j) = sqrt( ( x(i) - x(j) )^2 + ( y(i)-y(j) )^2 );
                    s.factor  (i,j) =  (x_t(i)*x_t(j) + y_t(i)*y_t(j))/s.arclength(j); 
                end
            end
        end
        function graph(s)
            figure(); hold on; plot(s.x,s.y);hold off;
        end
    end
end