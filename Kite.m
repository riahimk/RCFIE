classdef Kite < Shape
   properties (SetAccess=private,GetAccess=public)
       radius;
   end
   methods
       function c=Kite(r,n)
          t   = discretize(n);
          title='Kite';
          r=0.65;%% Kress paper ! 
          x   = cos(t) +0.65*cos(2*t) - 0.65; y   = 1.5*sin(t);
          x_t = -sin(t)-1.3*sin(2*t)  ;       y_t = 1.5*cos(t);
          x_tt= -cos(t)-2.6*cos(2*t)  ;       y_tt=-1.5*sin(t);
          x_ttt= sin(t)+5.2*sin(2*t)  ;       y_ttt=-1.5*cos(t);
          c = c@Shape(title,n,t,x,x_t,x_tt,x_ttt,y,y_t,y_tt,y_ttt); % call superclass Shape constructor
          c.radius=r;
       end
       function graph(c)
           graph@Shape(c);
       end
   end
end