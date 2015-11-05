classdef Polygon < Shape
    methods
       function c=Polygon(gon,n)
          t   = discretize(n);
          title='Polygone';
          x=2*cos(pi/2/gon)*cos(.5*(t+pi/gon*(gon*t/pi+1)))-cos(pi/gon*(gon*t/pi+1));
          y=2*cos(pi/2/gon)*sin(.5*(t+pi/gon*(gon*t/pi+1)))-sin(pi/gon*(gon*t/pi+1));
          
          x_t  = x  ;       y_t = y ;
          x_tt = x  ;       y_tt= y ;
          x_ttt= x  ;       y_ttt=y ;
          c = c@Shape(title,n,t,x,x_t,x_tt,x_ttt,y,y_t,y_tt,y_ttt); % call superclass Shape constructor
          %c.gon=gon;
       end
       function graph(c)
           graph@Shape(c);
       end        
    end
    
end