classdef Lune < Shape
    
   properties (Access=public)
       radius;
   end
   methods 
       function c=Lune(r,n)% shild constructor
           t   = discretize(n);
           title='Lune';
           x   = r*cos(t) + 2*cos(2*t); y   = r*sin(t) + sin(2*t) + sin(3*t)/2;
           x_t =-r*sin(t) + 4*sin(2*t); y_t = r*cos(t) + 2*cos(2*t) + 3*cos(3*t)/2;
           x_tt=-r*cos(t) + 8*cos(2*t); y_tt=-r*sin(t) - 4*sin(2*t) - 9*sin(3*t)/2;
           
           ys   = -4*sin(t) + 7 * sin(2*t) - 6 * sin(3*t) + 2*sin(4*t);
           ys_t = -4*cos(t) + 14*cos(2*t) - 18*cos(3*t) + 8*cos(4*t);
           ys_tt=  4*sin(t) - 24*sin(2*t) + 54*sin(3*t) - 32*sin(4*t);
           
           y   = y   - ys/24;
           y_t = y_t - ys_t/24;
           y_tt= y_tt-ys_tt/24;
           
           c = c@Shape(title,n,t,x,x_t,x_tt,y,y_t,y_tt); % call superclass Shape constructor
           c.radius=r;
       end
       function graph(c)
           graph@Shape(c);
       end
           
   end
end
