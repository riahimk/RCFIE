classdef Ellipse < Shape
   properties (SetAccess=private,GetAccess=public)
       radius;
   end
   methods
       function c=Ellipse(r1,r2,n)           
           t =discretize(n);
           title='Ellipse';
           x =r1*cos(t);  x_t=-r1*sin(t); x_tt =-x;% -r*cos(t);
           y =r2*sin(t);  y_t= r2*cos(t); y_tt =-y;% -r*sin(t); 
           c = c@Shape(title,n,t,x,x_t,x_tt,y,y_t,y_tt); % call superclass Shape constructor
           c.radius=[r1,r2];
           
       end
       function graph(c)
           graph@Shape(c);
       end 
   end
end