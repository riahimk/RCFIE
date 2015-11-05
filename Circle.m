classdef Circle < Shape
   properties (SetAccess=private,GetAccess=public)
       radius;
   end
   methods
       function c=Circle(r,n)           
           t =discretize(n);
           title='Circle';
           x =r*cos(t);      y =r*sin(t);  
           x_t=-r*sin(t);    y_t= r*cos(t); 
           x_tt = -r*cos(t); y_tt =-r*sin(t); 
           x_ttt= r*sin(t); y_ttt = r*cos(t);
           c = c@Shape(title,n,t,x,x_t,x_tt,x_ttt,y,y_t,y_tt,y_ttt); % call superclass Shape constructor
           c.radius=r;
       end
       function graph(c)
           graph@Shape(c);
       end 
   end
end