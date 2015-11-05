classdef Single_Layer < Integral_Operator & handle
    %properties (SetAccess=private,GetAccess=public)
    %    K1;K2
    %end
    methods
        function SL=Single_Layer(shape,wavenumber)  %------ Class constructor call base class="superclass"
            SL=SL@Integral_Operator(shape);
            set_R(SL);
            set_Kernel(SL,shape,wavenumber); % setting the properties by the following handled-setter method
        end
         function set_Kernel(SL,shape,wavenumber)
             n=SL.n;
             t=discretize(n);
             bess_h_0=besselh(0, wavenumber * (shape.distance+eye(2*n)));
             for i =1:(2*n)
                 for j =1:(2*n)
                     if  shape.distance(i,j) < 1.e-09  % diagonal terms
                         SL.K1(i,j) = -1/2/pi* shape.arclength(j);
                         SL.K2(i,j) = (1i/2 - SL.gamma/pi -1/pi * log(wavenumber/2 * shape.arclength(j)))* shape.arclength(j);
                     else
                         SL.K1(i,j) = -1/2/pi * real(bess_h_0(i,j))* shape.arclength(j);
                         SL.K2(i,j) = 1i/2 * bess_h_0(i,j)* shape.arclength(j) ...
                                        - SL.K1(i,j)*log( 4*sin( t(i)/2-t(j)/2 )^2 + (i==j) );            
                     end
                 end
             end 
             
         end
    end
end
