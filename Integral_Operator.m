classdef Integral_Operator < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected, GetAccess=public)
        n; R; T;
        K1;K2;
    end
    properties (Constant)
        gamma = 0.5772156649015328606;
    end
    methods
        function I=Integral_Operator(shape)%---constructor
            I.n=shape.n;
        end
        
        function set_R(I)
            n=I.n;
            m=1:1:n-1;
            for i =1:(2*n)
                for j =1:(2*n)
                    R(i,j) = - 2 * pi/n * sum((cos(m *abs(i-j) *pi/n))./m) - (-1)^abs(i-j) * pi/n/n;
                    %R(i,j) = - 2 * pi/n * sum((cos(m *j *pi/n))./m) - (-1)^j * pi/n/n;
                end
            end
            I.R=R;
        end
        function set_T(I)
            n=I.n;
            for i =1:(2*n)
                for j =1:(2*n)
                    if abs(i-j)==0
                        T(i,j) =-n/2;
                    elseif rem(abs(i-j),2) == 0
                        T(i,j) = 0;
                    else
                        T(i,j) = 1/(2 * n * (sin(abs(i-j) * pi/2/n))^2);
                    end
                end
            end
            I.T=T;
        end
    end %------------- end methods
end