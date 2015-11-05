classdef Double_Layer < Integral_Operator
    %properties (SetAccess=private,GetAccess=public)
    %    K1;K2;
    %end
    methods 
        function DL=Double_Layer(shape,pb)
            DL=DL@Integral_Operator(shape);
            set_Kernel(DL,shape,pb);
        end
        function set_Kernel(DL,shape,ScParam)
            n = DL.n;
            t = discretize(n);
            bess_h_1 = besselh(1, ScParam.k * (shape.distance+eye(2*n)));
            for i =1:(2*n)
                for j =1:(2*n)
                    
                    if (shape.distance(i,j)) < 1.e-09  % diagonal terms 
                        DL.K1(i,j) = 0.;
                        DL.K2(i,j) = - ( shape.normal(:,i)'* [ shape.x_tt(i);shape.y_tt(i) ] ) ...
                                     /2/pi/shape.arclength(i)^2 ;
                    else
                        DL.K1(i,j) = ScParam.k/2/pi * ...
                            (shape.normal(:,j)'*[shape.x(i)-shape.x(j);shape.y(i)-shape.y(j)]) * ...
                            real(bess_h_1(i,j))/shape.distance(i,j);
                        DL.K2(i,j) = -1i*ScParam.k/2 * ...
                                    ( shape.normal(:,j)'*[shape.x(i)-shape.x(j);shape.y(i)-shape.y(j)] ) * ...
                                    bess_h_1(i,j)/shape.distance(i,j) ...
                                    - DL.K1(i,j) * log(4*sin(t(i)/2-t(j)/2 + (i==j) )^2  );
                    end
                end
            end
        end
    end
end