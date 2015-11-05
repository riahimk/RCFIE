classdef DnSingle_Layer < Integral_Operator
    %properties (SetAccess=private,GetAccess=public)
    %    K1;K2;
    %end
    methods
        function dnSL=DnSingle_Layer(shape,pb)
            dnSL=dnSL@Integral_Operator(shape);
            set_Kernel(dnSL,shape,pb);
        end
        function set_Kernel(dnSL,shape,ScParam)
            n = dnSL.n;
            t = discretize(n);
            %bess_h_0 = besselh(0, ScParam.k * (shape.distance+eye(2*n)));
            bess_h_1 = besselh(1, ScParam.k * (shape.distance+eye(2*n)));
            
            for i =1:(2*n)
                for j =1:(2*n)
                    if (shape.distance(i,j)) < 1.e-09
                        dnSL.K1(i,j) = 0.;
                        dnSL.K2(i,j) = dot(shape.normal(:,i),[shape.x_tt(i);shape.y_tt(i)])/(2 * pi * shape.arclength(i));
                    else
                        dnSL.K1(i,j) = ScParam.k/2/pi*dot(shape.normal(:,i),[shape.x(i)-shape.x(j);shape.y(i)-shape.y(j)]) * ...
                            real(bess_h_1(i,j))/shape.distance(i,j) * shape.arclength(j);
                        dnSL.K2(i,j) = -1i * ScParam.k/2 * dot(shape.normal(:,i),[shape.x(i)-shape.x(j);shape.y(i)-shape.y(j)]) ...
                            *bess_h_1(i,j)/shape.distance(i,j) * shape.arclength(j) ...
                            - dnSL.K1(i,j) * log(4*sin(t(i)/2 - t(j)/2 + (i==j) )^2  );
                    end;
                end % for j
            end % for i
            %H1 = dnSL.K1
            %H2 = dnSL.K2
        end % constructor 
    end % methods
end