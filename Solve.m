classdef Solve < handle
    properties
        density;
        Regdensity;
        %near_field;
        far_field;
    end
    methods
        function O=Solve(shape,lin_sys,ScParam)

            O.density    = lin_sys.matrix\(lin_sys.rhs);
            O.Regdensity = lin_sys.Reg   *  O.density;

            set_far_field(shape,ScParam,O);
            function  set_far_field(shape,ScParam,O)
                theta = 0 : pi/ScParam.n: 2*pi-pi/ScParam.n;
                ffld = .0;
                for i =1:(2*ScParam.n)
                    drc = [cos(theta(i)); sin(theta(i))];
                    
                        for l=1:(2*ScParam.n)
                            term = ( ...
                                     ScParam.k   * dot(drc,shape.normal(:,l)) * O.Regdensity(l) ...
                                   + ScParam.eta * shape.arclength(l)         * O.density(l)    ...
                                   ) * exp( -1i * ScParam.k * dot(drc,[shape.x(l);shape.y(l)] ) );
                            ffld = ffld + term;
                        end

                    O.far_field(i) = ffld*pi/ScParam.n * exp(-1i * pi/4)/sqrt(8 * pi * ScParam.k);
                    ffld=.0;
                end
            end
        end
    end % end methods
end % end class
