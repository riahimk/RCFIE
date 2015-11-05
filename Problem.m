classdef Problem < handle
    properties
        n;
        type;
        matrix;
        rhs;
        Reg;
        %  preconditionner;
    end
    methods
        function linear_system=Problem(shape,ScParam)
            linear_system.type = ScParam.type;
            linear_system.n    = shape.n;
            %----------------------------------------
            %fprintf('Assembling kernel for SL\n');tic
            SL=Single_Layer(shape,ScParam.k); % instanciate single layer kernels
            %----------------------------------------
            %fprintf('Assembling kernel for DL\n');tic
            DL=Double_Layer(shape,ScParam); % instanciate double layer kernels
            %----------------------------------------
            %fprintf('Assembling kernel for dn_DL\n');tic
            dnDL=DnDouble_Layer(shape,ScParam);
            %----------------------------------------
            %fprintf('Assembling kernel for dn_SL\n');tic
            dnSL=DnSingle_Layer(shape,ScParam);
            %----------------------------------------
            if(ScParam.regu==1)
                SLik              = Single_Layer(shape,ScParam.ik);
                Exp               = exp(-abs(ScParam.ik)*(shape.distance.^4)) ;
                REG               = Exp.*(SL.R .* SLik.K1 + pi/ScParam.n*SLik.K2) ;
                ExpC              = eye(size(Exp)) - Exp;               
                besselH           = pi/ScParam.n * besselh(0, ScParam.ik * ( shape.distance + eye(2*ScParam.n)));
                linear_system.Reg = REG + .25*1i * ExpC .* besselH .* shape.ARCL;
                else
                linear_system.Reg = eye(2*ScParam.n);
            end
            
            %----------------------------------------
            
            %fprintf('Assembling regularization kernels.. \n');
            
            if     (ScParam.type=='Dirichlet ')
                %fprintf('We solve the Scattering problem with::>> Dirichlet boundary condition !\n');
                set_Matrix_Dirichlet(linear_system);
                set_rhs_Dirichlet   (linear_system);
                
            elseif (ScParam.type=='Neumann   ')
                %fprintf('We solve the Scattering problem with::>> Neumann boundary condition !\n');
                set_Matrix_Neumann(linear_system);
                set_rhs_Neumann   (linear_system);
                
            elseif (ScParam.type=='Impedance ')
                %fprintf('We solve the Scattering problem with Impedance boundary condition !\n');
                set_Matrix_Impedance(linear_system);
                set_rhs_Impedance   (linear_system);
                
            elseif (ScParam.type=='GImpedance')
                %fprintf('We solve the Scattering problem with General Impedance boundary condition !\n');
                set_Matrix_GImpedance(linear_system);
                set_rhs_GImpedance   (linear_system);
            end
            
            %--------------------------%
            %  Matrices setting %
            %--------------------------%
            
            function set_Matrix_Dirichlet(linear_system)
                    linear_system.matrix= - SL.R.*( DL.K1*linear_system.Reg + 1i*ScParam.eta*SL.K1 ) ...
                                  - pi/ScParam.n *( DL.K2*linear_system.Reg + 1i*ScParam.eta*SL.K2 ) ...
                                  + linear_system.Reg;
            end
            %-----------------------------------------
            function set_Matrix_Neumann(linear_system)
              Neumann11 = ScParam.k^2 * SL.K1 .* shape.factor - dnDL.K1 ;
              Neumann12 =  - 1i*ScParam.eta * dnSL.K1;
              Neumann21 = ScParam.k^2 * SL.K2 .* shape.factor - dnDL.K2;
              Neumann22 = - 1i*ScParam.eta * dnSL.K2;

             linear_system.matrix =  (dnDL.T + SL.R .* Neumann11 + pi/ScParam.n * Neumann21) ...
                                     *linear_system.Reg ...
                                 + SL.R .* Neumann12 + pi/ScParam.n * Neumann22 ...
                                 + 1i*ScParam.eta * diag(shape.arclength) ;
            end

            %-------------------------------------------
            function set_Matrix_Impedance(linear_system)
              Neumann11 = ScParam.k^2 * SL.K1 .* shape.factor - dnDL.K1 ;
              Neumann12 =  - 1i*ScParam.eta * dnSL.K1;
              Neumann21 = ScParam.k^2 * SL.K2 .* shape.factor - dnDL.K2;
              Neumann22 = - 1i*ScParam.eta * dnSL.K2;
              % - - - - - - - - - - - - - - - - 
              Dirichl11 = - SL.R.*( DL.K1*linear_system.Reg + 1i*ScParam.eta*SL.K1 ) ...
                                  - pi/ScParam.n *( DL.K2*linear_system.Reg + 1i*ScParam.eta*SL.K2 ) ...
                                  + linear_system.Reg;              
             linear_system.matrix =  (dnDL.T + SL.R .* Neumann11 + pi/ScParam.n * Neumann21) ...
                                     *linear_system.Reg ...
                                 + SL.R .* Neumann12 + pi/ScParam.n * Neumann22 ...
                                 + ScParam.zed*Dirichl11 ...                      
                                 + 1i*ScParam.eta * diag(shape.arclength) ;
            end
            %--------------------------%
            %  rigth hand side setting %
            %--------------------------%
            
            function set_rhs_Dirichlet(linear_system)
                %Plane-wave Homogeneous boundary condition u^s + u^i = 0
                %<-> u^s = -u^i
                for i=1:(2*ScParam.n)
                    linear_system.rhs(i,1)= - 2. * exp(1i * ScParam.k * dot([shape.x(i);shape.y(i)], ScParam.d));
                end
            end % end sub-function
            function set_rhs_Neumann(linear_system)
                for i=1:(2*ScParam.n)
                    linear_system.rhs(i,1) = - 2i * ScParam.k * shape.arclength(i) ...
                        * exp(1i * ScParam.k * dot([shape.x(i);shape.y(i)], ScParam.d)) ...
                        * dot(ScParam.d,shape.normal(:,i));
                end
            end
            function set_rhs_Impedance(linear_system)
                for i=1:(2*ScParam.n)
                    linear_system.rhs(i,1) = - 2 *( ...
                        1i* ScParam.k * shape.arclength(i)* dot(ScParam.d,shape.normal(:,i)) ...
                        + ScParam.zed             ) ...
                        * exp(1i * ScParam.k * dot([shape.x(i);shape.y(i)], ScParam.d)) ;
                        
                end
            end
        end% end function constructor
    end
end