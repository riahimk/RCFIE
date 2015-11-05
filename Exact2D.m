classdef Exact2D
    % have a look on http://www.applet-magic.com/helmholtz.htm
    properties
        bv;
    end
    methods
        function ss=Exact2D(Scat_Pb,g,N)
            %
            C=Circle(1,Scat_Pb.n);
            FFT_g = fft(g);
            if(Scat_Pb.type=='Dirichlet ')
            ss.bv  =  FFT_g;
            Coefs  =  - FFT_g./besselh(0,1,Scat_Pb.k);
            
            elseif (Scat_Pb.type=='Neumann   ')
                    %
            else 
                fprintf('\n Boundary conditions not yet supported  !! \n');
            end
        end
    end
end