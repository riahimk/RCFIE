function [w] = Plane_wave(C,Scat_pb)
    d_dot_xy = [C.x;C.y]'*Scat_pb.d ;
 w = exp(1i*2*pi*Scat_pb.k * d_dot_xy);
end