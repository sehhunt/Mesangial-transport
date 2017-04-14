function [ cap_Q, cap_Qc, cap_c ] = calculate_capillary_values(c0, dz, ...
    nnodes, dQ_func, dQc_func, direction, varargin)
%calculate_capillary_values Returns two matrices (cap_Q, cap_c) containing
% the volumetric flow rate as a function of z-position in the capillary and
% the concentration of albumin as a function of position in the capillary
%   Q0 is the inlet volumetric flow
%   c0 is the inlet albumin concentration
%   dz is the z-node spacing
%   nnodes is the number of nodes in the z-direction
%   dQ_function is a function of capillary concentration and z position
%   that gives dQ/dz - this is dimensionless
%   Cvector gives concentration as a function of zposition along the
%   capillary - no dimesionless
%
%   dQc_func is a function of z position that
%   gives d(Qc)/dz - dimensionless

cap_c = zeros(nnodes,1);
cap_Qc = zeros(nnodes,1);
cap_Q = zeros(nnodes,1);

if length(varargin)==1
    isIgA = true;
    calb_func = varargin{1};
else
    isIgA = false;
end
% Need to know where the inlet is

if direction > 0
    cap_c(1) = 1;  % Initialize inlet values
    cap_Qc(1) = 1; % Equal to Q/Q0*c/c0
    cap_Q(1) = 1;
 
  for i=2:nnodes
      if ~isIgA
          cap_Q(i) = dQ_func(cap_c(i-1),i-1)*dz+cap_Q(i-1);
          cap_Qc(i) = dz*dQc_func(i-1)+cap_Qc(i-1);
          cap_c(i) = (cap_Qc(i)/cap_Q(i));
          
      else
          cap_Q(i) = dQ_func(calb_func(i-1),i-1)*dz+cap_Q(i-1);
          cap_Qc(i) = dz*dQc_func(i-1)+cap_Qc(i-1);
          cap_c(i) = (cap_Qc(i)/cap_Q(i));
      end
   end  
else
    cap_c(nnodes) = 1;
    cap_Qc(nnodes) = 1;
    cap_Q(nnodes) = 1;
    
    for i=nnodes-1:-1:1
        if ~isIgA
            cap_Q(i) = -dQ_func(cap_c(i+1),i+1)*dz+cap_Q(i+1);
            cap_Qc(i) = -dz*dQc_func(i+1)+cap_Qc(i+1);
            cap_c(i) = (cap_Qc(i)/cap_Q(i));
        else
            cap_Q(i) = -dQ_func(calb_func(i+1),i+1)*dz+cap_Q(i+1);
            cap_Qc(i) = -dz*dQc_func(i+1)+cap_Qc(i+1);
            cap_c(i) = (cap_Qc(i)/cap_Q(i));
        end
    end
end



end

