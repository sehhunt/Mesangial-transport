function [QL, QR, P, Calbumin, CIg, calb_capL, calb_capR,loop_counter, outer_counter] = ...
    cross_capillary(x,y,z, ...
    P_cap_funcL, P_cap_funcR, RL, RR, kL, kR, m, h, Q0R, c0R, Q0L, c0L, ...
    cR_IgA, cL_IgA, direction_left, direction_right, bm_width, ...
    mu, km, kbm, Dalb, DIgA, capzs, capdz, issinglecap, flag)
% Solves the problem including flow into and out of a capillary or adjacent
% capillaries.
% x: not from meshgrid, just the vector of x-values
% y: vector of y values
% z: vector of z values
% P_cap_funcL: a function handle that gives pressure in the left-hand 
% capillary as a function of z
% P_cap_funcR: function handle that gives pressure in right-hand capillary
% as a function of z, in mmHg
% RL is the radius of the lefthand capillary
% RR is the radius of the righthand capilary
% k is the permeability of the capillary wall
% m is the total mesangial width
% Q0R is the volumetric flow rate at entrance of the right capillary
% c0 is the concentration at the entrance of the right capillary
% Q0L is the volumetric flow rate at the left capillary entrance
% c0L is the concentration at the left capillary entrance
% flag is whether the right capillary starts where the left leaves off
% (true if that's the case)
nxnodes = length(x);
nynodes = length(y);
nznodes = length(z);

xbounds = [x(1) x(nxnodes)];
ybounds = [y(1) y(nynodes)];
zbounds = [z(1) z(nznodes)];

dz = z(2)-z(1);
Lz =z(nznodes)-z(1);
%dzstar = dz/Lz;
cap_dz_star = capdz/Lz;
ncap_nodes = length(capzs);

M0L= Q0L*c0L; % The units need to match on these
M0R = Q0R*c0R;

[xmesh,ymesh,zmesh] = meshgrid(x,y,z);

% Fill a vector with the capillary pressure at position i, by
% non-dimensional z
P_capR = [P_cap_funcR(capzs)]; 
P_capL = [P_cap_funcL(capzs)];

%Decide what pressure value to turn to zero
P_maxR = max(P_capR);
P_maxL = max(P_capL);
P_max = max(P_maxR, P_maxL);

% To convert from mmHg to Pascals, and make highest pressure zero
P_capR = (P_capR-P_max)*133.3; % Transformed, and now in Pascals
P_capL = (P_capL-P_max)*133.3;
P_urinary = 1.3-P_max; % mmHg

% First assume that mesangium is a wall - no flow into or out of capillary
% Calculate the function that is equal to dc/dz in this situation

% We need this for both the albumin concentration and the IgA concentration
dcdzL = no_mes_concentration(RL,kL,Q0L,c0L,m, Lz,P_capL, P_urinary, ...
    P_max, direction_left);
dcIgAdzL = no_mes_concentration(RL,kL,Q0L,cL_IgA,m,Lz,P_capL,P_urinary,...
    P_max, direction_left, true);

c_func_zL = calculate_cap_concentration_no_mes(cap_dz_star,dcdzL,ncap_nodes,direction_left); %dimless, albumin
c_IgA_func_zL = calculate_cap_concentration_no_mes(cap_dz_star, dcIgAdzL, ncap_nodes,direction_left,c_func_zL*c0L); % dimless, IgA
%*********************************************************
% Returns and requires DIMLESS c (as input and output)
%*********************************************************

if flag
    dQdz_L = calculate_capillary_Q_no_mes(kL,RL,Lz,Q0L,m,direction_left,c0L,...
        P_capL,P_urinary,P_max);
    Q_L = calculate_capillary_flow(dQdz_L,c_func_zL,cap_dz_star,ncap_nodes,direction_left); %dimless
    Q0R = Q_L(ncap_nodes)*Q0L;
    c0R = c_func_zL(ncap_nodes)*c0L; 
end

if flag
    cR_IgA = c_IgA_func_zL(ncap_nodes)*cL_IgA;
end

if ~issinglecap

    % Calculate concentration function in right capillary for albumin and
    % IgA
    dcdzR = no_mes_concentration(RR,kR,Q0R,c0R,m,Lz,P_capR, P_urinary, ...
        P_max, direction_right);
    dcIgAdzR = no_mes_concentration(RR,kR,Q0R,cR_IgA,m,Lz,P_capR, ...
        P_urinary,P_max, direction_right,true);
end


% Now, use dc/dz to calculate c(z) down the length of the capillary
if ~issinglecap
    c_func_zR = calculate_cap_concentration_no_mes(cap_dz_star,dcdzR,ncap_nodes,direction_right);%dimless, albumin
    c_IgA_func_zR = calculate_cap_concentration_no_mes(cap_dz_star,dcIgAdzR,ncap_nodes,direction_right,c_func_zR*c0R); % dimless, IgA
end



% These are vectors that give concentration in capillary at z position by
% nodes along the z-direction

% We need to fill in these concentrations and pressures to be passed to
% convection_diffusion in a vector filled by global node number

mes_conc_L = zeros(nynodes*nznodes,1); % because lefthand boundary is nodes 1 - ny*nz
IgA_mes_conc_L = zeros(nynodes*nznodes,1);
P_mes_L = zeros(nynodes*nznodes,1);

row_nums = nynodes*nznodes*(nxnodes-1)+1:1:nynodes*nznodes*nxnodes;
col_nums = ones(size(row_nums));
conc_filler = zeros(size(row_nums));
IgA_filler_R = zeros(size(row_nums));

PmesR_filler = zeros(size(row_nums));
multiplier = floor(dz/capdz);
for i=1:nznodes
    mes_conc_L((i-1)*nynodes+1:nynodes*i) = c0L*c_func_zL(1+(i-1)*multiplier);
    IgA_mes_conc_L((i-1)*nynodes+1:nynodes*i) = ...
        cL_IgA*c_IgA_func_zL(1+(i-1)*multiplier);
    if ~issinglecap
        conc_filler(nynodes*(i-1)+1:...
            nynodes*i) = c0R*c_func_zR(1+(i-1)*multiplier);
        IgA_filler_R(nynodes*(i-1)+1:nynodes*i) = ...
            cR_IgA*c_IgA_func_zR(1+(i-1)*multiplier);
        PmesR_filler(nynodes*(i-1)+1:nynodes*i)=P_capR(1+(i-1)*multiplier);
    end
    
    P_mes_L((i-1)*nynodes+1:nynodes*i) = P_capL(1+(i-1)*multiplier);
end

if ~issinglecap
    mes_conc_R = sparse(row_nums,col_nums,conc_filler);
    P_mes_R = sparse(row_nums,col_nums,PmesR_filler);
    IgA_conc_R = sparse(row_nums,col_nums,IgA_filler_R);
    IgA_conc_R = full(IgA_conc_R);
    mes_conc_R = full(mes_conc_R);
    P_mes_R = full(P_mes_R);
end

% Initialize the boundary conditions to pass to convection_diffusion
% What's the boundary condition for pressure and concentration on boundary 1?
bc_p1 = {'N' 0};
bc_c1 = {'N' 0};
bc_Ig1 = {'N' 0};
% Same for boundaries 2 - 6

% Create the boundary condition function for total pressure that depends on
% osmotic pressure

% Set constants for osmotic pressure function
a1 = 1.63; % mmHg / (g/dL)
a2 = 0.294;

% This creates an anonymous function, total_pressure, that takes one
% argument (conc). Concentration should be a vector of concentration
% values. total_pressure is a function handle

total_pressure = @(conc) (P_urinary+...
    ((a1.*conc + a2.*conc.^2)))*133.3; % Pa
% 1.333*10^-4 is a conversion factor from mmHg to kg/(um *s^2)

% Now use this function in the boundary condition
%bc_p2 = {'N' 0};
bc_p2 = {'N' total_pressure};
%bc_p2 = {'D' 1.3*1.333*10^(-4)}; % kg/(um s^2)
bc_c2 = {'N' 0};
bc_Ig2 = {'N' 0};

%bc_p3 = {'D' 35*133.3}; % Pa 
bc_p3 = {'D' P_mes_L};
bc_c3 = {'D' mes_conc_L}; % g/100 mL (g/dL)
bc_Ig3 = {'D' IgA_mes_conc_L}; % g/100 mL

if ~issinglecap
    bc_p4 = {'D' P_mes_R};
    bc_c4 = {'D' mes_conc_R};
    bc_Ig4 = {'D' IgA_conc_R};
else
    bc_p4 = {'N' 0};
    bc_c4 = {'N' 0};
    bc_Ig4 = {'N' 0};
end


bc_p5 = {'N' 0};
bc_c5 = {'N' 0};
bc_Ig5 = {'N' 0};

bc_p6 = {'N' 0};
bc_c6 = {'N' 0};
bc_Ig6 = {'N' 0};

pressure_bcs = {bc_p1; bc_p2; bc_p3; bc_p4; bc_p5; bc_p6};
albumin_bcs = {bc_c1; bc_c2; bc_c3; bc_c4; bc_c5; bc_c6};
Ig_bcs = {bc_Ig1; bc_Ig2; bc_Ig3; bc_Ig4; bc_Ig5; bc_Ig6};

NNODES = nynodes*nxnodes*nznodes;

error = 1;
limit = 0.001;

% Now run the old problem to get mesangial pressure and concentration
% distribution

[P,Calbumin,CIg,loop_counter] = convection_diffusion(xmesh,ymesh,zmesh,...
    xbounds,ybounds,zbounds, pressure_bcs, albumin_bcs, Ig_bcs, bm_width,...
    mu,kbm,km,Dalb,DIgA,issinglecap);

outer_counter = 1;

while error > limit
    
    % Make vectors from the matrices
    Pvector = make_Vector(P);
    Calb_vector = make_Vector(Calbumin);
    CIg_vec = make_Vector(CIg);
    
    %****************************************************************
    % The above vectors are in terms of the mesh global node number
    % Now we need to interpolate these values onto the capillary zmesh
    %******************************************************************
    
    % Make vectors a function of z-position, not global node number
    Pcap_zfuncL = fill_boundary_function(ncap_nodes,...
        Pvector(1:nynodes:nynodes*nznodes),capdz,dz,z,capzs,false);
%     Calb_zfuncL = Calb_vector(1:nynodes:nynodes*nznodes);
%     CIg_zfuncL = CIg_vec(1:nynodes:nynodes*nznodes);
    
    Pcap_zfuncR = fill_boundary_function(ncap_nodes,...
        Pvector(nynodes*nznodes*(nxnodes-1)+1:...
        nynodes:nynodes*nxnodes*nznodes), capdz, dz,z,capzs,false);
%     Calb_zfuncR = Calb_vector(nynodes*nznodes*(nxnodes-1)+1:...
%         nynodes:nynodes*nxnodes*nznodes);
%     CIg_zfuncR = CIg_vec(nynodes*nznodes*(nxnodes-1)+1:...
%         nynodes:nynodes*nxnodes*nznodes);
    plasma_fluxL = calculate_plasma_flux(Pvector,km,mu,xbounds(1),true,...
        xmesh,ymesh,zmesh);
    plasma_fluxR = calculate_plasma_flux(Pvector,km,mu,xbounds(2),false,...
        xmesh,ymesh,zmesh);
    % These are two vectors that have km/mu*dp/dx at the cell centers
    % moving down the the z-coordinate (or -km/mudpdx)
    
    % This fills in vectors to convert the plasma_fluxes to functions of
    % z-position
    plasma_fluxL_z = fill_boundary_function(ncap_nodes,plasma_fluxL,...
        capdz,dz,z,capzs,true);
    plasma_fluxR_z = fill_boundary_function(ncap_nodes,plasma_fluxR,...
        capdz,dz,z,capzs,true);
    
    dQdzL = with_mes_concentration(RL,kL,h,Lz,Q0L,c0L,m,Pcap_zfuncL,...
        plasma_fluxL_z,P_urinary,P_max,direction_left);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test alternative organization of concentration flux calculation     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dQcdzL = calculate_concentration_flux(Pvector,Calb_vector,h,Lz,Q0L,c0L,km,mu,Dalb,...
    %    xbounds(1),direction_left, true,xmesh,ymesh,zmesh);
    %dQIgdzL = calculate_concentration_flux(Pvector,CIg_vec, h,Lz,Q0L,...
    %    cL_IgA,km,mu,DIgA,xbounds(1),direction_left, true,xmesh,ymesh,zmesh);
    
%     dQcdzL = fill_boundary_function(ncap_nodes,dQcdzL,capdz,dz,z,capzs,...
%         true);
%     dQIgdzL = fill_boundary_function(ncap_nodes,dQIgdzL,capdz,dz,z,capzs,...
%         true);
    
    concentration_flux_alb_L = calculate_concentration_flux2(Pvector,Calb_vector,...
        km,mu,Dalb,xbounds(1),true, xmesh,ymesh,zmesh);
    concentration_flux_IgA_L = calculate_concentration_flux2(Pvector, CIg_vec,...
        km,mu,DIgA,xbounds(1),true,xmesh,ymesh,zmesh);
    
    alb_flux_L = fill_boundary_function(ncap_nodes,concentration_flux_alb_L,...
        capdz,dz,z,capzs,true);
    IgA_flux_L = fill_boundary_function(ncap_nodes,concentration_flux_IgA_L,...
        capdz, dz,z,capzs,true);
    
    dQcdzL = calculate_capillary_flux(Q0L,c0L,Lz,h,alb_flux_L);
    dQIgdzL = calculate_capillary_flux(Q0L,cL_IgA,Lz,h,IgA_flux_L);
    
    [QL, cap_QcL, calb_capL ] = calculate_capillary_values(c0L,cap_dz_star,ncap_nodes,dQdzL,...
        dQcdzL, direction_left);
    [~,~,cIgA_capL ] = calculate_capillary_values(cL_IgA,cap_dz_star,...
        ncap_nodes,dQdzL,dQIgdzL, direction_left, calb_capL);
    
    if flag
        c0R = calb_capL(ncap_nodes)*c0L;
        cR_IgA = cIgA_capL(ncap_nodes)*cL_IgA;
        Q0R = QL(ncap_nodes)*Q0L;
    end
    
    
    dQdzR = with_mes_concentration(RR,kR,h,Lz,Q0R,c0R,m,Pcap_zfuncR,...
        plasma_fluxR_z,P_urinary,P_max,direction_right);
    %*****************************************************
    % These functions of DIMLESS c
    % 
    %******************************************************
    
    % Use alternate concentration flux calculations
%     dQcdzR = calculate_concentration_flux(Pvector,Calb_vector,h,Lz, Q0R,c0R,km,mu,Dalb,...
%         xbounds(2),direction_right, false,xmesh,ymesh,zmesh);
%     
%     
%     dQIgdzR = calculate_concentration_flux(Pvector,CIg_vec,h,Lz,Q0R,...
%         cR_IgA,km,mu,DIgA,xbounds(2),direction_right, false,xmesh,ymesh,zmesh);

conc_flux_alb_R = calculate_concentration_flux2(Pvector,Calb_vector,km,mu,Dalb,...
    xbounds(2),false,xmesh,ymesh,zmesh);


conc_flux_IgA_R = calculate_concentration_flux2(Pvector,CIg_vec,km,mu,DIgA,...
    xbounds(2),false,xmesh,ymesh,zmesh);

alb_flux_R = fill_boundary_function(ncap_nodes,conc_flux_alb_R,...
    capdz,dz,z,capzs,true);
IgA_flux_R = fill_boundary_function(ncap_nodes,conc_flux_IgA_R,...
    capdz,dz,z,capzs,true);
    %**************************************************************
    % Qc function above takes is only a function of z
    % These functions return fluxes at center of elements in z direction, so
    % we have to convert them to fill the capillary z-nodes
    %**************************************************************
    
    % The above are functions of zelement, not z-postion
    % Convert to functions of z-position
   
    dQcdzR = calculate_capillary_flux(Q0R,c0R,Lz,h,alb_flux_R);
    
    
%     dQcIgdzL = calculate_concentration_flux(Pvector,CIg_vec,h,km,mu,DIgA,...
%         xbounds(1),true,xmesh,ymesh,zmesh);
%     dQcIgdzR = with_mes_concentration(Pvector,CIg_vec,h,km,mu,DIgA,xbounds(2),...
%         false,xmesh,ymesh,zmesh);
%     
    
    dQIgdzR = calculate_capillary_flux(Q0R,cR_IgA,Lz,h,IgA_flux_R);
    
    
    
    [QR, cap_QcR, calb_capR ] = calculate_capillary_values(c0R,cap_dz_star,ncap_nodes,dQdzR,...
        dQcdzR, direction_right);
    
    [~,~,cIgA_capR] = calculate_capillary_values(cR_IgA,cap_dz_star,...
        ncap_nodes,dQdzR,dQIgdzR, direction_right, calb_capR);
    %**********************************************************
    % The above lines return DIMENSIONLESS cs
    % These are returned at all the CAPILLARY nodes, not the mesangial
    % nodes
    %**************************************************************
    
    % Fill in updated concentrations from the capillary solve above
    for i=1:nznodes
        mes_conc_L2((i-1)*nynodes+1:i*nynodes)=c0L*calb_capL(1+(i-1)*multiplier);
        mes_conc_R2(nynodes*nznodes*(nxnodes-1)+1+(i-1)*nynodes:...
            nynodes*nznodes*(nxnodes-1)+i*nynodes) = c0R*calb_capR(1+(i-1)*multiplier);
        
        mes_conc_IgAL((i-1)*nynodes+1:i*nynodes) = cL_IgA*cIgA_capL(1+(i-1)*multiplier);
        mes_conc_IgAR(nynodes*nznodes*(nxnodes-1)+1+(i-1)*nynodes:...
            nynodes*nznodes*(nxnodes-1)+i*nynodes) = cR_IgA*cIgA_capR(1+(i-1)*multiplier);
    end
    
    % Use new concentrations to update boundary conditions
    bc_c3 = {'D' mes_conc_L2}; % g/100 mL (g/dL)
    bc_Ig3 = {'D' mes_conc_IgAL};
    if ~issinglecap
        bc_c4 = {'D' mes_conc_R2};
        bc_Ig4 = {'D' mes_conc_IgAR};
    end
    
    albumin_bcs = {bc_c1; bc_c2; bc_c3; bc_c4; bc_c5; bc_c6};
    Ig_bcs = {bc_Ig1; bc_Ig2; bc_Ig3; bc_Ig4; bc_Ig5; bc_Ig6};
    
    [P2, Calb2, CIg2, loop_counter] = convection_diffusion(xmesh,ymesh,...
        zmesh,xbounds,ybounds,zbounds,pressure_bcs,albumin_bcs,Ig_bcs,...
        bm_width,mu,kbm,km,Dalb,DIgA,issinglecap);
    
    Cvec2 = make_Vector(Calb2);
    error = (Cvec2 - Calb_vector).^2;
    error = max(error);
    Calbumin = Calb2;
    P = P2;
    CIg = CIg2;
    
    outer_counter= outer_counter+1;
end

QL = QL.*10^(-6); % converts to nL/s
QR = QR.*10^(-6);
end
