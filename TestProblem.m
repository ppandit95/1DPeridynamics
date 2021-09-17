%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Test Program on Longitudinal Vibration of a bar            %%
%%%                                       --Pushkar Anirudha Pandit      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                                       %Starts counting the process time
clc                                       %Clears the Screen
clear                                     %Clears the memory
%Simulation Input Parameters
N = 1000;                                %Number of Material Points
NT = 26000;                               %Number of time steps
e_moduli = 200e9;                         %Elastic Moduli of the Bar
length = 1.0;                             %Length of the Bar
density = 7850;                           %Density of the Bar
dx = length/N;                            %Spacial Distance between the Material Points
pos = dx/2 : dx : length-dx/2;            %Setting up position vector which contains coordinate of each material point in the solid body
bndpos = [-dx/2,-3*dx/2,-5*dx/2];         %Positions of fictitious material points which enforce dirichlet BC
pos = cat(2,pos,bndpos);                  %Concatenating both arrays so as to have values of coordinates for all the material points whether real or ficticious
disp = 0.001*pos;                         %Setting the displacement of each material point for a uniformly stretched bar(BC has not been enforced)
delta = 3.015*dx;                         %Setting up the Horizon Size in terms of spacing
bond_const = 2*e_moduli/(dx^2 * delta^2); %Calculating bond constant for isotropic homogeneous 1D bar
dt = 0.8*sqrt(2*density*dx/(2*delta*(dx^2)*bond_const)); %Calculating the time intervel for stable time marching as per the Von Neumann Criteria
%Matrix Initialization
vel = zeros(1,(size(pos,2)-3));           %Setting up a vector to store velocity at each material point
HIndex = zeros(7,size(pos,2));            %Stores the index of material points in pos vector that are present in Horizon of each material point corresponding to each column
Xstate = zeros(7,size(pos,2));            %Stores the relative poistions of points in the horizon of each material point corresponding to each column
Ustate = zeros(7,size(pos,2));            %Stores the relative displacements of points in the horizon of each material point corresponding to each column
Ystate = zeros(7,size(pos,2));            %Stores the relative deformed positions of points in the horizon of each material points corresponding to each column
Tstate = zeros(7,size(pos,2));            %Stores the force contribution for each bond in the horizon of each material point corresponding to each column
forcevec = zeros(1,size(pos,2));          %Stores the net force contribution for each material point in a vector
strain_energy_CCM = 0.5*e_moduli*1e-6*ones(1,size(pos,2));%Stores the strain energy density corresponding to each material point calculated using classical continuum mechanics
strain_energy_PD = zeros(1,size(pos,2));  %Stores the strain energy density corresponding to each material point calculated using peridynamics theory
surface_correction = ones(1,size(pos,2)); %Stores the surface correction factor corresponding to each material point 
surface_correction_pairwise = zeros(7,size(pos,2));%Stores the suface correction factor corresponding to each bond force in the horizon of each material point corresponding to each column 
volume_correction = zeros(7,size(pos,2)); %Stores the volume correction factor corresponding to each neighbor in the horizon of each material point corresponding to each column
PDDisp = zeros(2,NT);                     %A post processing variable which stores the displacement variation and corresponding time for a specific point in a 1D bar

%Determination of Horizon for each Material Point

for i = 1:size(pos,2)
    horizon = find(pos <= (pos(1,i) + delta) & pos >= (pos(1,i) - delta));
    horizon = horizon(horizon ~= i);
    HIndex(1:size(horizon',1),i) = horizon';
end

%Setting up X and U states 
for j=1:size(pos,2)
    for i = 1:6
        if HIndex(i,j) == 0
            break;
        else
            Xstate(i,j) = pos(1,HIndex(i,j)) - pos(1,j);
            Ustate(i,j) = disp(1,HIndex(i,j)) - disp(1,j);
        end
    end 
end

%Calculating Y state and stretch
Ystate = Xstate + Ustate;
stretch = abs(Ystate) - abs(Xstate);
stretch = stretch ./ abs(Xstate);



%Calculation of Volume Correction Factor
volume_correction = ((abs(Xstate) <= (delta - dx/2)) & (abs(Xstate) ~= 0));
dummy =  ((delta+dx/2) - abs(Xstate).*((abs(Xstate) <= (delta+dx/2))&(abs(Xstate) > (delta - dx/2)))) ./ dx;
dummy(dummy ==(delta+dx/2)/dx) = 0;
volume_correction = volume_correction + dummy;

%Calculation of Surface Correction Factor
dummy1 = (0.25 * bond_const * dx^3) .*(stretch .^ 2).* abs(Xstate) .* volume_correction;
dummy1(isnan(dummy1)) = 0;
strain_energy_PD = sum(dummy1,1);
surface_correction = (strain_energy_CCM ./ strain_energy_PD)  ;

for j=1:size(pos,2)
    for i = 1:6
        if HIndex(i,j) == 0 
            break;
        else
            surface_correction_pairwise(i,j) = (surface_correction(1,HIndex(i,j)) + surface_correction(1,j))/2;
        end
    end 
end

%Setting Up Boundary Condition
disp((size(pos,2)-2):size(pos,2)) = 0;
%Marching in Time
for tt=1:NT
    %Setting up X and U states 
    for j=1:size(pos,2)
        for i = 1:6
            if HIndex(i,j) == 0
                break;
            else
                Ustate(i,j) = disp(1,HIndex(i,j)) - disp(1,j);
            end
        end 
    end

    %Calculating Y state and stretch
    Ystate = Xstate + Ustate;
    stretch = abs(Ystate) - abs(Xstate);
    stretch = stretch ./ abs(Xstate);

    %Calculating force contribution for each bond pair in each horizon
    p = Ystate ./ abs(Ystate);
    Tstate =  bond_const .* stretch .* p .* volume_correction .* surface_correction_pairwise;
    Tstate(isnan(Tstate)) = 0;
    forcevec = sum(Tstate,1).*dx^3 ;
    forcevec((size(pos,2)-2):size(pos,2)) = []; % Trimming the last 3 elements since these are contributions for ficticious material points


    acc = forcevec ./ density;

    %Update displacement and velocity
    vel = vel + acc .* dt;
    disp(1,1:(size(pos,2)-3)) = disp(1,1:(size(pos,2)-3)) + vel .* dt;
    PDDisp(1,tt) = disp(1,floor(N/2));
    PDDisp(2,tt) = tt*dt;

end
plot(PDDisp(2,:),PDDisp(1,:))
title('Displacement Variation with Time at x = 0.5 m')
xlabel('Time')
ylabel('Displacement')
toc
