%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          PD Program to Simulate Vibrations in a 1D Stepped Bar        %%
%%           Steel+Brass                      --Pushkar Anirudha Pandit  %%
%%  E_brass/E_steel = 0.485,rho_brass/rho_steel = 1.083                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc 
clear

tic 
%%Program Input%%
ndivx = 1000;
nbnd = 3;
totnode = ndivx + nbnd;
nt = 26000;
maxfam = 100;

%%Matrix Initialisation%%
coord = zeros(totnode,1);
numfam = zeros(totnode,1);
pointfam = zeros(totnode,1);
pforce = zeros(totnode,1);
bforce = zeros(totnode,1);
stendens = zeros(totnode,1);
fncst = ones(totnode,1);
fncstold = ones(totnode,1);
disp = zeros(totnode,1);
vel = zeros(totnode,1);
acc = zeros(totnode,1);

andisp = zeros(nt,1);
pddisp = zeros(nt,1);
pddisp2 = zeros(nt,1);
pdtime = zeros(nt,1);

nodefam = zeros(totnode*totnode,1);

%%Initializing Parameters%%
length = 1.0;
dx = length / ndivx;
delta = 3.015 * dx;
dens = 7850;
emod = 200e9;
area = dx * dx;
vol = area * dx;
bc = 2*emod / (area*delta*delta);
sedload1 = 0.5*emod*1.0e-6;
dt = 0.8*sqrt((2*dens*dx)/(2*delta*area*bc));
totime = nt*dt;
ctime = 0.0;
idist = 0.0;
fac = 0.0;
radij = dx / 2.0;
cnode = 0;
nlength = 0;
dforce1 = 0.0;
ntotrao = 20;
cwave = sqrt(emod/dens);

%Specifying location of Material Points
coord(1:ndivx,1) = dx/2 + ((0:ndivx - 1).*dx)';


%Material points of the constrained region
coord(ndivx+1:end,1) = (-0.5*dx) - ((0:nbnd-1).*dx)';

%Determination of material points inside horizon of each material point
for i=1:totnode
    if i==1
        pointfam(i,1) = 1;
    else
        pointfam(i,1) = pointfam(i-1,1)+numfam(i-1,1);
    end
    for j=1:totnode
        idist = sqrt((coord(j,1)-coord(i,1))^2);
        if i~=j
            if idist <= delta
                numfam(i,1) = numfam(i,1)+1;
                nodefam(pointfam(i,1)+numfam(i,1)-1,1) = j;
            end
        end
    end
end

%Determination of Surface Correction Factors
for i=1:totnode
    disp(i,1) = 0.001*coord(i,1);
end 

for i=1:totnode
    stendens(i,1) = 0.0;
    for j=1:numfam(i,1)
        cnode = nodefam(pointfam(i,1)+j-1,1);
        idist = sqrt((coord(cnode,1)-coord(i,1))^2);
        nlength = sqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))^2);
        if idist <= (delta-radij)
            fac = 1;
        elseif idist <= delta+radij
            fac = (delta+radij-idist)/(2*radij);
        else
            fac = 0;
        end
%         if i > 500
%             stendens(i,1) = stendens(i,1) + 0.25*0.485*bc*((nlength-idist)/idist)^2*idist*vol*fac;
%         else
            stendens(i,1) = stendens(i,1) + 0.25*bc*((nlength-idist)/idist)^2*idist*vol*fac;
%         end
    end
    %Calculation of surface correction factor in x direction
    %by finding the ratio of the analytical strain energy density value
    %to the strain energy density value obtained from PD theory
%     if i <= 500
        fncst(i,1) = sedload1/stendens(i,1);
        
%     else
%          fncst(i,1) = 0.485*sedload1/stendens(i,1);
%     end
end


%Initialization of displacements and velocities
vel = zeros(totnode,1);
disp = zeros(totnode,1);

%Initial Condition
for i=1:ndivx
    vel(i,1) = 0.0;
    disp(i,1) = 0.001*coord(i,1);
end

%Boundary Condition - Zero displacement at x=0
for i=(ndivx+1):totnode
    vel(i,1) = 0;
    disp(i,1) = 0;
end


%Time Integration
for tt=1:nt
    fprintf('tt = %d\n',tt);
    ctime = tt*dt;
    for i=1:ndivx
        pforce(i,1) = 0;
        for j=1:numfam(i,1)
            cnode = nodefam(pointfam(i,1)+j-1,1);
            idist = sqrt((coord(cnode,1)-coord(i,1))^2);
            nlength = sqrt((coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))^2);
            
            %Volume Correction
            if idist <= (delta-radij)
                fac = 1;
            elseif idist <= (delta+radij)
                fac = (delta+radij-idist)/(2*radij);
            else
                fac = 0;
            end
            
            %Determination of surface correction between two material
            %points
            scr = (fncst(i,1)+fncst(cnode,1))/2;
            
            %Calculation of the peridynamic Force in the x direction
            %acting on a material point i due to the material point j
%             if i > 500
%                 dforce1 = 0.485*bc*(nlength-idist)/idist*vol*scr*fac*(coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))/nlength;
%             else
                dforce1 = bc*(nlength-idist)/idist*vol*scr*fac*(coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1))/nlength;
%             end
            pforce(i,1) = pforce(i,1)+dforce1;
          
        end
        
    end
    
    for i=1:ndivx
        %Calculate the acceleration of material point i
%         if i > 500
%              acc(i,1) = (pforce(i,1)+bforce(i,1))/(dens*1.083);
%         else
            acc(i,1) = (pforce(i,1)+bforce(i,1))/dens;
            
%         end
        %acc(i,1) = (pforce(i,1)+bforce(i,1))/dens;
        
        %Calculate the velocity of material point i
        %by integrating the acceleration of material point i
        vel(i,1) = vel(i,1) + acc(i,1)*dt;
        
        %Calculate the displacement of material point i
        %by integrating the velocity of material point i
        disp(i,1) = disp(i,1)+vel(i,1)*dt;
    end
    
    %Store the displacement and time information for the material point at
    %the center of the bar
    pddisp(tt,1) = disp(400,1);
    pddisp2(tt,1) = disp(600,1);
    pdtime(tt,1) = ctime;
    
    %Calculate the analytical displacement solution of the material point
    %at the center
%     for nrao=0:ntotrao
%         andisp(tt,1) = andisp(tt,1) + ((-1)^nrao)/((2*nrao+1)^2)*sin((2*nrao+1)*pi*coord(600,1)/2)*cos((2*nrao+1)*pi*sqrt(10)*cwave*ctime/2);
%     end
%     andisp(tt,1) =  8*0.001*1/(pi^2)*andisp(tt,1);
end

%printing results
writeoutput=fopen('D:\Personal Files\SRF-IITH\Peridynamics\BenchmarkProblems\MATLAB_codes\ProjectCodes\coord.txt','w');
for i=1:nt
    %fprintf(writeoutput,'%d\t%d\t%d\n',pdtime(i,1),pddisp(i,1),andisp(i,1));
    fprintf(writeoutput,'%d\t%d\t%d\n',pdtime(i,1),pddisp(i,1),pddisp2(i,1));
end

% data = importdata('D:\Personal Files\SRF-IITH\Peridynamics\BenchmarkProblems\MATLAB_codes\ProjectCodes\coord.txt');
% figure;
% plot(data(:,1),data(:,2),data(:,1),data(:,3));
% %plot(data(:,1),data(:,2));
% xlabel('time');
% ylabel('Displacement');
% %legend('Peridynamic Solution','Analytical Solution');
% legend('At x=0.4(Steel)','At x=0.6(Brass)');
% %legend('Peridynamics Solution')
% title('Variation in Displacement with Time at x = 0.6 m  ');
toc
            
    
    



        
     
        
                
