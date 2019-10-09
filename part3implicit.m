function [t,ymid,vmid] = part3implicit(n)

clc; 
%initialize
%n = 21; %number of nodes
if mod(n,2) == 0
    error ('Must be an odd number of nodes'); 
end
mid = n; %index of middle node
u = 1000; %viscosity (Pa-s)
l = 0.10; %total length of beam (m) 
dl = l/(n-1); %length of each partition of beam (m)
pmetal = 7000; %density of metal (kg/m^3)
pfluid = 1000; %density of fluid (kg/m^3)
r0 = 0.001; %cross-sectional radius of beam (m)
E = 1*(10^9); %Young modulus (Pa)
g = 9.8; %acceleration due to gravity (m/s^2)
tfinal = 50;
dt = 10^(-2);
tsteps = round(tfinal/dt + 1);
EA = E*pi*(r0)^2; %stretching stiffness
EI = E*pi*((r0)^4)/4; %bending stiffness
v = zeros((n*2),1);
ymid = zeros(tsteps,1);
vmid = zeros(tsteps,1); %velocity to be stored
tol = EI/l^2 * 0.001;
err = 1;

%create R (radius of each node, m)
R = zeros((n*2),1); 
R(1:(mid-1)) = dl/10;
R(mid:mid+1) = 0.025; 
R((mid+2):(n*2)) = dl/10;

%create W (net weight of each sphere, N)
W = zeros((n*2),1);
W(2:2:(n*2)) = (4/3)*pi.*(R(2:2:(n*2)).^3)*(pmetal-pfluid)*g; % evens (y-coordinates) have a weight force in that direction
W(1:2:(n*2)) = zeros(((n*2)/2),1); %no weight force in x-direction (evens are x-coordinates)

m = (4/3)*pi*(R.^3)*(pmetal-pfluid)*g; %mass vector for each node

%create q-vector (position of each node)
qold = zeros((n*2),1);
for k = 1:n
    qold((2*k)-1) = (k-1)*dl;
end
qnew = qold;

%set up videowriter
% vid = VideoWriter('implicitmovie2.avi'); %uncomment to create video
% vid.FrameRate = 100; %uncomment to create video
% open(vid); %uncomment to create video

for i = 1:tsteps
    %store
    ymid(i) = qold(mid+1);
    vmid(i) = v(mid+1);
    t(i) = i*dt - dt;
    
    qnew = qold;
    
    while err > tol

        dE = create_dE(n,qnew);
        hE = create_hE(n,qnew);
        
        %implicit equations
        f = m / dt .* ( (qnew - qold) / dt - v ) + 6*pi*u*R.*(qnew-qold)/dt + W + dE;
        
        %create Jacobian
        J = zeros(n*2,n*2);
            for c=1:n*2
                J(c,c) = m(c)/dt^2 + 6*pi*u*R(c)/dt;
            end
        J = J + hE;
        
        dpos = J \ f;

        qnew = qnew - dpos;
        err = sum(abs(f));  

    end
    
    if mod(i,6) == 0
     plot(qnew(1:2:(n*2)), qnew(2:2:(n*2)), 'bo-');
     ylim([-.32,0]);
     drawnow
    end
%     saveas(gca, 'frame1.jpg'); %uncomment to create video
%     frame = imread('frame1.jpg'); %uncomment to create video
%     writeVideo(vid, frame); %uncomment to create video
    
    v = (qnew-qold)/dt;
    qold = qnew;
    err = 1;
end
% close(vid); %uncomment to create video
    
end




        
        
        
    