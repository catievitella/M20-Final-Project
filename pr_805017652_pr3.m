%% MAE M20 Final Project Problem 3

%% M20 Final Project Problem 3
%Catie Vitella

%Plots the positions and velocities of three balls attached to an elastic
%beam falling in viscous fluid using explicit and implicit methods

clc; clear all;

%initialize
n = 3;
%make sure n is odd
if mod(n,2) == 0
    error('Number of nodes must be odd.');
end

u = 1000; %viscosity (Pa-s)git@github.com:catievitella/M20-Final-Project.git
pmetal = 7000; %density of metal (kg/m^3)
pfluid = 1000; %density of fluid (kg/m^3)
l = 0.10; %length of beam (m)
dl = l/(n-1); %length of partitions of beam
r0 = 0.001; %cross-sectional radius of beam (m)
E = 1*(10^9); %Young modulus (Pa)
g = 9.8; %acceleration due to gravity (m/s^2)
tfinal = 10;
dt = 10^(-5);
tsteps = tfinal/dt + 1;
R = createR(dl,n); %creates radius vector, m
W = (4/3).*pi.*(R^3)*(pmetal-pfluid)*g; %net weight of each sphere vector
m = (4/3)*pi*(R^3)*pmetal; %mass of each sphere vector, kg
%preallocating
qkx = zeros(n,1); %x-positions of balls at tk
qky = zeros(n,1); %y-positions of balls at tk
qkxp1 = zeros(n,1); %x-positions of balls at tk-1
qkyp1 = zeros(n,1); %y-positions of balls at tk-1
vkx = zeros(n,1); %x-velocities of balls at tk
vky = zeros(n,1); %y-velocities of balls at tk
dEs1k = zeros(4,1); %partial derivatives of stretching moment for left spring
dEs2k = zeros(4,1); %partial derivatives of stretching moment for right spring
dEb = zeros(6,1); %partial derivatives of bending moment for system
dEx = zeros(3,1); %partial derivatives of total elastic energy in x
dEy = zeros(3,1); %partial derivatives of total elastic energy in y

%define stiffnesses
EA = E*pi*(r0)^2; %stretching stiffness
EI = E*pi*((r0)^4)/4; %bending stiffness

for n = 1:tsteps
    
    %store
    x1(n) = qkx(1); y1(n) = qky(1); x2(n) = qkx(2); y2(n) = qky(2); x3(n) = qkx(3); y3(n) = qky(3);
    vx1(n) = vkx(1); vy1(n) = vky(1); vx2(n) = vkx(2); vy2(n) = vky(2); vx3(n) = vkx(3); vy3(n) = vky(3);
    t(n) = n*dt - dt;
    
    %update turning angle
    theta(n) = (((qkx(3)-qkx(2))^2)+((qkx(1)-qkx(2))^2) + ((qky(3)-qkx(2))^2)+ ((qky(1)-qkx(2))^2))/(2*sqrt((qkx(3)-qkx(2))^2)+((qky(3)-qkx(2))^2))*sqrt(((qkx(1)-qkx(2))^2)+((qky(1)-qkx(2))^2));
    
    %compute next position values
    qkxp1 = qkx + (((m.*vkx) - (dt.*dEx))./((m./dt)+(6*pi*u*R))); 
    vkxp1 = (qkxp1-qkx)./dt;
    qkyp1 = qky + (((m.*vky) - (dt.*dEy) - (dt.*W))./((m./dt)+(6*pi*u*R)));
    vkyp1 = (qkyp1-qky)./dt;
   
    %update
    qkx = qkxp1;
    vkx = vkxp1;
    qky = qkyp1;
    vky = vkyp1;
    
    %gradES
    dEs1k = gradEs(qkx(1), qky(1), qkx(2), qky(2), dl, EA);
    dEs2k = gradEs(qkx(2), qky(2), qkx(3), qky(3), dl, EA);
    dEb = gradEb(qkx(1), qky(1), qkx(2), qky(2), qkx(3), qky(3), dl, EI);
    
    %energy gradients
    dEx = [dEs1k(1) + dEb(1); dEs1k(3) + dEs2k(1) + dEb(3); dEs2k(3) + dEb(5)];
    dEy = [dEs1k(2) + dEb(2); dEs1k(4) + dEs2k(2) + dEb(4); dEs2k(4) + dEb(6)];

end

%plot
figure (1)
plot(t,y2,'b-');
title('y-position of middle node');
xlabel('t');
ylabel('y');

figure (2)
plot(t,vy2,'b-');
title('y-velocity of middle node');
xlabel('t');
ylabel('v');

figure (3)
plot(t,theta,'b-');
title('turning angle over time');
xlabel('t');
ylabel('theta');
