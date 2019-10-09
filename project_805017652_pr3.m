%% MAE M2O Final Project Problem 3
% Generates position, velocity, and turning angle data implicitly and
% explicitly for any number of nodes
%Catie Vitella


clc; clear all;
n = 3;

[timplicit, yimplicit, vimplicit] = part3implicit(n);

%plot
figure (1)
plot(timplicit,yimplicit,'b-');
hold on
plot(timplicit,vimplicit,'r-');
str = sprintf('Position and Veloity for Middle Node of %d -Node System, Implicit', n);
title(str);
legend('y(t)', 'v(t)');
xlabel('t');




    