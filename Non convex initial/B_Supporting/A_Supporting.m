
%this code will generate all the necessary f iles for S shape plot at time
%0.25, 0.75 and 1.25

clear; clc; close all;

%front trackinf code: It will generate two files: 
%1: solution at each time 
%2: solution at specific time ie., T=0.25, 0.75 and 1.25
%don't change the number of points representing front i.e., N=100;

disp('Front tracking code in progress!')
F_front_tracking

%Bonnetier kernel
%this code have option to create movie for both front tracking and bbc
disp('BBC code in progress!')
Bon


%EJZ PD kernel
%this code have option to create movie for both front tracking and EJZ_PD
disp('EJZ PD code in progress!')

EJZ_PD


%EJZ_FD kernel
%this code have option to create movie for both front tracking and EJZ_FD
%  this code have adjusted the time scalling because this kernel does not
%  represent exact mobility
disp('EJZ_FD code in progress!')

EJZ_FD
 
 
 
 
 
 
 






