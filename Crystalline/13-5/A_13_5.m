
%plotting for the crystalline case



tic;
clc, clear, close all;

    
gpi=13;
dti=5;
       
% Initial setup for Numerical result    
N=2^gpi; L=10; dx = L/N;
x = -L/2:dx:L/2-dx; y = x; [X, Y]= meshgrid(x,y);

%analytical solution in sign distance form 
A = 2 - abs(X+Y) - abs(X-Y);

1
%Bon
load(sprintf('%s_%d %d','Bon_U_gpi_dti',gpi,dti),'U'); 
contour(X,Y,U,[0.5 0.5],'--r','DisplayName','a','linewidth',4.5), hold on
2
%EE eps=0.01
epsi=0.01;
load(sprintf('%s_%d %d %d','EE_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],':b','DisplayName','b','linewidth',4.5),
3
%EE eps=0.05
epsi=0.05;
load(sprintf('%s_%d %d %d','EE_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],'-.g','DisplayName','c','linewidth',4.5),
4
%EE eps=0.1
epsi=0.1;
load(sprintf('%s_%d %d %d','EE_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],'--y','DisplayName','d','linewidth',4.5),
5
%EJZ FD
save(sprintf('%s_%d %d','EJZ_FD_U_gpi_dti',gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],':c','DisplayName','e','linewidth',4.5),
6
%EJZ PD
epsi=0.01;
save(sprintf('%s_%d %d %d','EJZ_PD_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],'-.m','DisplayName','f','linewidth',4.5),
7
%Analytical solution
contour(X,Y,A,[0 0],'-k','DisplayName','g','linewidth',4.5), 
8

set(gca,'FontSize',55),  
axis equal, axis([0.75 1.1 0.75 1.1])
title('\fontsize{40} dx=0.0012 & \deltat = 0.078');

 disp('small plot')

text(0.98,0.98,'B','Fontsize',40)
text(0.795,0.9,'A \rightarrow','Fontsize',40)
text(0.92,0.96,'C \rightarrow','Fontsize',40)
text(0.83,0.94,'D, E, F \rightarrow','Fontsize',40)
text(1,0.96,'\leftarrow G','Fontsize',40)

%small box
axes('position',[.57 .68 .24 .24]),



%Analytical solution
1
%Bon
load(sprintf('%s_%d %d','Bon_U_gpi_dti',gpi,dti),'U'); 
contour(X,Y,U,[0.5 0.5],'--r','DisplayName','a','linewidth',3), hold on
2
%EE eps=0.01
epsi=0.01;
load(sprintf('%s_%d %d %d','EE_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],':b','DisplayName','b','linewidth',3),
3
%EE eps=0.05
epsi=0.05;
load(sprintf('%s_%d %d %d','EE_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],'-.g','DisplayName','c','linewidth',3),
4
%EE eps=0.1
epsi=0.1;
load(sprintf('%s_%d %d %d','EE_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],'--y','DisplayName','d','linewidth',3),
5
%EJZ FD
save(sprintf('%s_%d %d','EJZ_FD_U_gpi_dti',gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],':c','DisplayName','e','linewidth',3),
6
%EJZ PD
epsi=0.01;
save(sprintf('%s_%d %d %d','EJZ_PD_eps_gpi_dti',epsi*100,gpi,dti),'U');
contour(X,Y,U,[0.5 0.5],'-.m','DisplayName','f','linewidth',3),
7
contour(X,Y,A,[0 0],'-k','DisplayName','g','linewidth',3), 


plot([0.6 1.25 1.25 0.6 0.6],[0.6 0.6 1.25 1.25 0.6],'k','Linewidth',2);
 

axis equal, 
set(gca,'Visible','off')
axis([-1.5 1.5 -1.5 1.5]); 
title('\fontsize{40} dx=0.0012 & \deltat=0.078');
hold off


toc




 
 
 
 
 
 
 






