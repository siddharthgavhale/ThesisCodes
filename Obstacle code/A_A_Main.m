

%It will generate error plot for each kernel and evolution plot using BBC
%kernel and at the end contact angle(left)  for all the kernels. 
%plots will be saved automatically
%keep the  gpi and dti range of all the "A_***"  codes as follows
%range of grid points %gpiSV=9; gpiEV=12;
 %range of diffusion time %dtiSV=0; dtiEV=7;
clear;close all;clc

%BBC
A_BBC

%EE
 A_EE

%EJZ FD
A_EJZFD

%EJZ PD
A_EJZPD


%contact angle plot
%line type
linS = {'-',':','-.','--','--*',':o','-.v','--^'};

%Bon
load('CA_BBC.mat');
plot(DTT,CA_BBC(1,:),'linestyle',linS{1},'Linewidth',7.5), hold on

%EE
load('CA_EE.mat');
plot(DTT,CA_EE(1,:),'linestyle',linS{2},'Linewidth',7.5), hold on

%EJZ FD
load('CA_EJZ_FD.mat');
plot(DTT,CA_EJZ_FD(1,:),'linestyle',linS{3},'Linewidth',7.5), hold on

%EJZ PD
load('CA_EJZ_PD.mat');
plot(DTT,CA_EJZ_PD(1,:),'linestyle',linS{4},'Linewidth',7.5), hold on

 text(0.004,120,'\leftarrow BBC','Fontsize',40)
  text(0.019,95,'\leftarrow EE (\epsilon = 0.01)','Fontsize',40)
   text(0.05,105,'EJZ(positive physical) \rightarrow','Fontsize',40)
     text(0.042,118,'EJZ(positive Fourier) \rightarrow','Fontsize',40)
        text(0.095,130,'\downarrow','Fontsize',40)
          text(0.085,133,'Analytical','Fontsize',40)  
 

hold on, set(gca,'FontSize',60)
xlabel('\deltat','FontSize',55),ylabel('Contact angle','FontSize',55),

set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
line([0,0.25],[127.37,127.37],'DisplayName','Analytical','Color','k','LineWidth',5);
xlim([0 0.125]); ylim([ 87.37 147.37]);

%maximise the plot and save
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
saveas(gcf,'Fig6_5.png')



% title('\fontsize{40} Contact angle');

hold off





