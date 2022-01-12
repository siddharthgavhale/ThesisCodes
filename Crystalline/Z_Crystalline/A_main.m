% for figure 4:
close all;clear;clc

%Yes=1; No=0;
plot1=1; %want 1st plot? Yes=1; No=0;
plot2=1; %want 2nd plot?
plot3=1; %want 3rd plot?
movie=0; %Want movie?


%
if Plot1 ==1
    gpi=13; dti=5; 

    A1_BBC(gpi,dti,movie);%BBC
    A2_EJZ_FD(gpi,dti,movie);%EJZ_FD
    A3_EJZ_PD(gpi,dti,movie);%EJZ_PD

    eps=0.1;
    A4_EE(gpi,dti,eps,movie);%EE

    eps=0.05;
    A4_EE(gpi,dti,eps,movie);%EE

    eps=0.01;
    A4_EE(gpi,dti,eps,movie);%EE
end



if Plot2== 1
    gpi=13; dti=7;

    A1_BBC(gpi,dti,movie);%BBC
    A2_EJZ_FD(gpi,dti,movie);%EJZ_FD
    A3_EJZ_PD(gpi,dti,movie);%EJZ_PD

    eps=0.1;
    A4_EE(gpi,dti,eps,movie);%EE

    eps=0.05;
    A4_EE(gpi,dti,eps,movie);%EE

    eps=0.01;
    A4_EE(gpi,dti,eps,movie);%EE
end



if Plot3 == 1
    gpi=13; 

    dti=9;
    A1_BBC(gpi,dti,movie);%BBC

    dti=11;
    A2_EJZ_FD(gpi,dti,movie);%EJZ_FD

    dti=6;
    A3_EJZ_PD(gpi,dti,movie);%EJZ_PD

    dti=0;
    eps=0.1;
    A4_EE(gpi,dti,eps,movie);%EE

    eps=0.05;
    A4_EE(gpi,dti,eps,movie);%EE

    eps=0.01;
    A4_EE(gpi,dti,eps,movie);%EE
end



