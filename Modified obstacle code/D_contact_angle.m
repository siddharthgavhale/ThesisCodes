
%Determine contact angle 

function [CA_L, CA_R]=D_contact_angle...
                    (CP_upper_range,CP_lower_range,M,delta,...
                          dx,X,Y,ays1,ays2,sub_level)


cs=round(CP_lower_range/dx); %exclude first THIS many points
%if number of points is less than 2 then adjust it to 3
if cs < 2
    cs=3;
end

w=sub_level;

%Seting 1 by 1 matrix to avoid MAtlab warning 
%(If we don't define then it wont change anything)
x_cont_L=zeros(1,1); y_cont_L=zeros(1,1);
x_cont_R=zeros(1,1); y_cont_R=zeros(1,1);


%Counting parameter (keep it local)
condition=0; i=1;

    while CP_upper_range > condition
        % check points on diffused matrix at thresholding height (delta)
        temp= find(diff((M(w+cs,:)) <= delta));

        %left contact angle (using interpolatin, finding exact intterface, 
        %It will find interface even if particle is not exactly on grid point)
        dist_excl_l = dx*(delta - M(w+cs,temp(1)+1))/(M(w+cs,temp(1)) - M(w+cs,temp(1)+1));
        x_cont_L(i)=X(w+cs,temp(1)+1) - dist_excl_l;
        y_cont_L(i)=Y(w+cs,temp(1));

        %right contact angle
        dist_excl_l = dx*(delta - M(w+cs,temp(2)+1))/(M(w+cs,temp(2)) - M(w+cs,temp(2)+1));
        x_cont_R(i)=X(w+cs,temp(2)+1) - dist_excl_l;
        y_cont_R(i)=Y(w+cs,temp(2));

        condition= y_cont_R(i); % did we cross the uper limit i.e., CP_upper_range?
        i=i+1; cs=cs+1;

    end
%     figure()
%     contour(X,Y,M,[delta,delta],'b','LineWidth',2); hold on
%     plot(x_cont_R,y_cont_R,'-or','LineWidth',2)
%     plot(x_cont_L,y_cont_L,'-ok','LineWidth',2)
%     axis equal, grid on, xlabel('X'), 
%     axis([-3.5 3.5 0 3.5]) 

%Line regression  %Using polyfit command to find  contact  angle

%left contact angle (numerical)
p1 = polyfit(x_cont_L,y_cont_L,1);
CA_L=atand(p1(1));

%atan gives value between [-90,90], therefore adjusting it
% if abs(ays2) > 90
%     CA_L = CA_L+ sign(ays2)* 180;
% end


% %Right contact
 p2 = polyfit(x_cont_R,y_cont_R,1);
 CA_R=atand(p2(1));
  
%  %atan gives value between [-90,90], therefore adjusting it
% if abs(ays1) > 90
%     CA_R = CA_R + sign(ays1)* 180;
% end
 
end
