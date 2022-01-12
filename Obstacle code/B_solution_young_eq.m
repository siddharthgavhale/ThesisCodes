

% find solution of anisotropic young equaion (Y_eq)
% Find parametric form of Wulff shape (Without inflation of area) (X_wulff & Y_wulff)
% Convert Wulff shape to  close curve by adding substrate  part in it
% Find area of Wulff shape

   
% m fold= 1, crystalline=2, ellipse = 3, ellipse with without smooth end =4,
function [ays1,ays2,xw,yw,Wint_area] = B_solution_young_eq(anisotropy,beta,phi,m,th,gammaSP,gammaSV)
    
    %Anisotropic young equaiton for given anisotropy
    if anisotropy == 1
       y_eq = @(th) -(-beta*m*sin(m*(phi + th))).*sin(th) + ...
          (beta*cos(m*(phi + th)) + 1).*cos(th) + gammaSP - gammaSV;

    %Exact solution/ analytical contact angle
    ays1=fzero(y_eq, -pi);
    ays2=fzero(y_eq, pi);
    
    %area of Wulff shape (Without triangle, which has vertices at Wulff point, cont_x_l & cont_x_r)
    %Actual formula as follows:
    % Wint_area=int(((G+ddG)*G),ani_young_sol_1,ani_young_sol_2); 
    % Wint_area=double(Wint_area*0.5)

    %To avoid symbolic toolbox, using following adjustment
    int_lim_1 = ays1 - m*((sin(2*m*(phi + ays1))*beta^2)/4 + sin(m*(phi + ays1))*beta) +...
        ((sin(2*m*(phi + ays1))*beta^2)/4 + 2*sin(m*(phi + ays1))*beta)/m + ...
        (beta^2*ays1)/2 - (beta^2*m^2*ays1)/2;

    int_lim_2 = ays2 - m*((sin(2*m*(phi + ays2))*beta^2)/4 + sin(m*(phi + ays2))*beta) +...
        ((sin(2*m*(phi + ays2))*beta^2)/4 + 2*sin(m*(phi + ays2))*beta)/m + ...
        (beta^2*ays2)/2 - (beta^2*m^2*ays2)/2;

    Winterbottom_area_pre = (int_lim_2-int_lim_1)/2;

    %find contact angle (on x-Axis)
    cont_x_l = double(subs((-(1+ beta*cos(m*(th+phi))).*sin(th)) - ...
             ((-beta*m*sin(m*(th+phi))).*cos(th)), th, ays1));
         
    cont_x_r = double(subs((-(1+ beta*cos(m*(th+phi))).*sin(th)) - ...
             ((-beta*m*sin(m*(th+phi))).*cos(th)), th, ays2));
         
    %area of triangle subtract or add
    Tri_area=abs((cont_x_r-cont_x_l)*(gammaSP-gammaSV))/2;

    if (gammaSP-gammaSV) <= 0
        Wint_area = Winterbottom_area_pre-Tri_area; %disp('sub')
    else
        Wint_area = Winterbottom_area_pre+Tri_area; %disp('add')
    end
    
    % Find parametric form of Wulff shape (Without inflation of area) (X_wulff & Y_wulff)  
        xw = double(subs((-(1+ beta*cos(m*(th+phi))).*sin(th)) - ...
             ((-beta*m*sin(m*(th+phi))).*cos(th)), th, ays1:0.01:ays2));
    
        yw = double(subs(-(-beta*m*sin(m*(phi + th))).*sin(th) + ...
               (beta*cos(m*(phi + th)) + 1).*cos(th) +...
                      gammaSP - gammaSV,th,ays1:0.01:ays2));   

    % make it a close curve (by adding substrate region)
    sub=xw(length(xw)):(xw(1)-xw(length(xw)))/1000:xw(1); %take end points on x axis
    xw=[xw sub]; %add that line
    yw=[yw zeros(1,length(sub))]; % same as xw, do with yw
    end
    
    %radian to degree
    ays1=rad2deg(ays1);
    ays2=rad2deg(ays2);
      
 end