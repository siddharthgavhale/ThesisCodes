
%%Initialization: particle (P), substrate(S) & vapor(V)

function [P,S,V,area] = C_initialization(initial_shape,N,dx,X,Y)
    
P=zeros(N,N); S=zeros(N,N);

%substrate
 for i=1:N
        for j=1:N
            if Y(i,j) <= 0
                 S(i,j) = 1;
            end
        end
 end

 %particle: semi-circle
 if initial_shape == 1
     radius=2;
     for i=1:N
        for j=1:N
            if (X(i,j)^2 + Y(i,j)^2 - radius^2 <= 0) && Y(i,j) >= 0
                P(i,j) = 1;
            end
        end
     end
     
     %Analytical area of initial semicircle
      area=round(0.5*radius*radius*pi/dx^2);%area of half circle 
 end
 
 %particle: square
  if initial_shape == 2
     side=2.5;
     for i=1:N
        for j=1:N
            if abs(X(i,j)) <= side/2 && abs(Y(i,j)) <= side && Y(i,j) >= 0
                P(i,j) = 1;
            end
        end
     end
     
     %Analytical area of initial semicircle
      area=round(side*side/dx^2);%area of half circle 
  end
  
 %particle: rectangle
  if initial_shape ==  3
     length = 2*pi;
     breadth = pi/2;
     for i=1:N
        for j=1:N
            if abs(X(i,j)) <= length/2 && abs(Y(i,j)) <= breadth && Y(i,j) >= 0
                P(i,j) = 1;
            end
        end
     end
     
     %Analytical area of initial semicircle
      area=round(length*breadth/dx^2);%area of half circle 
  end
   
  %particle: half astroid
  if initial_shape ==  4
      a = pi; %length of axis
     for i=1:N
        for j=1:N
            if (X(i,j)^2 + Y(i,j)^2 - a^2)^3 + 27*(X(i,j)*Y(i,j)*a)^2 <=0 && Y(i,j) >= 0
                P(i,j) = 1;
            end
        end
     end
     
     %Analytical area of initial semicircle
      area=round((3*pi*a^2/16)/dx^2);%area of half circle 
  end
 
  % (X(i,j)/2)^2 + (Y(i,j)/1.5)^2 -1 <= 0
  
  
  
   %particle: half ellipse
  if initial_shape ==  5
      a=3; %major axis 'a'
      b=1.5; %major axis 'b'
     for i=1:N
        for j=1:N
            if (X(i,j)/a)^2 + (Y(i,j)/b)^2 -1 <= 0 && Y(i,j) >= 0
                P(i,j) = 1;
            end
        end
     end
     
     %Analytical area of initial semicircle
      area=round((pi*a*b/2)/dx^2);%area of half circle 
  end
  
  
%   flower shape
%  %sqrt(X(i,j)^2 + Y(i,j)^2)-sin(5*atan2(Y(i,j),X(i,j))) <= 1.69 && Y(i,j) >= 0
    
 

V=ones(N,N)-(P+S);



end