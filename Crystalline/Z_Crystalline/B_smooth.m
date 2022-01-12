
%smoothening function

function A = B_smooth(A,N)

AC = A; % a copy of A is needed

for i=2:N-1
    for j=2:N-1
        if sign(A(i,j)) == sign(A(i,j+1)) == sign(A(i+1,j)) == sign(A(i,j-1)) ... 
                == sign(A(i-1,j)) == sign(A(i-1,j-1)) == sign(A(i-1,j+1)) ...
                == sign(A(i+1,j-1)) == sign(A(i+1,j+1)) 
            A(i,j)= 0.5*sign(A(i,j));
        else % compAte the area ratio of positive phase on 4 grid cells centered at (i,j)
            
    b=10^-7;   
w1=AC(i-1,j-1); w2=AC(i-1,j); w3=AC(i,j-1);w4=AC(i,j);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=areacell;

w1=AC(i-1,j); w2=AC(i-1,j+1); w3=AC(i,j);w4=AC(i,j+1);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=Areacell+areacell;

w1=AC(i,j-1); w2=AC(i,j); w3=AC(i+1,j-1);w4=AC(i+1,j);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=Areacell+areacell;


w1=AC(i,j); w2=AC(i,j+1); w3=AC(i+1,j);w4=AC(i+1,j+1);
if (w1 > 0)
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                areacell = 1.0;
            else
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            end
        else
            if (w4 > 0)
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            end
        else
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
                
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end
            else
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    end
else
    if (w2 > 0)
        if (w3 > 0)
            if (w4 > 0)
                u1 = w3; u2 = w1; u3 = w4; u4 = w2;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = 1.0 - (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = 1.0 - ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end

            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
                
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = areacell + (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = areacell + ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w1; u2 = w2; u3 = w3; u4 = w4;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        end
    else
        if (w3 > 0)
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                   if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                       areacell= (u2*(u4-u3) + u4*(u2-u1))/(2.0*(u4-u3)*(u2-u1)); 
                   else
                       areacell = (u2-u4+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                        log((u2-u1)/(u4-u3)))/(-u1+u2+u3-u4);
                   end
            else
                u1 = w4; u2 = w3; u3 = w2; u4 = w1;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            end
        else
            if (w4 > 0)
                u1 = w2; u2 = w4; u3 = w1; u4 = w3;
                  if -b < -u1+u2+u3-u4  && -u1+u2+u3-u4 < b
                         areacell = (u2*u2)/(2.0*(u2-u1)*(u2-u4));
                  else
                         areacell = ((u2+(u2*u3-u1*u4)/(-u1+u2+u3-u4)* ...
                           log((u1-u2)*(u4-u2)/(u1*u4-u2*u3)))/(-u1+u2+u3-u4));
                  end 
            else
                areacell = 0.0;
            end
        end
    end
end

Areacell=Areacell+areacell;            
            A(i,j) = Areacell/4.0 -0.5;
        end
    end
end
A=A+0.5;

%AdjAsting boAndary
A(1,:) = 0; A(N,:) = 0; 
A(:,1) = 0; A(:,N) = 0;

end
