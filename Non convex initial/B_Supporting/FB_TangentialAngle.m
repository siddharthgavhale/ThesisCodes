
function [nu , nud] = FB_TangentialAngle(tv,N)

%Tangle angle 'nu'
 nu = zeros(1,N+3);
    if tv(2,1) < 0
        nu(1) = 2*pi - acos(tv(1,1));
    else
        nu(1) = acos(tv(1,1));
    end
    
    for i = 1:N
        II = dot(tv(:,i),tv(:,i+1));
        DD = det([tv(:,i)';tv(:,i+1)']);
        if II > 0
            nu(i+1) = nu(i) + asin(DD);  
        elseif DD > 0
            nu(i+1) = nu(i) + acos(II);
        else
            nu(i+1) = nu(i) - acos(II);
        end
    end
    nu(N+3) = nu(1) - (nu(N+1) - nu(N));
    nu(N+2) = nu(N+1) + (nu(2) - nu(1));
    
    nud = zeros(1,N+2);
    nud(1:N+2) = (nu(1:N+2) + nu(2:N+3))/2;
end