function [d, f] = gTRandomiser(dRange, fRange)
    % Find random values for d and f in specified ranges
    
    np = length(dRange);
    d = rand(np,1).*(dRange(:,2)-dRange(:,1)) + dRange(:,1);
    f = zeros(np,1);
    f(1) = rand()*(fRange(1,2)-fRange(1,1)) + fRange(1,1);
    
    for i=2:np-1
        ul = 100-fRange(np,1)-sum(f(1:i-1));    % upper limit for f(2)/f_inter including constraint for minimum f(3)/f_fast
        ll = fRange(i,1);
        if ul > fRange(i,2)                     % upper limit for f(2) reduced if f(1) too small
           ul = fRange(i,2); 
        end             
        if ll < 100 - f(1) - fRange(i+1,2)      % lower limit for f(2) including constraintr for maximum f(3)
            ll = 100 - f(1) - fRange(i+1,2);
        end
        f(i) = rand()*(ul-ll) + ll;
    end
    f(np) = 100 - sum(f(1:np-1)); % sum(f_i) != 1
    
    d=d';
    f=f';
end