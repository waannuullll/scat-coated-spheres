function [poch] = poc(x,n)

if x<0 
    poch = ((-1).^n).*(gamma(abs(x)+1)./gamma(abs(x)-n+1));
else if x == 0
    if n ~= 0
        poch = 0;
    else
        poch = 1;
    end
else
    poch = gamma(x+n)./gamma(x);
end
end