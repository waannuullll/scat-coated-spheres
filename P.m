function y = P(m,n,x)
p = legendre(n,x);
% if m >= 0
%     y = p(m+1);
% else
%     y = (-1)^(-m)*factorial(n+m)*p(-m+1)/factorial(n-m);
% end

if m >= 0
    y = p(m+1)/(-1)^m;
else
    y = factorial(n+m)*p(-m+1)/factorial(n-m);
end
end