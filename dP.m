function y = dP(m,n,x)
    Y = -(n+1)*x*P(m,n,x) + (n-m+1)*P(m,n+1,x);
    X = Y/(x^2-1);
    y = -X*sqrt(1-x^2);
end