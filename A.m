function At = A(m,n,mu,nu,kd,theta0,phi0)
format long

j = @(n,x) sqrt(pi/(2*x))*besselj(n+0.5,x); %spherical bessel
h = @(n,x) sqrt(pi/(2*x))*(besselj(n+0.5,x) + 1i*bessely(n+0.5,x)); %spherical Hankel
qmax = min([n,nu,0.5*(n+nu-abs(m-mu))]);
Qmax = min([n+1,nu,0.5*(n+1+nu-abs(m-mu))]);
%calculate Gaunt coefficient 0-->1
%calculate the denominator of a1q
%a1(1) = (pochhammer(n+1,n)*pochhammer(nu+1,nu)*factorial(n+nu--m-mu))...
%    /(pochhammer(n+nu+1,n+nu)*factorial(n--m)*factorial(nu-mu));
a1(1) = 1;
a0 = a1(1);
%calculate the normalized a1q
for q=2:Qmax+1
    p = n+nu-2*(q-1);
    n4 = n+nu--m-mu;
    S1 = 0;
    for k=1:q
        s1 = (poc(-m-n,2*k-2)*poc(mu-nu,2*q-2*k))...
            /(factorial(k-1)*factorial(q-k)*poc(-n+0.5,k-1)*poc(-nu+0.5,q-k));
        S1 = S1 + s1;
    end
    S2 = 0;
    for j=1:q-1
        s2 = a1(j)*poc(-p-q+j+0.5,q-j)/factorial(q-j);
        S2 = S2 + s2;
    end
    a1(q) = poc(p+0.5,2*q-2)*S1/poc(-n4,2*q-2)-S2;
end
%calculate A coefficient
A1 = exp(1i*(mu-m)*(phi0))*((-1)^m*(1i)^(nu+n)*poc(n+2,n-1)*poc(nu+2,nu+1)*factorial(n+nu+m-mu))...
     /(4*n*poc(n+nu+1,n+nu)*factorial(n-m)*factorial(nu+mu));
%A1 = (-1)^(m+n)*a0*((2*n+1)/(2*n*(n+1)))*exp(1i*(mu-m)*phi0);
A2 = 0;
for q=1:qmax+1
    p = n+nu-2*(q-1);
    sA = (-1)^(q-1)*(n*(n+1)+nu*(nu+1)-p*(p+1))*a1(q)*h(p,kd)*P(mu-m,p,cos(theta0));
    A2 = A2 + sA;
end
At = A1*A2;
