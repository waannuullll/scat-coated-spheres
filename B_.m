function Bt = B(m,n,mu,nu,kd,theta0,phi0)

J = @(n,x) sqrt(pi/(2*x))*besselj(n+0.5,x); %spherical bessel
h = @(n,x) sqrt(pi/(2*x))*(besselj(n+0.5,x) + 1i*bessely(n+0.5,x));; %spherical Hankel
qmax = min([n,nu,0.5*(n+nu-abs(m-mu))]);
Qmax = min([n+1,nu,0.5*(n+1+nu-abs(m-mu))]);
% if m == -mu & n == nu
%     Bt = 0;
% elseif m == -n & mu == -nu
%     Bt = 0;
% elseif m == n & mu == nu
%     Bt = 0;
% else
    %calculate Gaunt coefficient 0-->1
    %calculate the denominator of a2q and a3q
    %a2(1) = (pochhammer(n+1+1,n+1)*pochhammer(nu+1,nu)*factorial(n+1+nu-(-m-1)-(mu+1)))...
    %    /(pochhammer(n+1+nu+1,n+1+nu)*factorial(n+1-(-m-1))*factorial(nu-(mu+1)));
    %a3(1) = (pochhammer(n+1+1,n+1)*pochhammer(nu+1,nu)*factorial(n+1+nu--m-mu))...
    %    /(pochhammer(n+1+nu+1,n+1+nu)*factorial(n+1--m)*factorial(nu-mu));
    a2(1) = 1; a3(1) = 1;
    %b0 = ((2*n+1)*(n+nu+m-mu+1))/((2*n+2*nu+1)*(n+m+1));
    %calculate the normalized a2q
    for q=2:Qmax+1
        p = n+1+nu-2*(q-1);
        n4 = n+1+nu-(-m-1)-(mu+1);
        S1 = 0;
        for k=1:q
            s1 = (poc(-m-1-(n+1),2*k-2)*poc(mu+1-nu,2*q-2*k))...
                /(factorial(k-1)*factorial(q-k)*poc(-(n+1)+0.5,k-1)*poc(-nu+0.5,q-k));
            S1 = S1 + s1;
        end
        S2 = 0;
        for j=1:q-1
            s2 = a2(j)*poc(-p-q+j+0.5,q-j)/factorial(q-j);
            S2 = S2 + s2;
        end
        a2(q) = poc(p+0.5,2*q-2)*S1/poc(-n4,2*q-2)-S2;
    end
    %calculate the normalized a3q
    for q=2:Qmax+1
        p = n+1+nu-2*(q-1);
        n4 = n+1+nu--m-mu;
        S1 = 0;
        for k=1:q
            s1 = (poc(-m-(n+1),2*k-2)*poc(mu-nu,2*q-2*k))...
                /(factorial(k-1)*factorial(q-k)*poc(-(n+1)+0.5,k-1)*poc(-nu+0.5,q-k));
            S1 = S1 + s1;
        end
        S2 = 0;
        for j=1:q-1
            s2 = a3(j)*poc(-p-q+j+0.5,q-j)/factorial(q-j);
            S2 = S2 + s2;
        end
        a3(q) = poc(p+0.5,2*q-2)*S1/poc(-n4,2*q-2)-S2;
    end
    %calculate B coefficient

    B1 = exp(1i*(mu-m)*(phi0))*((-1)^m*(1i)^(nu+n+1)*poc(n+2,n+1)*poc(nu+2,nu+1)*factorial(n+nu+m-mu+1))...
        /(4*n*(n+1)*(n+m+1)*poc(n+nu+2,n+nu+1)*factorial(n-m)*factorial(nu+mu));
    %B1 = (-1)^(m+n)*1i*a0*b0*((2*n+1)/(2*n*(n+1)))*exp(1i*(mu-m)*phi0);
    B2 = 0;
    for q=1:Qmax+1
        p = n+nu-2*(q-1);
        sB = (-1)^(q-1)*(2*(n+1)*(nu-mu)*a2(q)-(p*(p+3)-nu*(nu+1)-n*(n+3)-2*mu*(n+1))*a3(q))...
            *J(p+1,kd)*P(mu-m,p+1,cos(theta0));
        B2 = B2 + sB;
    end
    Bt = B1*B2;
end
