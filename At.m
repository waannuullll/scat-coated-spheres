function Atot = At(m,n,mu,nu,k2d,tetalj,philj)
format long

%Bessel Functions
sphj= @(n,kr) sqrt(pi./(2*kr)).*besselj(n+0.5,kr);
sphy=@(n,kr) (sqrt(pi./(2*kr)).*(bessely(n+0.5,kr)));
sphh=@(n,kr) sphj(n,kr)+(1i.*sphy(n,kr));

%Translation coefficient function
Amn=@(m,n,mu,nu,phi)((((-1).^m).*((1i).^(nu+n)).*(gamma((2.*n)+1)./gamma(n+2)).*(gamma((2.*nu)+3)./gamma(nu+2)).*factorial(n+nu+m-mu))./(4.*n.*(gamma((2.*n)+(2.*nu)+1)./gamma(n+nu+1)).*factorial(n-m).*factorial(nu+mu)))*exp((1i).*(mu-m).*(philj+pi));
amn=@(n,nu,d,p,q,A,leg) ((-1).^q).*((n.*(n+1))+(nu.*(nu+1))-(p.*(p+1))).*A.*sphh(p,d).*leg;

%Gaunt coefficient function (Version1)
a1 = @(p,q) gamma(p+(2.*q)+0.5)./gamma(p+0.5);
a2 = @(m,n,mu,nu,q) gamma(n+nu+m-mu+1)./(((-1).^(2.*q)).*gamma(n+nu+m-mu-(2.*q)+1));
a3 = @(m,n,k) gamma(m+n+1)./(((-1).^(2.*k)).*gamma(m+n-2.*k+1));
a4 = @(mu,nu,q,k) gamma(nu-mu+1)./(((-1).^((2.*q)-(2.*k))).*gamma(nu-mu-(2.*q)+(2.*k)+1));
a5 = @(n,k) gamma(n+0.5)./(((-1).^k).*gamma(n-k+0.5));
a6 = @(nu,q,k)  gamma(nu+0.5)./(((-1).^(q-k)).*gamma(nu-q+k+0.5));
a7 = @(p,q,j) gamma(p+q-j+0.5)./(((-1).^(q-j)).*gamma(p+0.5));

%Addition Coefficients
M = n+nu;
qmax = (n+nu-abs(m-mu))/2;
for i =0:1:qmax
    pp(i+1) = M-(2.*i);
end
A = zeros(floor(qmax)+1,1);
A(1)=1;
for q=1:1:qmax
    a12= a1(pp(q+1),q)./a2(m,n,mu,nu,q);
    aa = 0;
    bb = 0;
    for k=0:1:q
        aa = aa+((a3(m,n,k).*a4(mu,nu,q,k))./(factorial(k).*factorial(q-k).*a5(n,k).*a6(nu,q,k)));
    end
    for j=0:1:q-1
        bb = bb +(((a7(pp(q+1),q,j))./factorial(q-j)).*A(j+1));
    end
    aq=(a12.*aa)-bb;
    A(q+1) = aq; 
end

asum = 0;
for q=0:1:qmax
    leg = legendrePP(pp(q+1),mu-m,cos(tetalj));
    asum = asum+ amn(n,nu,k2d,pp(q+1),q,A(q+1),leg);
end

Atot = Amn(m,n,mu,nu,philj).*asum;
