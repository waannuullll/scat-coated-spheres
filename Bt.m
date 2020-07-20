function Btot = Bt(m,n,mu,nu,k2d,tetalj,philj)
format long

%Bessel Functions
sphj=@(n,kr) sqrt(pi./(2*kr)).*besselj(n+0.5,kr);
sphy=@(n,kr) (sqrt(pi./(2*kr)).*(bessely(n+0.5,kr)));
sphh=@(n,kr) sphj(n,kr)+(1i.*sphy(n,kr));

%Translation coefficient function
Bmn = @(m,n,mu,nu,phi) ((((-1).^m).*((1i).^(nu+n+1)).*poc(n+2,n+1).*poc(nu+2,nu+1).*factorial(n+nu+m-mu+1))./(4.*n.*(n+1).*(n+m+1).*poc(n+nu+2,n+nu+1).*factorial(n-m).*factorial(nu+mu))).*exp((1i).*(mu-m).*(philj+pi));
bmn = @(n,mu,nu,d,p,q,E,G,leg) ((-1).^q).*((2.*(n+1).*(nu-mu).*E)-(((p.*(p+3))-(nu.*(nu+1))-(n.*(n+3))-(2.*mu.*(n+1))).*G)).*sphh(p+1,d).*leg;
if m == -mu & n == nu
    Btot = 0;
elseif m == -n & mu == -nu
    Btot = 0;
elseif m == n & mu == nu
    Btot = 0;
else
a1 = @(p,q) poc(p+1.5,2*q);
    a2 = @(m,n,mu,nu,q) poc(-n-nu-m+mu-1,2*q);
    a3 = @(m,n,k) poc(-m-n-2,2*k);
    a4 = @(mu,nu,q,k) poc(mu-nu+1,(2*q)-(2*k));
    a5 = @(n,k) poc(-n+0.5-1,k);
    a6 = @(nu,q,k)  poc(-nu+0.5,q-k);
    a7 = @(p,q,j) poc(-p-q+j+0.5-1,q-j);

    M = n+nu;
    qmax = (n+nu-abs(m-mu))/2;
    Qmax=(n+nu+1-abs(m-mu))/2;
    for i =0:1:floor(Qmax)
        pp(i+1) = M-(2.*i);
    end

    E = zeros(floor(Qmax)+1,1);
    E(1)=1;
    G = zeros(floor(Qmax)+1,1);
    G(1)=1;

    for q=0:1:floor(Qmax)
        a12= a1(pp(q+1),q)./a2(m,n,mu,nu,q);
        aa = 0;
        bb = 0;
        for k=0:1:q
            aa = aa+((a3(m,n,k).*a4(mu,nu,q,k))./(factorial(k).*factorial(q-k).*a5(n,k).*a6(nu,q,k)));
        end
        for j=0:1:q-1
            bb = bb +(((a7(pp(q+1),q,j))./factorial(q-j)).*E(j+1));
        end
        aq=(a12.*aa)-bb;
        E(q+1) = aq;
    end
    %Gaunt coefficient function (Version1)
    a11 = @(p,q) poc(p+1.5,2*q);
    a21 = @(m,n,mu,nu,q) poc(-n-nu-m+mu-1,2*q);
    a31 = @(m,n,k) poc(-m-n-1,2*k);
    a41 = @(mu,nu,q,k) poc(mu-nu,(2*q)-(2*k));
    a51 = @(n,k) poc(-n+0.5-1,k);
    a61 = @(nu,q,k)  poc(-nu+0.5,q-k);
    a71 = @(p,q,j) poc(-p-q+j+0.5-1,q-j);

    for q=0:1:floor(Qmax)
        a12= a11(pp(q+1),q)./a21(m,n,mu,nu,q);
        aa = 0;
        bb = 0;
        for k=0:1:q
            aa = aa+((a31(m,n,k).*a41(mu,nu,q,k))./(factorial(k).*factorial(q-k).*a51(n,k).*a61(nu,q,k)));
        end
        for j=0:1:q-1
            bb = bb +(((a71(pp(q+1),q,j))./factorial(q-j)).*G(j+1));
        end
        aq=(a12.*aa)-bb;
        G(q+1) = aq;
    end

    bsum = 0;
    for q=0:1:floor(Qmax);
        leg = legendrePP(pp(q+1)+1,mu-m,cos(tetalj));
        bsum = bsum+ bmn(n,mu,nu,k2d,pp(q+1),q,E(q+1),G(q+1),leg);
    end

    Bmn(m,n,mu,nu,philj);
    Btot =Bmn(m,n,mu,nu,philj).*bsum;
end
