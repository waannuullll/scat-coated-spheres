% --- COATED MULTISPHERE SCATTERING --- (x polarized)

% This program can be completed after going through two stages. This 
% program was created by:
% - (May 2019) Misael Natanael: made multisphere scattering program 
% benchmarked with Mühlig's FDTD.
% - (May 2020) Ikhwanul Muslimin: made coated multisphere scattering 
% program based on Misael's program benchmarked with Crut's experiment.
% Both are from Department of Physics, Faculty of Mathematics and Natural 
% Sciences, Bandung Institute of Technology.

% There are several variables that must be changed according to the 
% configuration being reviewed, namely:
% 1. Amplitude of the incoming field (Line 45)
% 2. Step in wavelength (Line 46)
% 3. Number of spheres (Line 47)
% 4. Number of layers each sphere (Line 48)
% 5. Number of Bessel mode (Line 49). Normally converge when N = 4 if the 
% configuration is simple.
% 6. Choose material reference (Line 53-69).
% 7. Range of lambda (Line 71). Begin with lambdamin and end with lambdamax.
% 8. Background's refractive index (Line 82).
% 9. Radius of each sphere (Line 83).
% 10. Distance between the spheres (Line 84).
% 11. Position of each sphere (Line 87-90).
% 12. Size parameter (Line 94). Just add/remove k*R(j) based on your number
% of spheres.
% 13. Choose material to be used (Line 96-118).
% 14. Choose multilayer parameter (radius, size parameter, relative
% refractive index) (Line 123-139).
% 15. What to plot (Line 370-439).
% Notice: DON'T CHANGE OTHER PART
% To comment use CTRL+R to uncomment use CTRL+T

% When you first download this program, it should work well. I gave an 
% example of a configuration of two gold spheres (R1 = 38 nm, Johnson 
% Christy material reference) coated with a dielectric (R2 = 40 nm, 
% constant parameter) that were separated by 5 nm positioned on 
% the x-axis in water.

clearvars; clc; tic; progressbar;
%% Definition of General Variables
c = 3E+08;                                                                  % speed of light 
mu0 = 4*pi*1E-07;                                                           % vacuum permittivity
E0 = 1;                                                                     % incoming field amplitude
step = 1;                                                                   % step in wavelength
global L; L = 2;                                                            % number of spheres
global L2; L2 = [2; 2];                                                     % number of layers each sphere [L1; L2; ...]
global N; N = 4;                                                            % number of Bessel mode
global sz; sz = 2*L*N*(N+2);                                                % number of equation
global cap; cap = N*(N+2);                                                  % number of mode for a certain N

% Material reference selection (lambda in nm, frequency in THz)
% - Metal, choose between Johnson Christy (J_C), Johnson, or Weaver. To
% choose it just uncomment the 'load' line and the line below it (lambda).
% And comment the other reference, like this example.
load J_C
lambda = 1240./J_C(:,1);
% load Johnson
% lambda = 1000*Johnson(:,1);
% load Weaver
% lambda = 1000*Weaver(:,1);
% - Dielectric, keep one of Si or SiO2 uncommented if you don't use it. If
% you are using Si or SiO2 uncomment it. If you choose constant parameter,
% don't change this part.
% load Si
% lambdadi = 1000*Si(:,1);
SiO2 = csvread('Malitson.csv');
lambdadi = 1000*SiO2(:,1);

lambdamin = 300; lambdamax = 800;                                           % range of lambda

dex = find(lambda>lambdamin & lambda<lambdamax);
dexdi = find(lambdadi>lambdamin & lambdadi<lambdamax);
Lambda = lambda(min(dex):max(dex));
Lambdadi = lambda(min(dexdi):max(dexdi));
lambda = lambdamin:step:lambdamax;
lambdadi = lambdamin:step:lambdamax;
[baris,kolom] = size(lambda);
f = c./(lambda*1E+3);
%% System Parameter
N0 = 1.33;                                                                  % background's refractive index
R = [40 40];                                                                % radius of each sphere (nm)
D = [5];                                                                   % distance between the spheres (nm)
D_ = 1E-05;                                                                 % use this if you want the position in [0 0 0]

% Position of each sphere in [x1 y1 z1; x2 y2 z2; ... ]
% You can't make the position of sphere in the center of the system [0 0 0]
% If you want to do it then use D_ instead [D_ 0 0].
global pos; pos = [-(R(1)+D(1)/2) 0 0; (R(2)+D(1)/2) 0 0];

cent = [0 0 0];                                                             % center of the system
k = 2*pi*N0./lambda;                                                        % wavenumber
x = [k*R(1); k*R(2)];                                                       % size parameter

% Choose material to be used (relative refractive index)
% If you use many material, just change the variable name to M1, M2, ...
% For Johnson & Christy
M1 = (J_C(:,6)+(1i*J_C(:,7)))/N0;                                           % for gold
% M = (J_C(:,4)+(1i*J_C(:,5)))/N0;                                            % for silver
% M = (J_C(:,2)+(1i*J_C(:,3)))/N0;                                            % for copper
% For Johnson
% M = (Johnson(:,2)+(1i*Johnson(:,3)))/N0;                                    % for gold
% M = (Johnson(:,4)+(1i*Johnson(:,5)))/N0;                                    % for silver
% For Weaver
% M = (Weaver(:,2)+(1i*Weaver(:,3)))/N0;                                      % for gold
% M = (Weaver(:,4)+(1i*Weaver(:,5)))/N0;                                      % for silver
% For Dielectric
% M = (Si(:,2)+(1i*Si(:,3)))/N0;                                              % for silicon
% M = SiO2(:,2)/N0;                                                           % for silica
Mc = 1.37*ones(1,kolom)/N0;                                                 % constant parameter
% Just add/remove m/M based on your configuration for data interpolation. 
% If you use dielectric (not constant), remember to change the
% lambda/Lambda to lambdadi/Lambdadi.
m1 = M1(min(dex):max(dex));
M1 = interp1(Lambda,m1,lambda,'pchip');
% m2 = M2(min(dex):max(dex));
% M2 = interp1(Lambdadi,m2,lambdadi,'pchip');

mu = 1;                                                                     % permeability
mur = ones(1,L);                                                            % relative permeability

% Multilayer Parameters
% 1. Radius of each layer of each sphere [R1 R2 ...]
Rj(1,:) = [38 40];
Rj(2,:) = [38 40];
% Rj(3,:) = ... just add more based on your configuration
% 2. Size parameter of each layer of each sphere [k*R k*R ...]
xj(1,:,:) = [k.*Rj(1,1); k.*Rj(1,2)];
xj(2,:,:) = [k.*Rj(2,1); k.*Rj(2,2)];
% xj(3,:,:) = ... just add more based on your configuration
M0 = ones(1,kolom);
M00 = zeros(1,kolom);
% 3. The relative refractive index compilation. You have to input M0 in the
% last index. The dimension of matrix of each ball (number of layers) have 
% to be same so use M00 to make it happen.
Mj(1,:,:) = [M1; Mc; M0];
Mj(2,:,:) = [M1; Mc; M0];
% Mj(3,:,:) = ... just add more based on your configuration
%% Functions to use
J = @(n,x) sqrt(pi./(2*x)).*besselj(n+0.5,x);                               % spherical Bessel
h = @(n,x) sqrt(pi/(2*x))*(besselj(n+0.5,x) + 1i*bessely(n+0.5,x));         % spherical Hankel
% Riccati-Bessel function
psi = @(n,x) x.*J(n,x);
dpsi = @(n,x) J(n,x) + (x/(2*n+1))*(n*J(n-1,x) - (n+1)*J(n+1,x));
xi = @(n,x) x.*h(n,x);
dxi = @(n,x) h(n,x) + (x/(2*n+1))*(n*h(n-1,x) - (n+1)*h(n+1,x));
% Multilayer functions
D1 = @(n,x) (dpsi(n,x)/psi(n,x));
D3 = @(n,x) (dxi(n,x)/xi(n,x));
Ra = @(n,x) (psi(n,x)/xi(n,x));
% Other functions
findja = @(i) ceil(i/cap);                                                  % what sphere?
findna = @(i,j) floor(sqrt(i-(j-1)*cap));                                   % what n?
findma = @(i,j,n) -n+(i-(j-1)*cap-n^(2));                                   % what m?
findjb = @(i) ceil((i-sz/2)/cap);                                           
findnb = @(i,j) floor(sqrt(i-sz/2-(j-1)*cap));
findmb = @(i,j,n) -n+(i-sz/2-(j-1)*cap-n^(2));
d = @(i,j) sqrt((pos(j,1)-pos(i,1))^2 + (pos(j,2)-pos(i,2))^2 +...
    (pos(j,3)-pos(i,3))^2);                                                 % distance between the spheres
the = @(i,j) acos((pos(j,3)-pos(i,3))/d(i,j));                              % polar angle between the spheres

%% Calculation (don't change this part)
% Step 1: Calculate the parameters between balls (distance, theta, phi)
for l1=1:L
    dc(l1) = sqrt((pos(l1,1)-cent(1))^2 + (pos(l1,2)-cent(2))^2 +...
        (pos(l1,3)-cent(3))^2);
    thetac(l1) = acos(-(pos(l1,3)+cent(3))/dc(l1));
    phic(l1) = atan2(-(pos(l1,2)+cent(2)),(-pos(l1,1)+cent(1)));
    for l2=1:L
        d0(l1,l2) = d(l1,l2); theta0(l1,l2) = the(l1,l2); phi0(l1,l2) = phi(l1,l2);
    end
end
for i=1:length(lambda)
    % Step 2: Calculate the scattering for each isolated ball and the
    % incoming wave coefficient
    progressbar(i/length(lambda))
    anj = zeros(N,L); bnj = zeros(N,L); pmnjj = zeros(N,2*N+1,L); qmnjj = zeros(N,2*N+1,L);
    for j=1:L
        for n=1:N
            AAA(n,1) = 0;
            BBB(n,1) = 0;
            for l=1:L2(j)
                Ha(n,l) = (Ra(n,Mj(j,l,i).*xj(j,l,i)).*D1(n,Mj(j,l,i).*xj(j,l,i))-AAA(n,l).*...
                    D3(n,Mj(j,l,i).*xj(j,l,i)))/(Ra(n,Mj(j,l,i).*xj(j,l,i))-AAA(n,l));
                Hb(n,l) = (Ra(n,Mj(j,l,i).*xj(j,l,i)).*D1(n,Mj(j,l,i).*xj(j,l,i))-BBB(n,l).*...
                    D3(n,Mj(j,l,i).*xj(j,l,i)))/(Ra(n,Mj(j,l,i).*xj(j,l,i))-BBB(n,l));
                AAA(n,l+1) = Ra(n,Mj(j,l+1,i).*xj(j,l,i)).*(Mj(j,l+1,i).*Ha(n,l)-...
                    Mj(j,1,i).*D1(n,Mj(j,l+1,i).*xj(j,l,i)))/(Mj(j,l+1,i).*Ha(n,l)-...
                    Mj(j,1,i).*D3(n,Mj(j,l+1,i).*xj(j,l,i)));
                BBB(n,l+1) = Ra(n,Mj(j,l+1,i).*xj(j,l,i)).*(Mj(j,1,i).*Hb(n,l)-...
                    Mj(j,l+1,i).*D1(n,Mj(j,l+1,i).*xj(j,l,i)))/(Mj(j,1,i).*Hb(n,l)-...
                    Mj(j,l+1,i).*D3(n,Mj(j,1,i).*xj(j,l,i)));
            end
            anj(n,j) = ((Ha(n,L2(j))/Mj(j,L2(j),i)+n/xj(j,L2(j),i))*psi(n,xj(j,L2(j),i))-psi(n-1,xj(j,L2(j),i)))/...
                ((Ha(n,L2(j))/Mj(j,L2(j),i)+n/xj(j,L2(j),i))*xi(n,xj(j,L2(j),i))-xi(n-1,xj(j,L2(j),i)));
            bnj(n,j) = ((Hb(n,L2(j))*Mj(j,L2(j),i)+n/xj(j,L2(j),i))*psi(n,xj(j,L2(j),i))-psi(n-1,xj(j,L2(j),i)))/...
                ((Hb(n,L2(j))*Mj(j,L2(j),i)+n/xj(j,L2(j),i))*xi(n,xj(j,L2(j),i))-xi(n-1,xj(j,L2(j),i)));
            for m=-n:1:n
                if m == -1
                    pmnjj(n,N+1+m,j) = -exp(1i*k(i)*pos(j,3))/(2*n*(n+1));
                    qmnjj(n,N+1+m,j) = -pmnjj(n,N+1+m,j);
                elseif m == 1
                    pmnjj(n,N+1+m,j) = exp(1i*k(i)*pos(j,3))/2;
                    qmnjj(n,N+1+m,j) = pmnjj(n,N+1+m,j);
                end
            end
        end
    end
    % Step 3: Calculate the translation coefficient of each ball
    AT = zeros(N,2*N+1,N,2*N+1,L,L);  BT = zeros(N,2*N+1,N,2*N+1,L,L);
    AT_ = zeros(N,2*N+1,N,2*N+1,L,L);  BT_ = zeros(N,2*N+1,N,2*N+1,L,L);
    for l1=1:L
        for l2=1:L
            if l1 ~= l2
                for n=1:N
                    for m=-n:1:n
                        for nu=1:N
                            for mu=-nu:1:nu
                                AT(n,N+m+1,nu,N+mu+1,l1,l2) = A(m,n,mu,nu,k(i)*d0(l1,l2),theta0(l1,l2),phi0(l1,l2));
                                BT(n,N+m+1,nu,N+mu+1,l1,l2) = B(m,n,mu,nu,k(i)*d0(l1,l2),theta0(l1,l2),phi0(l1,l2));
                                AT_(n,N+m+1,nu,N+mu+1,l1) = A_(m,n,mu,nu,k(i)*dc(l1),thetac(l1),phic(l1));
                                BT_(n,N+m+1,nu,N+mu+1,l1) = B_(m,n,mu,nu,k(i)*dc(l1),thetac(l1),phic(l1));
                            end
                        end
                    end
                end
            end
        end
    end
    % Step 4: Make scattering coefficient matrix to be solved
    AA = eye(sz);
    b = zeros(sz,1);
    for g=1:sz
        if g <= sz/2
            l = findja(g); 
            n = findna(g,l); 
            m = findma(g,l,n);
            b(g,1) = anj(n,l)*pmnjj(n,N+1+m,l);
            for j=1:sz
                if j <= sz/2
                    lu = findja(j); 
                    nu = findna(j,lu); 
                    mu = findma(j,lu,nu); 
                    if lu ~= l
                        AA(g,j) = anj(n,l)*AT(n,N+m+1,nu,N+mu+1,lu,l);
                    end
                else
                    lu = findjb(j); 
                    nu = findnb(j,lu); 
                    mu = findmb(j,lu,nu); 
                    if lu ~= l
                        AA(g,j) = anj(n,l)*BT(n,N+m+1,nu,N+mu+1,lu,l);
                    end
                end
            end
        else
            l = findjb(g); 
            n = findnb(g,l); 
            m = findmb(g,l,n);
            b(g,1) = bnj(n,l)*qmnjj(n,N+1+m,l);
            for j=1:sz
                if j <= sz/2
                    lu = findja(j); 
                    nu = findna(j,lu); 
                    mu = findma(j,lu,nu); 
                    if lu ~= l
                        AA(g,j) = bnj(n,l)*BT(n,N+m+1,nu,N+mu+1,lu,l);
                    end
                else
                    lu = findjb(j); 
                    nu = findnb(j,lu); 
                    mu = findmb(j,lu,nu); 
                    if lu ~= l
                        AA(g,j) = bnj(n,l)*AT(n,N+m+1,nu,N+mu+1,lu,l);
                    end
                end
            end
        end
    end
    % Step 5: Solve for amnj and bmnj
    W(:,1) = bicgstab(AA,b,1E-06,1000);
    % Step 6: Separate amnj and bmnj
    for g=1:sz
        if g <= sz/2
            j = findja(g);
            n = findna(g,j);
            m = findma(g,j,n);
            amnj(n,N+m+1,j,i) = W(g);
        else
            j = findjb(g);
            n = findnb(g,j);
            m = findmb(g,j,n);
            bmnj(n,N+m+1,j,i) = W(g);
        end
    end
    % Step 7: Calculate the total coefficient amn and bmn
    for n=1:N
        for m=-n:1:n
            sum1 = 0;
            sum2 = 0;
            for l=1:L
                for nu=1:N
                    for mu=-nu:1:nu
                        sum1 = sum1 + amnj(nu,N+1+mu,l,i)*AT_(n,N+m+1,nu,N+mu+1,l)...
                            + bmnj(nu,N+1+mu,l,i)*BT_(n,N+m+1,nu,N+mu+1,l);
                        sum2 = sum2 + amnj(nu,N+1+mu,l,i)*BT_(n,N+m+1,nu,N+mu+1,l)...
                            + bmnj(nu,N+1+mu,l,i)*AT_(n,N+m+1,nu,N+mu+1,l);
                    end
                end
            end
            amn(n,N+1+m,i) = sum1;
            bmn(n,N+1+m,i) = sum2;
            amn_(n,N+1+m,i) = -2*sqrt(pi)*(1i)^(-2*m)*sqrt((2*n+1)*factorial(n-m)/factorial(n+m))...
                *amn(n,N+1+m,i)/k(i)^2;
            bmn_(n,N+1+m,i) = -2*sqrt(pi)*(1i)^(-2*m)*sqrt((2*n+1)*factorial(n-m)/factorial(n+m))...
                *bmn(n,N+1+m,i)/k(i)^2;
        end
    end
	% Step 8: Calculate the cross section
    scn = zeros(2,N); se = 0; sc = 0;
    for n=1:N
        for m=-n:1:n
              sc = sc + n*(n+1)*(2*n+1)*((abs(amn(n,N+1+m,i)))^2 + (abs(bmn(n,N+1+m,i)))^2)...
                      *factorial(n-m)/factorial(n+m);
              se = se + n*(n+1)*(2*n+1)*real(conj(pmnjj(n,N+1+m,1))*amn(n,N+1+m,i)...
                   + conj(qmnjj(n,N+1+m,1))*bmn(n,N+1+m,i))*factorial(n-m)/factorial(n+m);
        end
    end
    Csca(i) = (4*pi/(k(i)^2))*sc;
    Cext(i) = (4*pi/(k(i)^2))*se;
    Cabs(i) = Cext(i) - Csca(i);
    % Step 9: Calculate the multipole expansion
    % Electric dipole
    px = amn_(1,N+1+1,i) - amn_(1,N-1+1,i); 
    py = 1i*(amn_(1,N+1+1,i) + amn_(1,N-1+1,i));
    pz = -sqrt(2)*amn_(1,N+0+1,i);
    % Magnetic dipole
    mx = bmn_(1,N+1+1,i) - bmn_(1,N-1+1,i); 
    my = 1i*(bmn_(1,N+1+1,i) + bmn_(1,N-1+1,i));
    mz = -sqrt(2)*bmn_(1,N+0+1,i);
    % Electric quadrupole
    Qxx = 1i*(amn_(2,N+2+1,i) + amn_(2,N-2+1,i)) - 0.5*1i*sqrt(6)*amn_(2,N+0+1,i);
    Qxy = amn_(2,N-2+1,i) - amn_(2,N+2+1,i);
    Qxz = 1i*(amn_(2,N-1+1,i) - amn_(2,N+1+1,i));
    Qyx = amn_(2,N-2+1,i) - amn_(2,N+2+1,i);
    Qyy = -1i*(amn_(2,N+2+1,i) + amn_(2,N-2+1,i)) - 0.5*1i*sqrt(6)*amn_(2,N+0+1,i);
    Qyz = amn_(2,N-1+1,i) + amn_(2,N+1+1,i);
    Qzx = 1i*(amn_(2,N-1+1,i) - amn_(2,N+1+1,i));
    Qzy = amn_(2,N-1+1,i) + amn_(2,N+1+1,i);
    Qzz = 1i*sqrt(6)*amn_(2,N+0+1,i);
    % Magnetic quadrupole
    Nxx = 1i*(bmn_(2,N+2+1,i) + bmn_(2,N-2+1,i)) - 0.5*1i*sqrt(6)*bmn_(2,N+0+1,i);
    Nxy = bmn_(2,N-2+1,i) - bmn_(2,N+2+1,i);
    Nxz = 1i*(bmn_(2,N-1+1,i) - bmn_(2,N+1+1,i));
    Nyx = bmn_(2,N-2+1,i) - bmn_(2,N+2+1,i);
    Nyy = -1i*(bmn_(2,N+2+1,i) + bmn_(2,N-2+1,i)) - 0.5*1i*sqrt(6)*bmn_(2,N+0+1,i);
    Nyz = bmn_(2,N-1+1,i) + bmn_(2,N+1+1,i);
    Nzx = 1i*(bmn_(2,N-1+1,i) - bmn_(2,N+1+1,i));
    Nzy = bmn_(2,N-1+1,i) + bmn_(2,N+1+1,i);
    Nzz = 1i*sqrt(6)*bmn_(2,N+0+1,i);
    % The cross-section
    Csde(i) = k(i)^2*(0.5*(abs(px+1i*py))^2 + (abs(pz))^2 + 0.5*(abs(px-1i*py))^2);
    Csdm(i) = k(i)^2*(0.5*(abs(mx+1i*my))^2 + (abs(mz))^2 + 0.5*(abs(mx-1i*my))^2);
    Csqe(i) = k(i)^2*(0.5*(abs(2*Qxx+2*1i*Qxy-Qzz))^2 + (abs(Qxz+1i*Qyz))^2 + (abs(1i*Qyz-Qxz))^2 ...
        +0.5*(abs(2*Qxx-2*1i*Qxy+Qzz))^2 + 2*(abs(Qzz))^2/3);
    Csqm(i) = k(i)^2*(0.5*(abs(2*Nxx+2*1i*Nxy-Nzz))^2 + (abs(Nxz+1i*Nyz))^2 + (abs(1i*Nyz-Nxz))^2 ...
        +0.5*(abs(2*Nxx-2*1i*Nxy+Nzz))^2 + 2*(abs(Nzz))^2/3);
 end
%% Output Plotting (Choose one or more)

% 1. Plot cross-section vs lambda
% figure
% plot(lambda,Csca,'Color','k','LineWidth',2)
% hold on
% plot(lambda,Cext,'LineWidth',2)
% hold on
% plot(lambda,Cabs,'Color','r','LineWidth',2)
% xlabel('\lambda (nm)'); ylabel('scattering cross section (nm^2)')
% legend ('C_{sca}','C_{ext}','C_{abs}')

% 2. Plot cross-section vs frequency
% figure
% plot(f,Csca,'Color','k','LineWidth',2)
% hold on
% plot(f,Cext,'Color','g','LineWidth',2)
% hold on
% plot(f,Cabs,'Color','r','LineWidth',2)
% xlabel('f (THz)'); ylabel('cross section (nm^2)')
% legend ('C_{sca}','C_{ext}','C_{abs}')

% 3. Plot multipole expansion
figure
plot(lambda,Csca,'LineWidth',2)
hold on
plot(lambda,Csde,'LineWidth',2)
hold on
plot(lambda,Csdm,'LineWidth',2)
hold on
plot(lambda,Csqe,'LineWidth',2)
hold on
plot(lambda,Csqm,'LineWidth',2)
xlabel('\lambda (nm)'); ylabel('cross section (nm^2)')
legend ('C_{sca}','C_{dip,e}','C_{dip,m}','C_{kuad,e}','C_{kuad,m}')

% 4. Polar Plot
% theta = 0:pi/360:2*pi; phi_ = PHI; r_ = 2000; th = 1;
% S1 = zeros(1,length(theta)); S2 = zeros(1,length(theta));
% for i=1:length(theta)
%     if theta(i) > pi & th == 1
%         th = th + 1;
%         phi_ = phi_ + pi;
%     end
%     sr = 0; st = 0; sp = 0;
%     for n=1:N
%         for m=-n:1:n
%             Mmnt = 1i*pimn(m,n,cos(theta(i)))*h(n,k(in)*r_)*exp(1i*m*phi_);
%             Mmnp = -taumn(m,n,cos(theta(i)))*h(n,k(in)*r_)*exp(1i*m*phi_);
%             Nmnr = n*(n+1)*P(m,n,cos(theta(i)))*h(n,k(in)*r_)*exp(1i*m*phi_)/(k(in)*r_);
%             Nmnt = taumn(m,n,cos(theta(i)))*dxi(n,k(in)*r_)*exp(1i*m*phi_)/(k(in)*r_);
%             Nmnp = 1i*pimn(m,n,cos(theta(i)))*dxi(n,k(in)*r_)*exp(1i*m*phi_)/(k(in)*r_);
%             Emn = (1i)^n*E0*(2*n+1)*factorial(n-m)/factorial(n+m);
%             Er(n,N+m+1) = 1i*Emn*amn(n,N+m+1,in)*Nmnr; 
%             Et(n,N+m+1) = 1i*Emn*(amn(n,N+m+1,in)*Nmnt + bmn(n,N+m+1,in)*Mmnt); 
%             Ep(n,N+m+1) = 1i*Emn*(amn(n,N+m+1,in)*Nmnp + bmn(n,N+m+1,in)*Mmnp);      
%             sr = sr +  Er(n,N+m+1); st = st + Et(n,N+m+1); sp = sp + Ep(n,N+m+1);
%         end
%     end
%     Etr(i) = sr; Ett(i) = st; Etp(i) = sp;
%     I(i) = (abs(Etr(i)))^2 + (abs(Ett(i)))^2 + (abs(Etp(i)))^2;
% end
% Ifwd = I(find(theta==0)); Ibck = I(find(theta==pi));
% ratio = Ifwd/Ibck
% figure
% hold on
% polarplot(theta,I,'LineWidth',1.5)
% ax = gca; % current axes
% ax.RTickLabel = [];
% hold on
toc
%% Other Functions
function y = phi(i,j) %mencari nilai phi
global pos;
d = sqrt((pos(j,1)-pos(i,1))^2 + (pos(j,2)-pos(i,2))^2 + (pos(j,3)-pos(i,3))^2); 
The = acos((pos(j,3)-pos(i,3))/d);
if The == 0
    The = The + 1E-10;
elseif The == pi
    The = The - 1E-10;
end
if (pos(j,2) - pos(i,2)) >= 0
    y = real(acos((pos(j,1)-pos(i,1))/(d*sin(The))));
else %harus bisa sampai 2 pi
    y = 2*pi - real(acos((pos(j,1)-pos(i,1))/(d*sin(The))));
end
end

function y = pimn(m,n,x)
if x == 1
    x = x - 1E-05;
elseif x == -1
    x = x + 1E-05;
end
y = m*P(m,n,x)/sqrt(1-x^2);
end

function y = taumn(m,n,x)
if x == 1
    x = x - 1E-05;
elseif x == -1
    x = x + 1E-05;
end
y = dP(m,n,x);
end

function y = dk(i,j)
if i == j
    y = 1;
else
    y = 0;
end
end