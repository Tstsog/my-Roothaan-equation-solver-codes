% This Matlab code computes the ground state energy for helium atom by solving the Roothaan equation using the
% Gaussian-type orbital (GTO) with n=2 (2 s-function)
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% July 12, 2023 & University of North Dakota 
%
function [] = he_scf_gto_2s
%
clear; clc; format long
itermax = 60; tol = 1e-12;
%
z_h = 2.; % nuclear charge for helium atom 
%
xi1 = 0.532149; d1 = 0.82559; % from S. Huzinaga, J. Chem. Phys. 42, 1293–1302 (1965), Table IX
xi2 = 4.097728; d2 = 0.28317;
%
N = 512; % number of grid points along radial distance r; you may change it
a = 0.0;  % starting point od coordinate r
b = 50.; % end point of cooridnate r; you may change it
%
[r,wr,D]=legDC2(N,a,b); % D is the differentation matrix of the first order
D1 = (2/(b-a))*D; rr = r;
%
% GTO: basis functions & chi = d * exp(-xi*r^2) * (2*xi/pi)^(3/4)
chi_1 = d1.*exp(-xi1.*r.^2).*(2.*xi1./pi).^(3/4);    
chi_2 = d2.*exp(-xi2.*r.^2).*(2.*xi2./pi).^(3/4);     

%%%%%%%%%%%%%%%%%%%%
% Hamiltonian matrix elements in atomic orbital basis functions 
[h11,s11] = H0_elements_ss(z_h,D1,r,wr,chi_1,chi_1);
[h12,s12] = H0_elements_ss(z_h,D1,r,wr,chi_1,chi_2);
[h22,s22] = H0_elements_ss(z_h,D1,r,wr,chi_2,chi_2);
%
%%%
dim = 2;
d_coef = zeros(1,dim);
d_coef(1) = d1*(2.*xi1./pi).^(3/4);
d_coef(2) = d2*(2.*xi2./pi).^(3/4);
%
xi_coef = zeros(1,dim);
xi_coef(1) = xi1;
xi_coef(2) = xi2;
%
H_core = [h11, h12; % the core hamiltonian: matrix elements
          h12, h22];
%
S_ov = [s11, s12; % overlap matrix elements  
        s12, s22];
%
P_old = 0.5 * ones(dim,dim); % initial charge population
%
for iter = 1:itermax
    iter
    P = P_old;
    %
    F = H_core;
    for p = 1:dim
        for q = 1:dim
            for r = 1:dim
                for s = 1:dim
                    F(p,q) = F(p,q) + P(r,s) * (tei_ssss(p,q,r,s, d_coef, xi_coef) - 0.5.*tei_ssss(p,r,q,s, d_coef, xi_coef));
                end
    
            end
    
        end
    end
    Ham_fock = F ;     % Fock matrix
    S_mat_fock = S_ov;

    [Vec,En] = eig(Ham_fock,S_mat_fock);                                     % Eigenvalue problem: F*c = En*S*c - Roothaan equation
    En = diag(En);
    [foo, ij] = sort(En);
    En = En(ij);
    [En(1), En(2)];  % orbital energies
    %
    Vec = Vec(:,ij);                       % expansion coefficients 
    %
    for i = 1:dim
        norm = 0.;
        for p = 1:dim
            for q = 1:dim
                norm = norm + Vec(p,i) * Vec(q,i) * S_ov(p,q);
            end
        end
        Vec(:,i) = Vec(:,i)/sqrt(norm);
    end
    %
    P_new = zeros(dim,dim);
    for i = 1:z_h/2
        for pp = 1:dim
            for qq = 1:dim
                P_new(pp,qq) = P_new(pp,qq) + 2*Vec(pp,i)*Vec(qq,i);
            end
        end
    end
    %
     if (abs(P_new-P_old) < tol)
            break 
     end
    %        
    P_old = P_new;

end
%%%

En_0 = (sum(0.5*diag(P(:,:)*(H_core(:,:) + F(:,:))))); % ground state energy in atomic unit
[En(1), En_0]
% [En(1), En_0] = -0.858910329320571  -2.747066128454680; vs [-0.858911, -2.7470661] from S. Huzinaga, J. Chem. Phys. 42, 1293–1302 (1965)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V1 = Vec(:,1); V2 = Vec(:,2);   % the expansion coefficients in basis functions
%%%
phi_1 = V1(1).*chi_1 + V1(2).*chi_2; % phi_1: molecular orbital basis function 
phi_2 = V2(1).*chi_1 + V2(2).*chi_2;
%%%
[h11,s11] = H0_elements_ss(z_h,D1,rr,wr,phi_1,phi_1);  % the core hamiltonian matrix elements in molecular orbital basis functions
[h12,s12] = H0_elements_ss(z_h,D1,rr,wr,phi_1,phi_2);
[h22,s22] = H0_elements_ss(z_h,D1,rr,wr,phi_2,phi_2);

%%%
P = Vec'; % charge density matrix
%
tei = zeros(dim,dim,dim,dim);  % two-electron integral (tei) in molecular orbital basis sets
for ii = 1:dim
    for jj = 1:dim
        for kk = 1:dim
            for ll = 1:dim
                for mm = 1:dim
                    for nn = 1:dim
                        for oo = 1:dim
                            for pp = 1:dim
                                tei(ii,jj,kk,ll) =  tei(ii,jj,kk,ll) + P(ii,mm)*P(jj,nn)*P(kk,oo)*P(ll,pp)*tei_ssss(mm,nn,oo,pp, d_coef, xi_coef);
                            end
                        end
                    end
                end
            end
        end
    end
end
%
E0_n = 2*h11 + tei(1,1,1,1); % the ground state energy E0_n = -2.747066128454680, from molecular orbital basis functions


%%%
return
end

%%%%%%%%%%%%%%%%
function [h11,s11] = H0_elements_ss(z_h,D1,r,wr,chi_i,chi_j)
% compute kinetic and potential energies and overlap matrix
%
T_11 = sum(wr.*(D1*chi_i).*(D1*chi_j).*r.*r) *(4*pi) ;
V_11 = sum(wr.*chi_i.*(-z_h).*chi_j.*r) * (4*pi) ;
h11 = 0.5*T_11 + V_11;
s11 = sum(wr.*chi_i.*chi_j.*r.*r) * (4*pi);
%
%%%
return
end
%%%

%%%
function [Q_pqrs] = tei_ssss(p,q,r,s, d_coef, xi_coef)
% analytical expression for the two-electron integral 
%
Q_pqrs_numer = 2.*pi.^(5/2);
Q_pqrs_denun = (xi_coef(p) + xi_coef(q))*(xi_coef(r) + xi_coef(s))*sqrt(xi_coef(p) + xi_coef(q) + xi_coef(r) + xi_coef(s));
Q_pqrs = Q_pqrs_numer/Q_pqrs_denun; 
%
Q_pqrs = d_coef(p) * d_coef(q) * d_coef(r) * d_coef(s) * Q_pqrs;
%
return
end




    function [xi,w,D]=legDC2(N,a,b)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % legDc.m
            %
            % Computes the Legendre differentiation matrix with collocation at the
            % Legendre-Gauss-Lobatto nodes.
            %
            % Reference:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            % Written by Greg von Winckel - 05/26/2004
            % Contact: gregvw@chtm.unm.edu
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % Truncation + 1
            N1=N+1;
            
            % CGL nodes
            xc=cos(pi*(0:N)/N)';
            
            % Uniform nodes
            xu=linspace(-1,1,N1)';
            
            % Make a close first guess to reduce iterations
            if N<3
                x=xc;
            else
                x=xc+sin(pi*xu)./(4*N);
            end
            
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1);
            
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            
            xold=2;
            while max(abs(x-xold))>eps
                
                xold=x;
                
                P(:,1)=1;    P(:,2)=x;
                
                for k=2:N
                    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                
                x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
            
            X=repmat(x,1,N1);
            Xdiff=X-X'+eye(N1);
            
            L=repmat(P(:,N1),1,N1);
            L(1:(N1+1):N1*N1)=1;
            D=(L./(Xdiff.*L'));
            D(1:(N1+1):N1*N1)=0;
            D(1)=(N1*N)/4;
            D(N1*N1)=-(N1*N)/4;
            
            % Linear map from[-1,1] to [a,b]
            xi=(a*(1-x)+b*(1+x))/2;        % added by Tsogbayar Tsednee
            
            % Compute the weights
            w=(b-a)./(N*N1*P(:,N1).^2);    % added by Tsogbayar Tsednee
            
    end

