% This Matlab code computes the ground state energy for helium (He) atom by solving the Roothaan equation using the
% cc-pVDZ basis set (restricted Hartree-Fock (rhf) calculation). 
%
% The core Hamiltonian matrix (H_core), overlap matrix (S_ov) and two-electron integrals (tei) (He_cc_pvdz_tei.txt) are computed 
% by my own developing code. An obtained total energy in atomic unit (au) is compared with
% that of the Psi4 computational chemistry package. 
% 
% Ref: A. Szabo and N. S. Ostlund "Modern Quantum Chemistry" book.  
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% April 11, 2024 & University of North Dakota 
%
function [] = He_rhf_cc_pvdz
clear; clc
%
format long
%
S_ov = [1.00000061 0.88884518 0.         0.         0.        ;
        0.88884518 1.         0.         0.         0.        ;
        0.         0.         1.         0.         0.        ;
        0.         0.         0.         1.         0.        ;
        0.         0.         0.         0.         1.        ];

%
H_core = [-1.94101369 -1.58026098  0.          0.          0.        ;
          -1.58026098 -1.29467114  0.          0.          0.        ;
           0.          0.          0.78499729  0.          0.        ;
           0.          0.          0.          0.78499729  0.        ;
           0.          0.          0.          0.          0.78499729];
%
%
dim = 5; % size of basis sets & (4s,1p) -> [2s,1p] = 2x1 + 1x3 = 5
%
N_el = 2.;               % number of electron 
%
itermax = 100; tol = 1e-8;
%
tei_n = 625;             % = 5^4, .i.e., all values of TEI
%
read_tei_data = fopen('He_cc_pvdz_tei.txt', 'r');               % data of two-electron integral in atomic basis set
tei_data_n5 = textscan(read_tei_data, '%d %d %d %d %f');
%
p = zeros(tei_n,1); q = zeros(tei_n,1); r = zeros(tei_n,1); s = zeros(tei_n,1); vals = zeros(tei_n,1);
p(1:tei_n) = tei_data_n5{1};
q(1:tei_n) = tei_data_n5{2};
r(1:tei_n) = tei_data_n5{3};
s(1:tei_n) = tei_data_n5{4};
vals(1:tei_n) = tei_data_n5{5};
for i = 1:tei_n
    tei(p(i),q(i),r(i),s(i)) = vals(i);
%    tei(q(i),p(i),r(i),s(i)) = vals(i);    
%    tei(p(i),q(i),s(i),r(i)) = vals(i);    
%    tei(q(i),p(i),s(i),r(i)) = vals(i);   
    %
%    tei(r(i),s(i),p(i),q(i)) = vals(i);    
%    tei(s(i),r(i),p(i),q(i)) = vals(i);        
%    tei(r(i),s(i),q(i),p(i)) = vals(i);        
%    tei(s(i),r(i),q(i),p(i)) = vals(i);            
end
%
P_old = 0.5 * ones(dim,dim); % initial charge population
%
for iter = 1:itermax
    iter;
    P = P_old;
    %
    F = H_core;
    for p = 1:dim
        for q = 1:dim
            for r = 1:dim
                for s = 1:dim
                    F(p,q) = F(p,q) + P(r,s) * (tei(p,q,r,s) - 0.5.*tei(p,r,q,s));
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
    for i = 1:N_el/2
        for pp = 1:dim
            for qq = 1:dim
                P_new(pp,qq) = P_new(pp,qq) + 2*Vec(pp,i)*Vec(qq,i);
            end
        end
    end
    %
     if (abs(sum(sum(P_new-P_old))) < tol)
            break 
     end
    %        
    P_old = P_new;

end
%%%
[En(1), En(2)]; % [-0.914147924737976   1.397441760376133]
 En_0 = (sum(0.5*diag(P(:,:)*(H_core(:,:) + F(:,:))))) % -2.855160473239254au vs -2.855160477242742au = Psi4



%%%
return
end
