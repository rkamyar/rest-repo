

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stability Analysis on polytopes using piecewise polynomial Lyapunov
%functions synthesized using Handelman forms

%By Chaitanya Murti and Reza Kamyar
%9/18/2014
%For any questions or comments, please contact us at cmurti@hawk.iit.edu or 
%rkamyar@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc 
clear all
close all
tic
disp('Setting up the polytope and Lyapunov functions');

disp('Initial setup')

V.deg = 16;                                          %Degree of Lyapunov function
V1.deg = V.deg;                                     %Degree of Lyapunov function
V2.deg = V.deg;                                     %Degree of Lyapunov function
V3.deg = V.deg;                                     %Degree of Lyapunov function
V4.deg = V.deg;                                     %Degree of Lyapunov function
A.deg = 3;                                          %Degree of Vec field
Vdot.deg = V.deg+A.deg-1;                           %Degree of Vdot
V1dot.deg = V.deg+A.deg-1;                          %Degree of Vdot
V2dot.deg = V.deg+A.deg-1;                          %Degree of Vdot
V3dot.deg = V.deg+A.deg-1;                          %Degree of Vdot
V4dot.deg = V.deg+A.deg-1;                          %Degree of Vdot
n_polytope = 3;                                     %Number of sides in the polytope
n_states = 2;                                       %No. of state variables
V.monom_coeffs = Coeff_total(n_states,V.deg);       %Coeffs of V(x) in monom. basis
V.poly_coeffs = Coeff_total(n_polytope,V.deg);     %Coeffs of V(x) in poly. basis
Vdot.monom_coeffs = Coeff_total(n_states,Vdot.deg); %Coeffs of Vdot(x) in monom. basis
Vdot.poly_coeffs = Coeff_total(n_polytope,Vdot.deg); %Coeffs of Vdot(x) in poly. basis
V1.monom_coeffs = V.monom_coeffs; 
V2.monom_coeffs = V.monom_coeffs;
V3.monom_coeffs = V.monom_coeffs;
V4.monom_coeffs = V.monom_coeffs;
V1.poly_coeffs = V.poly_coeffs;
V2.poly_coeffs = V.poly_coeffs;
V3.poly_coeffs = V.poly_coeffs;
V4.poly_coeffs = V.poly_coeffs;
V1dot.monom_coeffs = Vdot.monom_coeffs ;
V2dot.monom_coeffs = Vdot.monom_coeffs ;
V3dot.monom_coeffs = Vdot.monom_coeffs ;
V4dot.monom_coeffs = Vdot.monom_coeffs ;
V1dot.poly_coeffs = Vdot.poly_coeffs ;
V2dot.poly_coeffs = Vdot.poly_coeffs ;
V3dot.poly_coeffs = Vdot.poly_coeffs ;
V4dot.poly_coeffs = Vdot.poly_coeffs ;

%epsilon = 0.05;
epsilon = 0.01;
%Here, we choose the vertices of the polytope in which we wish to perform
%stability analysis. Alternatively, one can simply provide vectors of
%coefficients of the polynomials defining the facets.



P = 1.64*[-1.31,0.18];         Q = 1.64*[0.56,1.92];
       
           O = [0,0];
           
S = 1.64*[-0.56,-1.92];         R = 1.64*[1.31,-0.18];

%Here, we choose the LP solver

%solver='linprog';
solver='SeDuMi';
% solver='GUROBI';

%Here,we set up all the facets of the D-decomposition of the polytope
%A vector of coefficients is returned.

PO = handelman_linemaker_2d(P,O);
OQ = handelman_linemaker_2d(O,Q);
PQ = handelman_linemaker_2d(P,Q);
OR = handelman_linemaker_2d(O,R);
QR = handelman_linemaker_2d(Q,R);
SR = handelman_linemaker_2d(S,R);
PS = handelman_linemaker_2d(P,S);
SO = handelman_linemaker_2d(S,O);

%Set up the D-Decomposition of the polytope. Here we use 4 triangles. Each
%polytope has it's facets stored in a cell.

%Set up first polytope
l1.v{1} = -PQ;
l1.v{2} = PO;
l1.v{3} = OQ;

%Set up second polytope
l2.v{1} = -QR;
l2.v{2} = OR;
l2.v{3} = -OQ;

%Set up third polytope
l3.v{1} = SR;
l3.v{2} = -OR;
l3.v{3} = -SO;

%Set up fourth polytope
l4.v{1} = PS;
l4.v{2} = -PO;
l4.v{3} = SO;

%Here, we construct a matrix of coefficients of the facets for each
%polytope.
for j=1:3
    temp = l1.v{j};
    LLL1(j,:) = temp;
end

for j=1:3
    temp = l2.v{j};
    LLL2(j,:) = temp;
end

for j=1:3
    temp = l3.v{j};
    LLL3(j,:) = temp;
end

for j=1:3
    temp = l4.v{j};
    LLL4(j,:) = temp;
end

%%
%Set up the vector field to be analyzed. We create a vector containing coefficients of a polynomial of the degree of the chosen vector field.

disp('Setting up f(x)')

A.exp = lex_exps(n_states,A.deg);
A.coeff_total = Coeff_total(n_states,A.deg);
A.coeffs = zeros(n_states,Coeff_total(n_states,A.deg));

j1 =2;
j2 = 1;
% 


k = lex_index_nh([0 1]);
A.coeffs(1,k) = -1;
A.coeffs(2,k) = -1;
j = lex_index_nh([0 1]);
k = lex_index_nh([1 0]);
A.coeffs(j,k) = 1;
k = lex_index_nh([2 1]);
A.coeffs(j,k) = 1;
% 


%%
%Testing which facets have nonzero coeffs. This allows us to construct the map R_i(b_i,d) later.

facet_test_vec = [1,0,0].';


facet_rule_1 = zeros(size(l1.v));
for i=1:length(l1.v)
    I = l1.v{i}*facet_test_vec;
    if I==0
        facet_rule_1(i) = 1;
    end
end

facet_rule_2 = zeros(size(l2.v));
for i=1:length(l2.v)
    I = l2.v{i}*facet_test_vec;
    if I==0
        facet_rule_2(i) = 1;
    end
end

facet_rule_3 = zeros(size(l3.v));
for i=1:length(l3.v)
    I = l3.v{i}*facet_test_vec;
    if I==0
        facet_rule_3(i) = 1;
    end
end

facet_rule_4 = zeros(size(l4.v));
for i=1:length(l4.v)
    I = l4.v{i}*facet_test_vec;
    if I==0
        facet_rule_4(i) = 1;
    end
end



%%
%Setting up the expansion of V1 
disp('Setting up the expansion of V1')
 
V1.monom_exps = lex_exps(n_states,V.deg);
V1.poly_exps = lex_exps(n_polytope, V.deg);
G1.temp_coeffs = zeros(V.monom_coeffs,V1.poly_coeffs);
for i=1:V.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V1.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l1.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G1.temp_coeffs(1:K,i) = temp_coeff;
end
% [G1.coeffs,G1.exps,G1.coeff_total] = coeff_elim(G1.temp_coeffs,V1.poly_exps,[1 1 1]);
%Now, we remove the coefficients of those terms that are nonzero at the origin
 [G1.coeffs,G1.exps,G1.coeff_total] = coeff_elim(G1.temp_coeffs,V1.poly_exps,facet_rule_1);

%%

%Setting up the expansion of V2 
disp('Setting up the expansion of V2')
 
V2.monom_exps = lex_exps(n_states,V.deg);
V2.poly_exps = lex_exps(n_polytope, V.deg);
G2.temp_coeffs = zeros(V.monom_coeffs,V.poly_coeffs);
for i=1:V2.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V2.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l2.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G2.temp_coeffs(1:K,i) = temp_coeff;
end
%[G2.coeffs,G2.exps,G2.coeff_total] = coeff_elim(G2.temp_coeffs,V2.poly_exps,[1 1 1]);

%Now, we remove the coefficients of those terms that are nonzero at the origin
 [G2.coeffs,G2.exps,G2.coeff_total] = coeff_elim(G2.temp_coeffs,V2.poly_exps,facet_rule_2);

%%
%Setting up the expansion of V3 

disp('Setting up the expansion of V3 ')
 
V3.monom_exps = lex_exps(n_states,V.deg);
V3.poly_exps = lex_exps(n_polytope, V.deg);
G3.temp_coeffs = zeros(V.monom_coeffs,V.poly_coeffs);
for i=1:V3.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V3.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l3.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G3.temp_coeffs(1:K,i) = temp_coeff;
end
%[G3.coeffs,G3.exps,G3.coeff_total] = coeff_elim(G3.temp_coeffs,V3.poly_exps,[1 1 1]);

%Now, we remove the coefficients of those terms that are nonzero at the origin
 [G3.coeffs,G3.exps,G3.coeff_total] = coeff_elim(G3.temp_coeffs,V3.poly_exps,facet_rule_3);

%%

%Setting up the expnasion of V4

disp('Setting up the expnasion of V4')
 
V4.monom_exps = lex_exps(n_states,V.deg);
V4.poly_exps = lex_exps(n_polytope, V.deg);
G4.temp_coeffs = zeros(V.monom_coeffs,V.poly_coeffs);
for i=1:V4.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V4.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l4.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G4.temp_coeffs(1:K,i) = temp_coeff;
end
%[G4.coeffs,G4.exps,G4.coeff_total] = coeff_elim(G4.temp_coeffs,V4.poly_exps,[1 1 1]);
%Now, we remove the coefficients of those terms that are nonzero at the origin
 [G4.coeffs,G4.exps,G4.coeff_total] = coeff_elim(G4.temp_coeffs,V4.poly_exps,facet_rule_4);


%%
%Setting up the expansion of V1dot
disp('Setting up the expansion of V1dot')

V1dot.monom_exps = lex_exps(n_states,Vdot.deg);
V1dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L1.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V1dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V1dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l1.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L1.coeffs(1:K,i) = -temp_coeff;
end
[L1.coeffs,L1.exps,L1.coeff_total] = coeff_elim(L1.coeffs,V1dot.poly_exps,facet_rule_1);

%%
%Setting up the expansion of V2dot

disp('Setting up the expansion of V2dot')

V2dot.monom_exps = lex_exps(n_states,Vdot.deg);
V2dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L2.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V2dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V2dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l2.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L2.coeffs(1:K,i) = -temp_coeff;
end

[L2.coeffs,L2.exps,L2.coeff_total] = coeff_elim(L2.coeffs,V2dot.poly_exps,facet_rule_2);
%%
%Setting up the expansion of V3dot

disp('Setting up the expansion of V3dot')

V3dot.monom_exps = lex_exps(n_states,Vdot.deg);
V3dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L3.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V3dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V3dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l3.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L3.coeffs(1:K,i) = -temp_coeff;
end
[L3.coeffs,L3.exps,L3.coeff_total] = coeff_elim(L3.coeffs,V3dot.poly_exps,facet_rule_3);

%%
%Setting up the expansion of V4dot

disp('Setting up the expansion of V4dot')

V4dot.monom_exps = lex_exps(n_states,Vdot.deg);
V4dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L4.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V4dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V4dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l4.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L4.coeffs(1:K,i) = -temp_coeff;
end
[L4.coeffs,L4.exps,L4.coeff_total] = coeff_elim(L4.coeffs,V4dot.poly_exps,facet_rule_4);
%%
% Setting up the expansion of grad(V1)^T*f
disp('Setting up the expansion of grad(V1)^T*f')                                                                                                                                                                                                                                                                                                                              
H1.coeffs = zeros(V1dot.monom_coeffs,G1.coeff_total);
GradV1.exps = lex_exps(n_states,V.deg-1);
GradV1.coeff_total = Coeff_total(n_states,V.deg-1);

exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV1.coeff_total,G1.coeff_total);
    for j=1:V1.monom_coeffs
        temp_exp_1 = V1.monom_exps(j,:);
        temp_coeff_vec = G1.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            k=lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    temp_coeffs;
    I = 4-lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
%       
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        temp_A = A_temp(j);
        for l=1:GradV1.coeff_total
            temp_exp_4 = temp_exp_3+GradV1.exps(l,:);
            temp_coeff_h = temp_A*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H1.coeffs(k,:) = H1.coeffs(k,:)+temp_coeff_h;
        end
    end
end

%%
%
%Setting up the expansion of grad(V2)^T*f
disp('Setting up the expansion of grad(V2)^T*f')                                                                                                                                                                                                                                                                                                                               
H2.coeffs = zeros(Vdot.monom_coeffs,G2.coeff_total);
GradV2.exps = lex_exps(n_states,V.deg-1);
GradV2.coeff_total = Coeff_total(n_states,V.deg-1);
GradV2.prodcoeffs = ones(2,GradV2.coeff_total);
exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV2.coeff_total,G2.coeff_total);
    for j=1:V2.monom_coeffs
        temp_exp_1 = V2.monom_exps(j,:);
        temp_coeff_vec = G2.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            
            k = lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    
    I = 4 -lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
    
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        for l=1:GradV2.coeff_total
            temp_exp_4 = temp_exp_3+GradV2.exps(l,:);
            temp_coeff_h = A_temp(j)*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H2.coeffs(k,:) = H2.coeffs(k,:)+temp_coeff_h;
        end
    end
end


%%
%Setting up the expansion of grad(V3)^T*f                                                                                                                                                                                                                                            %Setting up the expansion of grad(V3)^T*f
disp('Setting up the expansion of grad(V3)^T*f')
                                                                                  
H3.coeffs = zeros(V3dot.monom_coeffs,G3.coeff_total);
GradV3.exps = lex_exps(n_states,V.deg-1);
GradV3.coeff_total = Coeff_total(n_states,V.deg-1);
GradV3.prodcoeffs = ones(2,GradV3.coeff_total);
exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV3.coeff_total,G3.coeff_total);
    for j=1:V3.monom_coeffs
        temp_exp_1 = V3.monom_exps(j,:);
        temp_coeff_vec = G3.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            k = lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    
    I = 4 -lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
    
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        for l=1:GradV3.coeff_total
            temp_exp_4 = temp_exp_3+GradV3.exps(l,:);
            temp_coeff_h = A_temp(j)*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H3.coeffs(k,:) = H3.coeffs(k,:)+temp_coeff_h;
        end
    end
end

%%
% 
%Setting up the expansion of grad(V4)^T*f 
disp('Setting up the expansion of grad(V4)^T*f')                                                                                                                                                                                                                                                              
H4.coeffs = zeros(Vdot.monom_coeffs,G4.coeff_total);
GradV4.exps = lex_exps(n_states,V.deg-1);
GradV4.coeff_total = Coeff_total(n_states,V.deg-1);
GradV4.prodcoeffs = ones(2,GradV4.coeff_total);
exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV4.coeff_total,G4.coeff_total);
    for j=1:V4.monom_coeffs
        temp_exp_1 = V4.monom_exps(j,:);
        temp_coeff_vec = G4.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            k = lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    
    I = 4 -lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
    
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        for l=1:GradV4.coeff_total
            temp_exp_4 = temp_exp_3+GradV4.exps(l,:);
            temp_coeff_h = A_temp(j)*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H4.coeffs(k,:) = H4.coeffs(k,:)+temp_coeff_h;
        end
    end
end


%%
%Setting up blocks for interface between P1 and P2

G12.coeffs = coeff_elim_zeros(G1.coeffs,G1.exps,[0 0 1]);
G21.coeffs = -coeff_elim_zeros(G2.coeffs,G2.exps,[0 0 1]);

%%
%Setting up blocks for interface between P2 and P3

G23.coeffs = coeff_elim_zeros(G2.coeffs,G2.exps,[0 1 0]);
G32.coeffs = -coeff_elim_zeros(G3.coeffs,G3.exps,[0 1 0]);

%%
%Setting up blocks for interface between P3 and P4

G34.coeffs = coeff_elim_zeros(G3.coeffs,G3.exps,[0 0 1]);
G43.coeffs = -coeff_elim_zeros(G4.coeffs,G4.exps,[0 0 1]);

%%
%Setting up blocks for interface between P4 and P1

G41.coeffs = coeff_elim_zeros(G4.coeffs,G4.exps,[0 1 0]);
G14.coeffs = -coeff_elim_zeros(G1.coeffs,G1.exps,[0 1 0]);

%%
%  %Set up Matrices for LP
% % 
% % 
% 

 disp('Setting up and solving the LP')


%Setting up equality constraints:

%Setting up outputs of maps G_i(b_i,d)
B1 = blkdiag(H1.coeffs,H2.coeffs,H3.coeffs,H4.coeffs);
%Setting up outputs of maps of F_i(c_i,d+d_f-1)
B2 = blkdiag(-L1.coeffs,-L2.coeffs,-L3.coeffs,-L4.coeffs);

%Setting up outputs of maps J_i(b_i,d,k)
zeroblock = zeros(V.monom_coeffs,G1.coeff_total);
B41 = [G12.coeffs,G21.coeffs];
B42 = [G23.coeffs,G32.coeffs];
B43 = [G34.coeffs,G43.coeffs];
B44 = [G14.coeffs,zeros(V.monom_coeffs,2*G1.coeff_total),G41.coeffs];
B4 = [B41,zeros(V.monom_coeffs,2*G1.coeff_total);zeros(V.monom_coeffs,G1.coeff_total),B42,zeros(V.monom_coeffs,G1.coeff_total);zeros(V.monom_coeffs,2*G1.coeff_total),B43;B44];

disp('size B4')



A1 = [B1;B4];

A2 = [B2;zeros(4*V.monom_coeffs,4*L1.coeff_total)];


% A2 = [B2;zeros(V1.monom_coeffs,4*Vdot.poly_coeffs)];
% A2 = B2;
Aeq = [A1,A2];
% Aeq = constraint_elim(Aeq);

A_ineq = eye(4*G1.coeff_total+4*L1.coeff_total); 
size(A_ineq);
Aeq_size = size(Aeq);
cmat = zeros(Aeq_size(2),1);
bmat = zeros(Aeq_size(1),1);
b_ineq = zeros(4*G1.coeff_total+4*L1.coeff_total,1);
% b_ineq(3) = epsilon;
% b_ineq(4) = epsilon;
% b_ineq(5) = epsilon;


%Epsilon that worked for r<1 := (3,4,5)
%Epsilon that worked for r=1 := (3,5)

for j=1:4
    b_ineq(3+G1.coeff_total*(j-1)) = epsilon;
    b_ineq(4+G1.coeff_total*(j-1)) = 0;
    b_ineq(5+G1.coeff_total*(j-1)) = epsilon;
end


% % % 
if (strcmp(solver,'linprog'))
    [XX,fval,exitflag] = linprog(cmat,-A_ineq,-b_ineq,Aeq,bmat);
    exitflag
    max(XX)
elseif (strcmp(solver,'SeDuMi'))
    [XX,YY,info] = sedumi(Aeq,bmat,cmat);
    max(XX)
    min(XX)
    info
elseif (strcmp(solver,'CPLEX'))
    [XX, fval, exitflag, output, lambda] = cplexlp (zeros(1,length(cmat)), -A_ineq, -b_ineq, Aeq, bmat);
    max(XX)
    min(XX)
    exitflag
elseif (strcmp(solver,'GUROBI'))
    solver_g
    XX=result1.x;
end
toc

%Obtaining the vector of coefficients of the Lyapunouv function
  temp = 0;
  XX1 = G1.coeffs*XX(temp+1:temp+G1.coeff_total);
  temp = temp+ G1.coeff_total;
  XX2 = G2.coeffs*XX(temp+1:G2.coeff_total+temp);
  temp = temp+G2.coeff_total;
  XX3 = G3.coeffs*XX(temp+1:G3.coeff_total+temp);
  temp = temp+G3.coeff_total;
  XX4 = G4.coeffs*XX(temp+1:G4.coeff_total+temp);

  Vexps = lex_exps(2,V.deg);


%Obtaining the vector of coefficients of the Lie derivative of the Lyapunouv function  
  temp = 0;
  DX1 = H1.coeffs*XX(temp+1:temp+G1.coeff_total);
  temp = temp+ G1.coeff_total;
  DX2 = H2.coeffs*XX(temp+1:G2.coeff_total+temp);
  temp = temp+G2.coeff_total;
  DX3 = H3.coeffs*XX(temp+1:G3.coeff_total+temp);
  temp = temp+G3.coeff_total;
  DX4 = H4.coeffs*XX(temp+1:G4.coeff_total+temp);
  temp = temp+G4.coeff_total;
  VDexps = lex_exps(2,Vdot.deg);

  
  


  LX1 = L1.coeffs*XX(temp+1:temp+L1.coeff_total);
    LL1 = XX(temp+1:temp+L1.coeff_total);
  temp = temp+ L1.coeff_total;
  LX2 = L2.coeffs*XX(temp+1:L2.coeff_total+temp);
  LL2 = XX(temp+1:L2.coeff_total+temp);
  temp = temp+L2.coeff_total;
  LX3 = L3.coeffs*XX(temp+1:L3.coeff_total+temp);
  LL3 = XX(temp+1:L3.coeff_total+temp);
  temp = temp+L3.coeff_total;
  LX4 = L4.coeffs*XX(temp+1:L4.coeff_total+temp);
  LL4 =XX(temp+1:L4.coeff_total+temp);


  
  temp = 0;
 
  GX12 = G12.coeffs*XX(temp+1:temp+G1.coeff_total);
  GX14 = G14.coeffs*XX(temp+1:temp+G1.coeff_total);
  temp = temp+ G1.coeff_total;
  GX21 = G21.coeffs*XX(temp+1:G2.coeff_total+temp);
  GX23 = G23.coeffs*XX(temp+1:G2.coeff_total+temp);
  temp = temp+G2.coeff_total;
  GX32 = G32.coeffs*XX(temp+1:G3.coeff_total+temp);
  GX34 = G34.coeffs*XX(temp+1:G3.coeff_total+temp);
  temp = temp+G3.coeff_total;
  GX43 = G43.coeffs*XX(temp+1:G4.coeff_total+temp);
  GX41 = G41.coeffs*XX(temp+1:G4.coeff_total+temp);
  temp = temp+G4.coeff_total;

  %Checking the coefficients: are they nonnegative? are they nonzero? if
  %each test gives a zero vector, the solution is acceptable.
  G12test = GX12+GX21;
  G23test = GX23+GX32;
  G34test = GX34+GX43;
  G41test = GX14+GX41;
