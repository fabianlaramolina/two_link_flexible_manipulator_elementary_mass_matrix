%Paper Title: Stochastic Finite-Element Modeling of Flexible Manipulators 
%Dynamic Model of the One-Link Flexible Manipulator by using Finite Element
%Method
%Fabian Andres Lara-Molina
%Federal University of Technology - Paraná
%Lara.f8@gmail.com
%APPENDIX B
%2021

function Elementary_mass_matrix_Two_Link_Flexible_Manipulator_App_B()
clear all
clc
syms rho A E I %Link Proprieties - equal for link1 and link 2
syms l  % Length element /Link 1 element length = Link 2 element length
syms l1  % Length link 1
syms u3_l1_end  %Last translational elastic dof of Length link 1
syms u4_l1_end  %Last rotational elastic dof of Length link 1
syms theta_1 theta_2 dtheta_1 dtheta_2 %joint angles and time derivartives
syms d % x coordinate of r2_j=[d y_2j]^T=[(j-1)l+x  y_2j]^T
syms x2 %x coordinate second element 0<=x2<=l
syms j_  %Number of elements of link 2
syms Phi Phi_T y2 % Phi - shape functions   / y2=Phi*[psi_{ij}]'=Phi*[u1 u2 u3 u4]'
syms u1 u2 u3 u4 %Elastic Degrees of Freedom of the element
syms n  %Number of elements for link
syms mp_2 % Tip mass
% ==========================================================================
%Second Link Massa Matrix-----
[M_2 M_p_2]=Mass_Matrix_second_link(theta_1,theta_2,u3_l1_end,u4_l1_end,x2,l,j_,rho,A,l1,y2,n,mp_2);

psi_2j=[u1 u2 u3 u4];
Me=Sub_non_linear_term_link_2(x2,psi_2j,M_2,y2,l); %Function to include non-linear terms
                                                   %within the mass elementary matrix

disp('Mass Elementary Matrix of the Second Link:')
Me


%Substituindo end tip em M_p_2
M_p_2=subs(M_p_2,y2,psi_2j(end-1));
disp('Tip Mass Elementary Matrix - Second Link:')
M_p_2




function [M_i M_p]=Mass_Matrix_second_link(theta_1,theta_2,u3_l1_end,u4_l1_end,x2,l,j_,rho,A,l1,y2,n,mp_2)
%Second Element
%Mass Matrix-----

d=(j_-1)*l+x2;

T1=[cos(theta_1) -sin(theta_1);sin(theta_1) cos(theta_1)]; %R(theta_1_Z1)


r21=T1*[l1;u3_l1_end];  %End Cartesian coodinate of fist link w.r.t. frame O
dr21=fulldiff(r21,{theta_1,u3_l1_end});

Ang_t_2=theta_2+u4_l1_end;
T2=[cos(Ang_t_2) -sin(Ang_t_2);sin(Ang_t_2) cos(Ang_t_2)]; %R(theta_1+theta_2+u4_l1_end,Z3) - Rotational matrix
r22=T1*T2*[d;y2];
dr22=fulldiff(r22,{theta_1,theta_2,u3_l1_end,u4_l1_end,y2});

dr2=simplify(dr21+dr22);

drm2 = [-(y2*cos(theta_1 + theta_2 + u4_l1_end) + sin(theta_1 + theta_2 + u4_l1_end)*(x2 + l*(j_ - 1))+u3_l1_end*cos(theta_1) + l1*sin(theta_1)),...
 - sin(theta_1),... 
 - (y2*cos(theta_1 + theta_2 + u4_l1_end) + sin(theta_1 + theta_2 + u4_l1_end)*(x2 + l*(j_ - 1))),...
 -(y2*cos(theta_1 + theta_2 + u4_l1_end) + sin(theta_1 + theta_2 + u4_l1_end)*(x2 + l*(j_ - 1))),...
 -sin(theta_1 + theta_2 + u4_l1_end);...
(l1*cos(theta_1) - u3_l1_end*sin(theta_1)-y2*sin(theta_1 + theta_2 + u4_l1_end) + cos(theta_1 + theta_2 + u4_l1_end)*(x2 + l*(j_ - 1))),... 
+ cos(theta_1),...
- (y2*sin(theta_1 + theta_2 + u4_l1_end) - cos(theta_1 + theta_2 + u4_l1_end)*(x2 + l*(j_ - 1))),...
- (y2*sin(theta_1 + theta_2 + u4_l1_end) - cos(theta_1 + theta_2 + u4_l1_end)*(x2 + l*(j_ - 1))),... 
+cos(theta_1 + theta_2 + u4_l1_end)];
drm_Txdrm=simplify(drm2.'*drm2);   
%Substituindo dy=Phi*du e dy'=du'*Phi'
Phi=[1-3*x2^2/l^2+2*x2^3/l^3, x2-2*x2^2/l+x2^3/l^2, 3*x2^2/l^2-2*x2^3/l^3, -x2^2/l+x2^3/l^2];
Phi_T=Phi.';


drm2_t= simplify([eye(4,5);zeros(4,4) Phi_T]*drm_Txdrm*[eye(5,4) [zeros(4,4);Phi]]);

M_pp=subs(drm_Txdrm,[j_,x2],[n,l]);
M_p=sym(zeros(8,8)); % size of elementary matrix
M_p(1:4,1:4)=M_pp(1:4,1:4);
M_p(7,1:4)=M_pp(5,1:4);
M_p(1:4,7)=M_pp(1:4,5);
M_p(7,7)=M_pp(5,5);
M_p=M_p*mp_2;

M_i=int(rho*A*drm2_t,x2,[0 l]);

function M_i=Sub_non_linear_term_link_2(x2,Q,M_1,y2,l)
Phi=[1-3*x2^2/l^2+2*x2^3/l^3, x2-2*x2^2/l+x2^3/l^2, 3*x2^2/l^2-2*x2^3/l^3, -x2^2/l+x2^3/l^2];
m_1j=int(Phi,[0 l]);
M_1j=int(Phi.'*Phi,[0 l]);
% M_1j=l/420*[156 22*l 54 -13*l;22*l 4*l^2 13*l -3*l^2;54 13*l 156 -22*l;-13*l -3*l^2 -22*l 4*l^2]
yp2=simplify(Q*M_1j*Q.');
yp1=simplify(m_1j*Q.');
M_i=subs(M_1,[y2^2,y2],[yp2,yp1]);




