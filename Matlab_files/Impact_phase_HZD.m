clear 
clc

%Use one absolute and four relative coordinates and their derivatives
syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 real

%Masses and inertia about com of the links
syms MT Mf Mt IT If It real

%Length parameters
syms lT lf lt PT Pf Pt real

MT=20;
Mf=6.8;
Mt=3.2;
lT=0.625;
lf=0.4;
lt=0.4;
IT=2.22;
If=1.08;
It=0.93;
PT=0.2;
Pf=0.163;
Pt=0.128;

%general parameters
syms g real

g=9.8;

qs = [q1;q2;q3;q4;q5];
qsd= [dq1;dq2;dq3;dq4;dq5];

qe = [q1;q2;q3;q4;q5;q6;q7];
qed= [dq1;dq2;dq3;dq4;dq5;dq6;dq7];

%pcm of each link
%t_stance
pt_x1= q6-(lt-Pt)*cos(q2+q4+q1-(pi/2));
%pt_x1= -Pt*sin(q2+q4+q1);
pt_y1= q7+(lt-Pt)*sin(q2+q4+q1-(pi/2));
%pt_y1=  -Pt*cos(q2+q4+q1);

%f_stance
pf_x1= q6-lt*cos(q2+q4+q1-(pi/2))-(lf-Pf)*cos(q2+q1-(pi/2));
pf_y1= q7+lt*sin(q2+q4+q1-(pi/2))+(lf-Pf)*sin(q2+q1-(pi/2));

%Torso
pT_x =q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+PT*sin(q1);
pT_y =q7+lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))+PT*cos(q1);

%f_swing
pf_x2=q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+Pf*cos(q3+q1-(pi/2));
pf_y2=q7+lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))-Pf*sin(q3+q1-(pi/2));

%t_swing
pt_x2=q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+lf*cos(q3+q1-(pi/2))+Pt*sin(q5+q3+q1);
pt_y2=q7+lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))-lf*sin(q3+q1-(pi/2))+Pt*cos(q5+q3+q1);



%velocity calculaiton

pt_x1d=jacobian(pt_x1,qe)*qed;
pt_y1d=jacobian(pt_y1,qe)*qed;

pf_x1d=jacobian(pf_x1,qe)*qed;
pf_y1d=jacobian(pf_y1,qe)*qed;

pT_xd=jacobian(pT_x,qe)*qed;
pT_yd=jacobian(pT_y,qe)*qed;

pt_x2d=jacobian(pt_x2,qe)*qed;
pt_y2d=jacobian(pt_y2,qe)*qed;

pf_x2d=jacobian(pf_x2,qe)*qed;
pf_y2d=jacobian(pf_y2,qe)*qed;

%theta abs
th_T=q1;
th_f1=q1+q2;
th_t1=q1+q2+q4;
th_f2=q1+q3;
th_t2=q1+q3+q5;

%omega abs
th_Td=jacobian(th_T,qe)*qed;
th_f1d=jacobian(th_f1,qe)*qed;
th_f2d=jacobian(th_f2,qe)*qed;
th_t1d=jacobian(th_t1,qe)*qed;
th_t2d=jacobian(th_t2,qe)*qed;

%potential energy calculation

V= (MT*pT_y + Mf*(pf_y1+pf_y2) + Mt*(pt_y1+pt_y2))*g;

%kinetic energy calculation 

ke_t1= (Mt/2)*((pt_x1d)'*(pt_x1d)+ (pt_y1d)'*(pt_y1d))+(It/2)*(th_t1d'*th_t1d);
ke_t2= (Mt/2)*((pt_x2d)'*(pt_x2d)+ (pt_y2d)'*(pt_y2d))+(It/2)*(th_t2d'*th_t2d);
ke_f1= (Mf/2)*((pf_x1d)'*(pf_x1d)+ (pf_y1d)'*(pf_y1d))+(If/2)*(th_f1d'*th_f1d);
ke_f2= (Mf/2)*((pf_x2d)'*(pf_x2d)+ (pf_y2d)'*(pf_y2d))+(If/2)*(th_f2d'*th_f2d);
ke_T= (MT/2)*((pT_xd)'*(pT_xd)+ (pT_yd)'*(pT_yd))+(IT/2)*(th_Td'*th_Td);

K= ke_t1+ke_t2+ke_T+ke_f1+ke_f2;

De = jacobian(jacobian(K,qed),qed);
De=simplify(De);

%gamma 2 - end points of the impact leg just before impact

p2_x=q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+lf*cos(q3+q1-(pi/2))-lt*sin(q5+q3+q1);
p2_y=q7+lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))-lf*sin(q3+q1-(pi/2))-lt*cos(q5+q3+q1);

p2=[p2_x;p2_y];

E2=jacobian(p2,qe);


list_q_e  = {'q1','q(1)'; 'q2','q(2)'; 'q3','q(3)'; 'q4','q(4)';'q5','q(5)'; 'q6','q(6)';'q7','q(7)';'dq1','q(8)';'dq2','q(9)';'dq3','q(10)';'dq4','q(11)';'dq5','q(12)';'dq6','q(13)';'dq7','q(14)'};
write_func('D_q_e.m',{'q'},[list_q_e],{De,'De'});
write_func('E2_q_e.m',{'q'},[list_q_e],{E2,'E2'});




