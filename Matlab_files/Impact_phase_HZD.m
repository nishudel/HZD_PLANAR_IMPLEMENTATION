clear 
clc

%Use one absolute and four relative coordinates and their derivatives
syms q1 q2 q3 q4 q5 q6 q7 q1d q2d q3d q4d q5d q6d q7d real

%Masses and inertia about com of the links
syms MT Mf Mt IT If It real

%Length parameters
syms lT lf lt PT Pf Pt real

%general parameters
syms g real

qs = [q1;q2;q3;q4;q5];
qsd= [q1d;q2d;q3d;q4d;q5d];
q  = [q1;q2;q3;q4;q5;q6;q7];
qd = [q1d;q2d;q3d;q4d;q5d;q6d;q7d];


%pcm of each link
%t_stance
pt_x1= q6-Pt*cos(q2+q4+q1-(pi/2));
%pt_x1= -Pt*sin(q2+q4+q1);
pt_y1= q7+ Pt*sin(q2+q4+q1-(pi/2));
%pt_y1=  -Pt*cos(q2+q4+q1);

%f_stance
pf_x1= q6-lt*cos(q2+q4+q1-(pi/2))-Pf*cos(q2+q1-(pi/2));
pf_y1= q7+ lt*sin(q2+q4+q1-(pi/2))+Pf*sin(q2+q1-(pi/2));

%Torso
pT_x =q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+PT*cos(q1);
pT_y =q7+ lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))+PT*sin(q1);

%f_swing
pf_x2=q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+Pf*cos(q3+q1-(pi/2));
pf_y2=q7+lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))-Pf*sin(q3+q1-(pi/2));

%t_swing
pt_x2=q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+lf*cos(q3+q1-(pi/2))-Pt*sin(q5+q3+q1);
p_y=q7+lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))-lf*sin(q3+q1-(pi/2))-Pt*cos(q5+q3+q1);


%velocity calculaiton

pt_x1d=jacobian(pt_x1,q)*qd;
pt_y1d=jacobian(pt_y1,q)*qd;

pf_x1d=jacobian(pf_x1,q)*qd;
pf_y1d=jacobian(pf_y1,q)*qd;

pT_xd=jacobian(pT_x,q)*qd;
pT_yd=jacobian(pT_y,q)*qd;

pt_x2d=jacobian(pt_x2,q)*qd;
pt_y2d=jacobian(p_y,q)*qd;

pf_x2d=jacobian(pf_x2,q)*qd;
pf_y2d=jacobian(pf_y2,q)*qd;

%theta abs
th_T=q1;
th_f1=q1+q2;
th_t1=q1+q2+q4;
th_f2=q1+q3;
th_t2=q1+q3+q5;

%omega abs
th_Td=jacobian(th_T,q)*qd;
th_f1d=jacobian(th_f1,q)*qd;
th_f2d=jacobian(th_f2,q)*qd;
th_t1d=jacobian(th_t1,q)*qd;
th_t2d=jacobian(th_t2,q)*qd;

%potential energy calculation

V= (MT*pT_y + Mf*(pf_y1+pf_y2) + Mt*(pt_y1+p_y))*g;

%kinetic energy calculation 

ke_t1= (Mt/2)*((pt_x1d)'*(pt_x1d)+ (pt_y1d)'*(pt_y1d))+(It/2)*(th_t1d'*th_t1d);
ke_t2= (Mt/2)*((pt_x2d)'*(pt_x2d)+ (pt_y2d)'*(pt_y2d))+(It/2)*(th_t2d'*th_t2d);
ke_f1= (Mf/2)*((pf_x1d)'*(pf_x1d)+ (pf_y1d)'*(pf_y1d))+(If/2)*(th_f1d'*th_f1d);
ke_f2= (Mf/2)*((pf_x2d)'*(pf_x2d)+ (pf_y2d)'*(pf_y2d))+(If/2)*(th_f2d'*th_f2d);
ke_T= (MT/2)*((pT_xd)'*(pT_xd)+ (pT_yd)'*(pT_yd))+(IT/2)*(th_Td'*th_Td);

K= ke_t1+ke_t2+ke_T+ke_f1+ke_f2;

D = jacobian(jacobian(K,qd),qd)

%gamma 2 - end points of the impact leg just before impact

p2_x=q6-lt*cos(q2+q4+q1-(pi/2))-lf*cos(q2+q1-(pi/2))+lf*cos(q3+q1-(pi/2))-lt*sin(q5+q3+q1);
p2_y=q7+lt*sin(q2+q4+q1-(pi/2))+lf*sin(q2+q1-(pi/2))-lf*sin(q3+q1-(pi/2))-lt*cos(q5+q3+q1);

p2=[p2_x;p2_y];

E2=jacobian(p2,q)





