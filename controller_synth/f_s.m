function f_sym = f_s(in1,in2,in3)
%F_S
%    F_SYM = F_S(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    09-Jun-2024 05:23:16

beta = in1(4,:);
delta = in2(2,:);
e_psi = in1(2,:);
e_y = in1(1,:);
k_ref = in3(1,:);
omega = in1(5,:);
u1 = in2(1,:);
v = in1(3,:);
t2 = cos(beta);
t3 = cos(e_psi);
t4 = sin(beta);
t5 = sin(e_psi);
t6 = beta+delta;
t7 = e_y.*k_ref;
t8 = v.^4;
t12 = -e_psi;
t13 = 1.0./v;
t20 = u1.*1.0e+4;
t22 = delta.*(2.74e+2./2.5e+1);
t23 = delta.*(4.11e+2./2.5e+1);
t28 = omega.*1.19016;
t29 = omega.*1.37484;
t31 = v.*1.127945952e+1;
t33 = v.*1.302972048e+1;
t9 = cos(t6);
t10 = sin(t6);
t11 = t4.*v;
t14 = 1.0./t2;
t15 = t2.*t3.*v;
t17 = beta+t12;
t18 = t7-1.0;
t24 = -t20;
t25 = -t23;
t35 = t8.*1.691918928;
t36 = t8.*1.954458072;
t16 = t5.*t11;
t19 = -t11;
t21 = cos(t17);
t32 = t11+t29;
t52 = t33+t36+5.863374216e+1;
t53 = t24+t31+t35+5.075756784e+1;
t26 = 1.0./t21;
t27 = t15+t16;
t34 = t19+t28;
t37 = t13.*t14.*t32;
t30 = 1.0./t27;
t38 = atan(t37);
t39 = t13.*t14.*t34;
t40 = atan(t39);
t43 = t38.*1.267e+1;
t44 = t38.*1.9005e+1;
t41 = t40.*(2.74e+2./2.5e+1);
t42 = t40.*(4.11e+2./2.5e+1);
t45 = atan(t43);
t47 = -t44;
t46 = -t41;
t48 = t45./2.0;
t49 = t22+t46;
t54 = t47+t48;
t50 = atan(t49);
t55 = atan(t54);
t51 = t50./2.0;
t56 = t55.*(1.3e+1./1.0e+1);
t57 = sin(t56);
t58 = t25+t42+t51;
t59 = atan(t58);
t60 = t59.*(1.3e+1./1.0e+1);
t61 = sin(t60);
mt1 = [t18.*t26.*sin(t17);-k_ref-omega.*t13.*t18.*t26;t18.*t30.*((t2.*t53)./1.239e+3-t4.*t57.*3.186287328490718+(t9.*t52)./1.239e+3-t10.*t61.*3.680710250201776);-t18.*t30.*(omega+(t13.*(t4.*t53+t2.*t57.*3.94781e+3+t10.*t52+t9.*t61.*4.5604e+3))./1.239e+3)];
mt2 = [t18.*t30.*(t57.*(-3.097949258219178)+t61.*cos(delta).*3.097948438356164+(sin(delta).*(t8.*2.32611781897152+v.*1.55074521264768e+1+6.97835345691456e+1))./1.752e+3);-t13.*t18.*t26];
f_sym = [mt1;mt2];
end
