function J = costFunction(in1)
%costFunction
%    J = costFunction(IN1)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    09-Jun-2024 05:23:35

q_arr1_1 = in1(3,:);
q_arr1_2 = in1(11,:);
q_arr1_3 = in1(19,:);
q_arr1_4 = in1(27,:);
q_arr1_5 = in1(35,:);
q_arr1_6 = in1(43,:);
q_arr1_7 = in1(51,:);
q_arr1_8 = in1(59,:);
q_arr1_9 = in1(67,:);
q_arr2_1 = in1(4,:);
q_arr2_2 = in1(12,:);
q_arr2_3 = in1(20,:);
q_arr2_4 = in1(28,:);
q_arr2_5 = in1(36,:);
q_arr2_6 = in1(44,:);
q_arr2_7 = in1(52,:);
q_arr2_8 = in1(60,:);
q_arr2_9 = in1(68,:);
q_arr3_1 = in1(5,:);
q_arr3_2 = in1(13,:);
q_arr3_3 = in1(21,:);
q_arr3_4 = in1(29,:);
q_arr3_5 = in1(37,:);
q_arr3_6 = in1(45,:);
q_arr3_7 = in1(53,:);
q_arr3_8 = in1(61,:);
q_arr3_9 = in1(69,:);
q_arr4_1 = in1(6,:);
q_arr4_2 = in1(14,:);
q_arr4_3 = in1(22,:);
q_arr4_4 = in1(30,:);
q_arr4_5 = in1(38,:);
q_arr4_6 = in1(46,:);
q_arr4_7 = in1(54,:);
q_arr4_8 = in1(62,:);
q_arr4_9 = in1(70,:);
q_arr5_1 = in1(7,:);
q_arr5_2 = in1(15,:);
q_arr5_3 = in1(23,:);
q_arr5_4 = in1(31,:);
q_arr5_5 = in1(39,:);
q_arr5_6 = in1(47,:);
q_arr5_7 = in1(55,:);
q_arr5_8 = in1(63,:);
q_arr5_9 = in1(71,:);
q_arr6_1 = in1(8,:);
q_arr6_2 = in1(16,:);
q_arr6_3 = in1(24,:);
q_arr6_4 = in1(32,:);
q_arr6_5 = in1(40,:);
q_arr6_6 = in1(48,:);
q_arr6_7 = in1(56,:);
q_arr6_8 = in1(64,:);
q_arr6_9 = in1(72,:);
q_arr1_10 = in1(75,:);
q_arr1_11 = in1(83,:);
q_arr1_12 = in1(91,:);
q_arr1_13 = in1(99,:);
q_arr1_14 = in1(107,:);
q_arr1_15 = in1(115,:);
q_arr1_16 = in1(123,:);
q_arr1_17 = in1(131,:);
q_arr1_18 = in1(139,:);
q_arr1_19 = in1(147,:);
q_arr2_10 = in1(76,:);
q_arr2_11 = in1(84,:);
q_arr2_12 = in1(92,:);
q_arr2_13 = in1(100,:);
q_arr2_14 = in1(108,:);
q_arr2_15 = in1(116,:);
q_arr2_16 = in1(124,:);
q_arr2_17 = in1(132,:);
q_arr2_18 = in1(140,:);
q_arr2_19 = in1(148,:);
q_arr3_10 = in1(77,:);
q_arr3_11 = in1(85,:);
q_arr3_12 = in1(93,:);
q_arr3_13 = in1(101,:);
q_arr3_14 = in1(109,:);
q_arr3_15 = in1(117,:);
q_arr3_16 = in1(125,:);
q_arr3_17 = in1(133,:);
q_arr3_18 = in1(141,:);
q_arr3_19 = in1(149,:);
q_arr4_10 = in1(78,:);
q_arr4_11 = in1(86,:);
q_arr4_12 = in1(94,:);
q_arr4_13 = in1(102,:);
q_arr4_14 = in1(110,:);
q_arr4_15 = in1(118,:);
q_arr4_16 = in1(126,:);
q_arr4_17 = in1(134,:);
q_arr4_18 = in1(142,:);
q_arr4_19 = in1(150,:);
q_arr5_10 = in1(79,:);
q_arr5_11 = in1(87,:);
q_arr5_12 = in1(95,:);
q_arr5_13 = in1(103,:);
q_arr5_14 = in1(111,:);
q_arr5_15 = in1(119,:);
q_arr5_16 = in1(127,:);
q_arr5_17 = in1(135,:);
q_arr5_18 = in1(143,:);
q_arr5_19 = in1(151,:);
q_arr6_10 = in1(80,:);
q_arr6_11 = in1(88,:);
q_arr6_12 = in1(96,:);
q_arr6_13 = in1(104,:);
q_arr6_14 = in1(112,:);
q_arr6_15 = in1(120,:);
q_arr6_16 = in1(128,:);
q_arr6_17 = in1(136,:);
q_arr6_18 = in1(144,:);
q_arr6_19 = in1(152,:);
u_arr2_1 = in1(2,:);
u_arr2_2 = in1(10,:);
u_arr2_3 = in1(18,:);
u_arr2_4 = in1(26,:);
u_arr2_5 = in1(34,:);
u_arr2_6 = in1(42,:);
u_arr2_7 = in1(50,:);
u_arr2_8 = in1(58,:);
u_arr2_9 = in1(66,:);
u_arr2_10 = in1(74,:);
u_arr2_11 = in1(82,:);
u_arr2_12 = in1(90,:);
u_arr2_13 = in1(98,:);
u_arr2_14 = in1(106,:);
u_arr2_15 = in1(114,:);
u_arr2_16 = in1(122,:);
u_arr2_17 = in1(130,:);
u_arr2_18 = in1(138,:);
u_arr2_19 = in1(146,:);
t2 = -q_arr4_1;
t3 = -q_arr4_2;
t4 = -q_arr4_3;
t5 = -q_arr4_4;
t6 = -q_arr4_5;
t7 = -q_arr4_6;
t8 = -q_arr4_7;
t9 = -q_arr4_8;
t10 = -q_arr4_9;
t11 = -q_arr4_10;
t12 = -q_arr4_11;
t13 = -q_arr4_12;
t14 = -q_arr4_13;
t15 = -q_arr4_14;
t16 = -q_arr4_15;
t17 = -q_arr4_16;
t18 = -q_arr4_17;
t19 = -q_arr4_18;
t20 = -q_arr4_19;
et1 = q_arr6_19.*1.0e+3+(q_arr3_1./1.0e+4-1.0./1.0e+3).*(q_arr3_1-1.0e+1)+(q_arr3_2./1.0e+4-1.0./1.0e+3).*(q_arr3_2-1.0e+1)+(q_arr3_3./1.0e+4-1.0./1.0e+3).*(q_arr3_3-1.0e+1)+(q_arr3_4./1.0e+4-1.0./1.0e+3).*(q_arr3_4-1.0e+1)+(q_arr3_5./1.0e+4-1.0./1.0e+3).*(q_arr3_5-1.0e+1)+(q_arr3_6./1.0e+4-1.0./1.0e+3).*(q_arr3_6-1.0e+1)+(q_arr3_7./1.0e+4-1.0./1.0e+3).*(q_arr3_7-1.0e+1)+(q_arr3_8./1.0e+4-1.0./1.0e+3).*(q_arr3_8-1.0e+1);
et2 = (q_arr3_9./1.0e+4-1.0./1.0e+3).*(q_arr3_9-1.0e+1)+(q_arr3_10./1.0e+4-1.0./1.0e+3).*(q_arr3_10-1.0e+1)+(q_arr3_11./1.0e+4-1.0./1.0e+3).*(q_arr3_11-1.0e+1)+(q_arr3_12./1.0e+4-1.0./1.0e+3).*(q_arr3_12-1.0e+1)+(q_arr3_13./1.0e+4-1.0./1.0e+3).*(q_arr3_13-1.0e+1)+(q_arr3_14./1.0e+4-1.0./1.0e+3).*(q_arr3_14-1.0e+1)+(q_arr3_15./1.0e+4-1.0./1.0e+3).*(q_arr3_15-1.0e+1)+(q_arr3_16./1.0e+4-1.0./1.0e+3).*(q_arr3_16-1.0e+1);
et3 = (q_arr3_17./1.0e+4-1.0./1.0e+3).*(q_arr3_17-1.0e+1)+(q_arr3_18./1.0e+4-1.0./1.0e+3).*(q_arr3_18-1.0e+1)+(q_arr3_19./1.0e+4-1.0./1.0e+3).*(q_arr3_19-1.0e+1)+t2.*(q_arr2_1-q_arr4_1.*(1.1e+1./1.0e+1))+t3.*(q_arr2_2-q_arr4_2.*(1.1e+1./1.0e+1))+t4.*(q_arr2_3-q_arr4_3.*(1.1e+1./1.0e+1))+t5.*(q_arr2_4-q_arr4_4.*(1.1e+1./1.0e+1))+t6.*(q_arr2_5-q_arr4_5.*(1.1e+1./1.0e+1))+t7.*(q_arr2_6-q_arr4_6.*(1.1e+1./1.0e+1))+t8.*(q_arr2_7-q_arr4_7.*(1.1e+1./1.0e+1))+t9.*(q_arr2_8-q_arr4_8.*(1.1e+1./1.0e+1))+t10.*(q_arr2_9-q_arr4_9.*(1.1e+1./1.0e+1))+t11.*(q_arr2_10-q_arr4_10.*(1.1e+1./1.0e+1));
et4 = t12.*(q_arr2_11-q_arr4_11.*(1.1e+1./1.0e+1))+t13.*(q_arr2_12-q_arr4_12.*(1.1e+1./1.0e+1))+t14.*(q_arr2_13-q_arr4_13.*(1.1e+1./1.0e+1))+t15.*(q_arr2_14-q_arr4_14.*(1.1e+1./1.0e+1))+t16.*(q_arr2_15-q_arr4_15.*(1.1e+1./1.0e+1))+t17.*(q_arr2_16-q_arr4_16.*(1.1e+1./1.0e+1))+t18.*(q_arr2_17-q_arr4_17.*(1.1e+1./1.0e+1))+t19.*(q_arr2_18-q_arr4_18.*(1.1e+1./1.0e+1))+t20.*(q_arr2_19-q_arr4_19.*(1.1e+1./1.0e+1))+q_arr2_1.*(q_arr2_1.*1.001+t2)+q_arr2_2.*(q_arr2_2.*1.001+t3)+q_arr2_3.*(q_arr2_3.*1.001+t4)+q_arr2_4.*(q_arr2_4.*1.001+t5)+q_arr2_5.*(q_arr2_5.*1.001+t6)+q_arr2_6.*(q_arr2_6.*1.001+t7);
et5 = q_arr2_7.*(q_arr2_7.*1.001+t8)+q_arr2_8.*(q_arr2_8.*1.001+t9)+q_arr2_9.*(q_arr2_9.*1.001+t10)+q_arr2_10.*(q_arr2_10.*1.001+t11)+q_arr2_11.*(q_arr2_11.*1.001+t12)+q_arr2_12.*(q_arr2_12.*1.001+t13)+q_arr2_13.*(q_arr2_13.*1.001+t14)+q_arr2_14.*(q_arr2_14.*1.001+t15)+q_arr2_15.*(q_arr2_15.*1.001+t16)+q_arr2_16.*(q_arr2_16.*1.001+t17)+q_arr2_17.*(q_arr2_17.*1.001+t18)+q_arr2_18.*(q_arr2_18.*1.001+t19)+q_arr2_19.*(q_arr2_19.*1.001+t20);
et6 = q_arr1_1.^2.*5.0e+1+q_arr1_2.^2.*5.0e+1+q_arr1_3.^2.*5.0e+1+q_arr1_4.^2.*5.0e+1+q_arr1_5.^2.*5.0e+1+q_arr1_6.^2.*5.0e+1+q_arr1_7.^2.*5.0e+1+q_arr1_8.^2.*5.0e+1+q_arr1_9.^2.*5.0e+1+q_arr5_1.^2.*1.0e+1+q_arr5_2.^2.*1.0e+1+q_arr5_3.^2.*1.0e+1+q_arr5_4.^2.*1.0e+1+q_arr5_5.^2.*1.0e+1+q_arr5_6.^2.*1.0e+1+q_arr5_7.^2.*1.0e+1+q_arr5_8.^2.*1.0e+1+q_arr5_9.^2.*1.0e+1+q_arr6_1.^2./1.0e+2+q_arr6_2.^2./1.0e+2+q_arr6_3.^2./1.0e+2+q_arr6_4.^2./1.0e+2+q_arr6_5.^2./1.0e+2+q_arr6_6.^2./1.0e+2+q_arr6_7.^2./1.0e+2+q_arr6_8.^2./1.0e+2+q_arr6_9.^2./1.0e+2+q_arr1_10.^2.*5.0e+1;
et7 = q_arr1_11.^2.*5.0e+1+q_arr1_12.^2.*5.0e+1+q_arr1_13.^2.*5.0e+1+q_arr1_14.^2.*5.0e+1+q_arr1_15.^2.*5.0e+1+q_arr1_16.^2.*5.0e+1+q_arr1_17.^2.*5.0e+1+q_arr1_18.^2.*5.0e+1+q_arr1_19.^2.*1.0005e+5+q_arr5_10.^2.*1.0e+1+q_arr5_11.^2.*1.0e+1+q_arr5_12.^2.*1.0e+1+q_arr5_13.^2.*1.0e+1+q_arr5_14.^2.*1.0e+1+q_arr5_15.^2.*1.0e+1+q_arr5_16.^2.*1.0e+1+q_arr5_17.^2.*1.0e+1+q_arr5_18.^2.*1.0e+1+q_arr5_19.^2.*1.0e+1+q_arr6_10.^2./1.0e+2+q_arr6_11.^2./1.0e+2+q_arr6_12.^2./1.0e+2+q_arr6_13.^2./1.0e+2+q_arr6_14.^2./1.0e+2+q_arr6_15.^2./1.0e+2+q_arr6_16.^2./1.0e+2+q_arr6_17.^2./1.0e+2+q_arr6_18.^2./1.0e+2;
et8 = q_arr6_19.^2./1.0e+2+u_arr2_1.^2./1.0e+1+u_arr2_2.^2./1.0e+1+u_arr2_3.^2./1.0e+1+u_arr2_4.^2./1.0e+1+u_arr2_5.^2./1.0e+1+u_arr2_6.^2./1.0e+1+u_arr2_7.^2./1.0e+1+u_arr2_8.^2./1.0e+1+u_arr2_9.^2./1.0e+1+u_arr2_10.^2./1.0e+1+u_arr2_11.^2./1.0e+1+u_arr2_12.^2./1.0e+1+u_arr2_13.^2./1.0e+1+u_arr2_14.^2./1.0e+1+u_arr2_15.^2./1.0e+1+u_arr2_16.^2./1.0e+1+u_arr2_17.^2./1.0e+1+u_arr2_18.^2./1.0e+1+u_arr2_19.^2./1.0e+1;
J = et1+et2+et3+et4+et5+et6+et7+et8;
end
