% NOTE THAT ::
% To run the code you need to add XSteam.m to the same folder
% that contain this m.file code


clc;clear,format long

P_c=0.1;
P_b=150;
P_R=30.54;
P_O=5.93;
P_c1=30.98;
P_c2=1.07;
P_c3=13.04;
T_max=500;
eta_T=0.85;
eta_P=0.95;

s=0.01;P_c4=(P_c+s):s:(P_b-s);c=length(P_c4);

h_1=XSteam('h_px',P_c,0);
s_1=XSteam('s_ph',P_c,h_1);

h_2s=XSteam('h_ps',P_O,s_1);
h_2=h_1+(h_2s-h_1)/eta_P;

h_7=XSteam('h_px',P_c2,0);
h_8=h_7;
T_7=XSteam('t_ph',P_c2,h_7);

T_6=T_7;
h_6=XSteam('h_pt',P_O,T_6);

h_9=XSteam('h_px',P_O,0);
s_9=XSteam('s_ph',P_O,h_9);

h_10s=XSteam('h_ps',P_b,s_9);
h_10=h_9+(h_10s-h_9)/eta_P;

h_12=XSteam('h_px',P_c3,0);
h_13=h_12;
T_12=XSteam('t_ph',P_c3,h_12);

T_11=T_12;
h_11=XSteam('h_pt',P_b,T_11);

h_15=XSteam('h_px',P_c1,0);
h_16=h_15;
T_15=XSteam('t_ph',P_c1,h_15);

T_14=T_15;
h_14=XSteam('h_pt',P_b,T_14);

h_17=XSteam('h_pt',P_b,T_max);
s_17=XSteam('s_ph',P_b,h_17);

h_20=XSteam('h_pt',P_R,T_max);
s_20=XSteam('s_ph',P_R,h_20);

for i=1:c

    h_4(i)=XSteam('h_px',P_c4(i),0);
    h_5(i)=h_4(i);
    T_4(i)=XSteam('t_ph',P_c4(i),h_4(i));

    T_3(i)=T_4(i);

    if P_c4(i)<P_c2

        h_3(i)=XSteam('h_pt',P_O,T_3(i));

        h_18s(i)=XSteam('h_ps',P_c1,s_17);
        h_18(i)=h_17-(h_17-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c1,h_18(i));

        h_19s(i)=XSteam('h_ps',P_R,s_18(i));
        h_19(i)=h_18(i)-(h_18(i)-h_19s(i))*eta_T;

        h_21s(i)=XSteam('h_ps',P_c3,s_20);
        h_21(i)=h_20-(h_20-h_21s(i))*eta_T;
        s_21(i)=XSteam('s_ph',P_c3,h_21(i));

        h_22s(i)=XSteam('h_ps',P_O,s_21(i));
        h_22(i)=h_21(i)-(h_21(i)-h_22s(i))*eta_T;
        s_22(i)=XSteam('s_ph',P_O,h_22(i));

        h_23s(i)=XSteam('h_ps',P_c2,s_22(i));
        h_23(i)=h_22(i)-(h_22(i)-h_23s(i))*eta_T;
        s_23(i)=XSteam('s_ph',P_c2,h_23(i));

        h_24s(i)=XSteam('h_ps',P_c4(i),s_23(i));
        h_24(i)=h_23(i)-(h_23(i)-h_24s(i))*eta_T;
        s_24(i)=XSteam('s_ph',P_c4(i),h_24(i));

        h_25s(i)=XSteam('h_ps',P_c,s_24(i));
        h_25(i)=h_24(i)-(h_24(i)-h_25s(i))*eta_T;

        m_1(i)=(h_14-h_11)/(h_18(i)-h_15);
        m_2(i)=(h_11-h_10+m_1(i)*(h_12-h_16))/(h_21(i)-h_12);
        m_3(i)=(h_9-h_6+(m_1(i)+m_2(i))*(h_6-h_13))/(h_22(i)-h_6);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i))*(h_6-h_3(i)))/(h_23(i)-h_7);
        m_5(i)=(m_4(i)*(h_4(i)-h_8)+(1-m_1(i)-m_2(i)-m_3(i))*(h_3(i)-h_2))/(h_24(i)-h_4(i));

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i)-m_5(i))*h_25(i)+(m_4(i)+m_5(i))*h_5(i)-(1-m_1(i)-m_2(i)-m_3(i))*h_1;
        Q_dot_add(i)=h_17-h_14+(h_20-h_19(i))*(1-m_1(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    elseif P_c4(i)<P_O 

        h_3(i)=XSteam('h_pt',P_O,T_3(i));

        h_18s(i)=XSteam('h_ps',P_c1,s_17);
        h_18(i)=h_17-(h_17-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c1,h_18(i));

        h_19s(i)=XSteam('h_ps',P_R,s_18(i));
        h_19(i)=h_18(i)-(h_18(i)-h_19s(i))*eta_T;

        h_21s(i)=XSteam('h_ps',P_c3,s_20);
        h_21(i)=h_20-(h_20-h_21s(i))*eta_T;
        s_21(i)=XSteam('s_ph',P_c3,h_21(i));

        h_22s(i)=XSteam('h_ps',P_O,s_21(i));
        h_22(i)=h_21(i)-(h_21(i)-h_22s(i))*eta_T;
        s_22(i)=XSteam('s_ph',P_O,h_22(i));

        h_24s(i)=XSteam('h_ps',P_c4(i),s_22(i));
        h_24(i)=h_22(i)-(h_22(i)-h_24s(i))*eta_T;
        s_24(i)=XSteam('s_ph',P_c4(i),h_24(i));

        h_23s(i)=XSteam('h_ps',P_c2,s_24(i));
        h_23(i)=h_24(i)-(h_24(i)-h_23s(i))*eta_T;
        s_23(i)=XSteam('s_ph',P_c2,h_23(i));

        h_25s(i)=XSteam('h_ps',P_c,s_23(i));
        h_25(i)=h_23(i)-(h_23(i)-h_25s(i))*eta_T;

        m_1(i)=(h_14-h_11)/(h_18(i)-h_15);
        m_2(i)=(h_11-h_10+m_1(i)*(h_12-h_16))/(h_21(i)-h_12);
        m_3(i)=(h_9-h_3(i)+(m_1(i)+m_2(i))*(h_3(i)-h_13))/(h_22(i)-h_3(i));
        m_5(i)=((1-m_1(i)-m_2(i)-m_3(i))*(h_3(i)-h_6))/(h_24(i)-h_4(i));
        m_4(i)=(m_5(i)*(h_7-h_5(i))+(1-m_1(i)-m_2(i)-m_3(i))*(h_6-h_2))/(h_23(i)-h_7);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i)-m_5(i))*h_25(i)+(m_4(i)+m_5(i))*h_8-(1-m_1(i)-m_2(i)-m_3(i))*h_1;
        Q_dot_add(i)=h_17-h_14+(h_20-h_19(i))*(1-m_1(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
        
    elseif P_c4(i)<P_c3 

        h_3(i)=XSteam('h_pt',P_b,T_3(i));

        h_18s(i)=XSteam('h_ps',P_c1,s_17);
        h_18(i)=h_17-(h_17-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c1,h_18(i));

        h_19s(i)=XSteam('h_ps',P_R,s_18(i));
        h_19(i)=h_18(i)-(h_18(i)-h_19s(i))*eta_T;

        h_21s(i)=XSteam('h_ps',P_c3,s_20);
        h_21(i)=h_20-(h_20-h_21s(i))*eta_T;
        s_21(i)=XSteam('s_ph',P_c3,h_21(i));

        h_24s(i)=XSteam('h_ps',P_c4(i),s_21(i));
        h_24(i)=h_21(i)-(h_21(i)-h_24s(i))*eta_T;
        s_24(i)=XSteam('s_ph',P_c4(i),h_24(i));

        h_22s(i)=XSteam('h_ps',P_O,s_24(i));
        h_22(i)=h_24(i)-(h_24(i)-h_22s(i))*eta_T;
        s_22(i)=XSteam('s_ph',P_O,h_22(i));

        h_23s(i)=XSteam('h_ps',P_c2,s_22(i));
        h_23(i)=h_22(i)-(h_22(i)-h_23s(i))*eta_T;
        s_23(i)=XSteam('s_ph',P_c2,h_23(i));

        h_25s(i)=XSteam('h_ps',P_c,s_23(i));
        h_25(i)=h_23(i)-(h_23(i)-h_25s(i))*eta_T;

        m_1(i)=(h_14-h_11)/(h_18(i)-h_15);
        m_2(i)=(h_11-h_3(i)+m_1(i)*(h_12-h_16))/(h_21(i)-h_12);
        m_5(i)=(h_3(i)-h_10+(m_1(i)+m_2(i))*(h_4(i)-h_13))/(h_24(i)-h_4(i));
        m_3(i)=((m_1(i)+m_2(i)+m_5(i))*(h_6-h_5(i))+h_9-h_6)/(h_22(i)-h_6);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*(h_6-h_2))/(h_23(i)-h_7);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i)-m_5(i))*h_25(i)+m_4(i)*h_8-(1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*h_1;
        Q_dot_add(i)=h_17-h_14+(h_20-h_19(i))*(1-m_1(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));  

    elseif P_c4(i)<P_R 

        h_3(i)=XSteam('h_pt',P_b,T_3(i));

        h_18s(i)=XSteam('h_ps',P_c1,s_17);
        h_18(i)=h_17-(h_17-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c1,h_18(i));

        h_19s(i)=XSteam('h_ps',P_R,s_18(i));
        h_19(i)=h_18(i)-(h_18(i)-h_19s(i))*eta_T;

        h_24s(i)=XSteam('h_ps',P_c4(i),s_20);
        h_24(i)=h_20-(h_20-h_24s(i))*eta_T;
        s_24(i)=XSteam('s_ph',P_c4(i),h_24(i));

        h_21s(i)=XSteam('h_ps',P_c3,s_24(i));
        h_21(i)=h_24(i)-(h_24(i)-h_21s(i))*eta_T;
        s_21(i)=XSteam('s_ph',P_c3,h_21(i));

        h_22s(i)=XSteam('h_ps',P_O,s_21(i));
        h_22(i)=h_21(i)-(h_21(i)-h_22s(i))*eta_T;
        s_22(i)=XSteam('s_ph',P_O,h_22(i));

        h_23s(i)=XSteam('h_ps',P_c2,s_22(i));
        h_23(i)=h_22(i)-(h_22(i)-h_23s(i))*eta_T;
        s_23(i)=XSteam('s_ph',P_c2,h_23(i));

        h_25s(i)=XSteam('h_ps',P_c,s_23(i));
        h_25(i)=h_23(i)-(h_23(i)-h_25s(i))*eta_T;

        m_1(i)=(h_14-h_3(i))/(h_18(i)-h_15);
        m_5(i)=(h_3(i)-h_11+m_1(i)*(h_4(i)-h_16))/(h_24(i)-h_4(i));
        m_2(i)=(h_11-h_10+(m_1(i)+m_5(i))*(h_12-h_5(i)))/(h_21(i)-h_12);
        m_3(i)=((m_1(i)+m_2(i)+m_5(i))*(h_6-h_13)+h_9-h_6)/(h_22(i)-h_6);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*(h_6-h_2))/(h_23(i)-h_7);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i)-m_5(i))*h_25(i)+m_4(i)*h_8-(1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*h_1;
        Q_dot_add(i)=h_17-h_14+(h_20-h_19(i))*(1-m_1(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    elseif P_c4(i)<P_c1 

        h_3(i)=XSteam('h_pt',P_b,T_3(i));

        h_18s(i)=XSteam('h_ps',P_c1,s_17);
        h_18(i)=h_17-(h_17-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c1,h_18(i));

        h_24s(i)=XSteam('h_ps',P_c4(i),s_18(i));
        h_24(i)=h_18(i)-(h_18(i)-h_24s(i))*eta_T;
        s_24(i)=XSteam('s_ph',P_c4(i),h_24(i));

        h_19s(i)=XSteam('h_ps',P_R,s_24(i));
        h_19(i)=h_24(i)-(h_24(i)-h_19s(i))*eta_T;

        h_21s(i)=XSteam('h_ps',P_c3,s_20);
        h_21(i)=h_20-(h_20-h_21s(i))*eta_T;
        s_21(i)=XSteam('s_ph',P_c3,h_21(i));

        h_22s(i)=XSteam('h_ps',P_O,s_21(i));
        h_22(i)=h_21(i)-(h_21(i)-h_22s(i))*eta_T;
        s_22(i)=XSteam('s_ph',P_O,h_22(i));

        h_23s(i)=XSteam('h_ps',P_c2,s_22(i));
        h_23(i)=h_22(i)-(h_22(i)-h_23s(i))*eta_T;
        s_23(i)=XSteam('s_ph',P_c2,h_23(i));

        h_25s(i)=XSteam('h_ps',P_c,s_23(i));
        h_25(i)=h_23(i)-(h_23(i)-h_25s(i))*eta_T;

        m_1(i)=(h_14-h_3(i))/(h_18(i)-h_15);
        m_5(i)=(h_3(i)-h_11+m_1(i)*(h_4(i)-h_16))/(h_24(i)-h_4(i));
        m_2(i)=(h_11-h_10+(m_1(i)+m_5(i))*(h_12-h_5(i)))/(h_21(i)-h_12);
        m_3(i)=((m_1(i)+m_2(i)+m_5(i))*(h_6-h_13)+h_9-h_6)/(h_22(i)-h_6);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*(h_6-h_2))/(h_23(i)-h_7);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i)-m_5(i))*h_25(i)+m_4(i)*h_8-(1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*h_1;
        Q_dot_add(i)=h_17-h_14+(h_20-h_19(i))*(1-m_1(i)-m_5(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    else

        h_3(i)=XSteam('h_pt',P_b,T_3(i));

        h_24s(i)=XSteam('h_ps',P_c4(i),s_17);
        h_24(i)=h_17-(h_17-h_24s(i))*eta_T;
        s_24(i)=XSteam('s_ph',P_c4(i),h_24(i));

        h_18s(i)=XSteam('h_ps',P_c1,s_24(i));
        h_18(i)=h_24(i)-(h_24(i)-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c1,h_18(i));

        h_19s(i)=XSteam('h_ps',P_R,s_18(i));
        h_19(i)=h_18(i)-(h_18(i)-h_19s(i))*eta_T;

        h_21s(i)=XSteam('h_ps',P_c3,s_20);
        h_21(i)=h_20-(h_20-h_21s(i))*eta_T;
        s_21(i)=XSteam('s_ph',P_c3,h_21(i));

        h_22s(i)=XSteam('h_ps',P_O,s_21(i));
        h_22(i)=h_21(i)-(h_21(i)-h_22s(i))*eta_T;
        s_22(i)=XSteam('s_ph',P_O,h_22(i));

        h_23s(i)=XSteam('h_ps',P_c2,s_22(i));
        h_23(i)=h_22(i)-(h_22(i)-h_23s(i))*eta_T;
        s_23(i)=XSteam('s_ph',P_c2,h_23(i));

        h_25s(i)=XSteam('h_ps',P_c,s_23(i));
        h_25(i)=h_23(i)-(h_23(i)-h_25s(i))*eta_T;

        m_5(i)=(h_3(i)-h_14)/(h_24(i)-h_4(i));
        m_1(i)=(h_14-h_11+m_5(i)*(h_15-h_5(i)))/(h_18(i)-h_15);
        m_2(i)=(h_11-h_10+(m_1(i)+m_5(i))*(h_12-h_16))/(h_21(i)-h_12);
        m_3(i)=((m_1(i)+m_2(i)+m_5(i))*(h_6-h_13)+h_9-h_6)/(h_22(i)-h_6);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*(h_6-h_2))/(h_23(i)-h_7);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i)-m_5(i))*h_25(i)+m_4(i)*h_8-(1-m_1(i)-m_2(i)-m_3(i)-m_5(i))*h_1;
        Q_dot_add(i)=h_17-h_3(i)+(h_20-h_19(i))*(1-m_1(i)-m_5(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
    
    end
end

[eta_max,n]=max(eta_th);
[W_max,m]=max(W_dot_net);

P_c4_eta=P_c4(n);
P_c4_W=P_c4(m);

disp('[Efficiency]')
fprintf('\nThe optimum CFWH_4 Pressure that maximizes the efficiency is : \n(%f)\n',P_c4_eta)
fprintf('\nThe maximum efficiency is : \n(%f)\n\n\n',eta_max)

disp('[Work Net]')
fprintf('\nThe optimum CFWH_4 Pressure that maximizes the net work is : \n(%f)\n',P_c4_W)
fprintf('\nThe maximum Work Net is : \n(%f)\n',W_max)

figure(1)
% yyaxis left
plot(P_c4,eta_th)
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-4} VS Pressure')
% yyaxis right
% plot(P_c4,W_dot_net)
grid on


N=1:5;
eta_th_opt=[0.401182  0.409879  0.417650  0.421975  0.426192 ];

figure(2)
plot(N,eta_th_opt)
