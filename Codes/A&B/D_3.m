% NOTE THAT ::
% To run the code you need to add XSteam.m to the same folder
% that contain this m.file code


clc;clear,format long

P_c=0.1;
P_b=150;
P_R=30.511;
P_O=5.93;
P_c1=30.98;
P_c2=1.07;
T_max=500;
eta_T=0.85;
eta_P=0.95;

s=0.01;P_c3=(P_c+s):s:(P_b-s);c=length(P_c3);

h_1=XSteam('h_px',P_c,0);
s_1=XSteam('s_ph',P_c,h_1);

h_2s=XSteam('h_ps',P_O,s_1);
h_2=h_1+(h_2s-h_1)/eta_P;

h_4=XSteam('h_px',P_c2,0);
h_5=h_4;
T_4=XSteam('t_ph',P_c2,h_4);

T_3=T_4;
h_3=XSteam('h_pt',P_O,T_3);

h_6=XSteam('h_px',P_O,0);
s_6=XSteam('s_ph',P_O,h_6);

h_7s=XSteam('h_ps',P_b,s_6);
h_7=h_6+(h_7s-h_6)/eta_P;

h_12=XSteam('h_px',P_c1,0);
h_13=h_12;
T_12=XSteam('t_ph',P_c1,h_12);

T_11=T_12;
h_11=XSteam('h_pt',P_b,T_11);

h_14=XSteam('h_pt',P_b,T_max);
s_14=XSteam('s_ph',P_b,h_14);

h_17=XSteam('h_pt',P_R,T_max);
s_17=XSteam('s_ph',P_R,h_17);

for i=1:c

    h_9(i)=XSteam('h_px',P_c3(i),0);
    h_10(i)=h_9(i);
    T_9(i)=XSteam('t_ph',P_c3(i),h_9(i));
    T_8(i)=T_9(i);   

    if P_c3(i)<P_c2

        h_8(i)=XSteam('h_pt',P_O,T_8(i));

        h_15s(i)=XSteam('h_ps',P_c1,s_14);
        h_15(i)=h_14-(h_14-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_c1,h_15(i));

        h_16s(i)=XSteam('h_ps',P_R,s_15(i));
        h_16(i)=h_15(i)-(h_15(i)-h_16s(i))*eta_T;

        h_19s(i)=XSteam('h_ps',P_O,s_17);
        h_19(i)=h_17-(h_17-h_19s(i))*eta_T;
        s_19(i)=XSteam('s_ph',P_O,h_19(i));

        h_20s(i)=XSteam('h_ps',P_c2,s_19(i));
        h_20(i)=h_19(i)-(h_19(i)-h_20s(i))*eta_T;
        s_20(i)=XSteam('s_ph',P_c2,h_20(i));

        h_18s(i)=XSteam('h_ps',P_c3(i),s_20(i));
        h_18(i)=h_20(i)-(h_20(i)-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c3(i),h_18(i));

        h_21s(i)=XSteam('h_ps',P_c,s_18(i));
        h_21(i)=h_18(i)-(h_18(i)-h_21s(i))*eta_T;

        m_1(i)=(h_11-h_7)/(h_15(i)-h_12);
        m_3(i)=(h_6-h_3+m_1(i)*(h_3-h_13))/(h_19(i)-h_3);
        m_4(i)=((1-m_1(i)-m_3(i))*(h_3-h_8(i)))/(h_20(i)-h_4);
        m_2(i)=((1-m_1(i)-m_3(i))*(h_8(i)-h_2)+m_4(i)*(h_9(i)-h_5))/(h_18(i)-h_9(i));

        Q_dot_add(i)=h_14-h_11+(h_17-h_16(i))*(1-m_1(i));
        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i))*h_21(i)+(m_2(i)+m_4(i))*h_10(i)-h_1*(1-m_1(i)-m_3(i));
        
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    elseif P_c3(i)<P_O

        h_8(i)=XSteam('h_pt',P_O,T_8(i));

        h_15s(i)=XSteam('h_ps',P_c1,s_14);
        h_15(i)=h_14-(h_14-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_c1,h_15(i));

        h_16s(i)=XSteam('h_ps',P_R,s_15(i));
        h_16(i)=h_15(i)-(h_15(i)-h_16s(i))*eta_T;

        h_19s(i)=XSteam('h_ps',P_O,s_17);
        h_19(i)=h_17-(h_17-h_19s(i))*eta_T;
        s_19(i)=XSteam('s_ph',P_O,h_19(i));

        h_18s(i)=XSteam('h_ps',P_c3(i),s_19(i));
        h_18(i)=h_19(i)-(h_19(i)-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c3(i),h_18(i));

        h_20s(i)=XSteam('h_ps',P_c2,s_18(i));
        h_20(i)=h_18(i)-(h_18(i)-h_20s(i))*eta_T;
        s_20(i)=XSteam('s_ph',P_c2,h_20(i));

        h_21s(i)=XSteam('h_ps',P_c,s_20(i));
        h_21(i)=h_20(i)-(h_20(i)-h_21s(i))*eta_T;

        m_1(i)=(h_11-h_7)/(h_15(i)-h_12);
        m_3(i)=(h_6-h_8(i)+m_1(i)*(h_8(i)-h_13))/(h_19(i)-h_8(i));
        m_2(i)=((1-m_1(i)-m_3(i))*(h_8(i)-h_3))/(h_18(i)-h_9(i));
        m_4(i)=((1-m_1(i)-m_3(i))*(h_3-h_2)+m_2(i)*(h_4-h_10(i)))/(h_20(i)-h_4);

        Q_dot_add(i)=h_14-h_11+(h_17-h_16(i))*(1-m_1(i));
        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i))*h_21(i)+(m_2(i)+m_4(i))*h_5-h_1*(1-m_1(i)-m_3(i));
        
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    elseif P_c3(i)<P_R

        h_8(i)=XSteam('h_pt',P_b,T_8(i));

        h_15s(i)=XSteam('h_ps',P_c1,s_14);
        h_15(i)=h_14-(h_14-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_c1,h_15(i));

        h_16s(i)=XSteam('h_ps',P_R,s_15(i));
        h_16(i)=h_15(i)-(h_15(i)-h_16s(i))*eta_T;

        h_18s(i)=XSteam('h_ps',P_c3(i),s_17);
        h_18(i)=h_17-(h_17-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c3(i),h_18(i));

        h_19s(i)=XSteam('h_ps',P_O,s_18(i));
        h_19(i)=h_18(i)-(h_18(i)-h_19s(i))*eta_T;
        s_19(i)=XSteam('s_ph',P_O,h_19(i));

        h_20s(i)=XSteam('h_ps',P_c2,s_19(i));
        h_20(i)=h_19(i)-(h_19(i)-h_20s(i))*eta_T;
        s_20(i)=XSteam('s_ph',P_c2,h_20(i));

        h_21s(i)=XSteam('h_ps',P_c,s_20(i));
        h_21(i)=h_20(i)-(h_20(i)-h_21s(i))*eta_T;

        m_1(i)=(h_11-h_8(i))/(h_15(i)-h_12);
        m_2(i)=(h_8(i)-h_7+m_1(i)*(h_9(i)-h_13))/(h_18(i)-h_9(i));
        m_3(i)=(h_6-h_10(i)*(m_1(i)+m_2(i))-h_3*(1-m_1(i)-m_2(i)))/(h_19(i)-h_3);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i))*(h_3-h_2))/(h_20(i)-h_4);

        Q_dot_add(i)=h_14-h_11+(h_17-h_16(i))*(1-m_1(i));
        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i))*h_21(i)+m_4(i)*h_5-h_1*(1-m_1(i)-m_2(i)-m_3(i));
        
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    elseif P_c3(i)<P_c1

        h_8(i)=XSteam('h_pt',P_b,T_8(i));

        h_15s(i)=XSteam('h_ps',P_c1,s_14);
        h_15(i)=h_14-(h_14-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_c1,h_15(i));

        h_18s(i)=XSteam('h_ps',P_c3(i),s_15(i));
        h_18(i)=h_15(i)-(h_15(i)-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c3(i),h_18(i));

        h_16s(i)=XSteam('h_ps',P_R,s_18(i));
        h_16(i)=h_18(i)-(h_18(i)-h_16s(i))*eta_T;

        h_19s(i)=XSteam('h_ps',P_O,s_17);
        h_19(i)=h_17-(h_17-h_19s(i))*eta_T;
        s_19(i)=XSteam('s_ph',P_O,h_19(i));

        h_20s(i)=XSteam('h_ps',P_c2,s_19(i));
        h_20(i)=h_19(i)-(h_19(i)-h_20s(i))*eta_T;
        s_20(i)=XSteam('s_ph',P_c2,h_20(i));

        h_21s(i)=XSteam('h_ps',P_c,s_20(i));
        h_21(i)=h_20(i)-(h_20(i)-h_21s(i))*eta_T;

        m_1(i)=(h_11-h_8(i))/(h_15(i)-h_12);
        m_2(i)=(h_8(i)-h_7+m_1(i)*(h_9(i)-h_13))/(h_18(i)-h_9(i));
        m_3(i)=(h_6-h_10(i)*(m_1(i)+m_2(i))-h_3*(1-m_1(i)-m_2(i)))/(h_19(i)-h_3);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i))*(h_3-h_2))/(h_20(i)-h_4);

        Q_dot_add(i)=h_14-h_11+(h_17-h_16(i))*(1-m_1(i)-m_2(i));
        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i))*h_21(i)+m_4(i)*h_5-h_1*(1-m_1(i)-m_2(i)-m_3(i));
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
        
    else 

        h_8(i)=XSteam('h_pt',P_b,T_8(i));

        h_18s(i)=XSteam('h_ps',P_c3(i),s_14);
        h_18(i)=h_14-(h_14-h_18s(i))*eta_T;
        s_18(i)=XSteam('s_ph',P_c3(i),h_18(i));

        h_15s(i)=XSteam('h_ps',P_c1,s_18(i));
        h_15(i)=h_18(i)-(h_18(i)-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_c1,h_15(i));

        h_16s(i)=XSteam('h_ps',P_R,s_15(i));
        h_16(i)=h_15(i)-(h_15(i)-h_16s(i))*eta_T;

        h_19s(i)=XSteam('h_ps',P_O,s_17);
        h_19(i)=h_17-(h_17-h_19s(i))*eta_T;
        s_19(i)=XSteam('s_ph',P_O,h_19(i));

        h_20s(i)=XSteam('h_ps',P_c2,s_19(i));
        h_20(i)=h_19(i)-(h_19(i)-h_20s(i))*eta_T;
        s_20(i)=XSteam('s_ph',P_c2,h_20(i));

        h_21s(i)=XSteam('h_ps',P_c,s_20(i));
        h_21(i)=h_20(i)-(h_20(i)-h_21s(i))*eta_T;

        m_2(i)=(h_8(i)-h_11)/(h_18(i)-h_9(i));
        m_1(i)=(h_11-h_7+m_2(i)*(h_12-h_10(i)))/(h_15(i)-h_12);
        m_3(i)=(h_6-h_13*(m_1(i)+m_2(i))-h_3*(1-m_1(i)-m_2(i)))/(h_19(i)-h_3);
        m_4(i)=((1-m_1(i)-m_2(i)-m_3(i))*(h_3-h_2))/(h_20(i)-h_4);

        Q_dot_add(i)=h_14-h_8(i)+(h_17-h_16(i))*(1-m_1(i)-m_2(i));
        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i)-m_4(i))*h_21(i)+m_4(i)*h_5-h_1*(1-m_1(i)-m_2(i)-m_3(i));
        
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
    end
end


[eta_max,n_2]=max(eta_th);
[W_max,n_1]=max(W_dot_net);


P_c3_eta=P_c3(n_2);
P_c3_W=P_c3(n_1);

disp('[Efficiency]')
fprintf('\nThe optimum CFWH_3 Pressure that maximizes the efficiency is : \n(%f)\n',P_c3_eta)
fprintf('\nThe maximum efficiency is : \n(%f)\n\n\n',eta_max)

disp('[Work Net]')
fprintf('\nThe optimum CFWH_3 Pressure that maximizes the net work is : \n(%f)\n',P_c3_W)
fprintf('\nThe maximum Work Net is : \n(%f)\n',W_max)


figure(1)
% yyaxis left
plot(P_c3,eta_th)
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-3} VS Pressure')
% yyaxis right
% plot(P_c3,W_dot_net)
grid on

