% NOTE THAT ::
% To run the code you need to add XSteam.m to the same folder
% that contain this m.file code


clc;clear,format long

P_c=0.1;
P_b=150;
P_R=30.54;
P_O=5.93;
P_c1=30.98;
T_max=500;
eta_T=0.85;
eta_P=0.95;

s=0.01;P_c2=(P_c+s):s:(P_b-s);c=length(P_c2);

h_1=XSteam('h_px',P_c,0);
s_1=XSteam('s_ph',P_c,h_1);

h_2s=XSteam('h_ps',P_O,s_1);
h_2=h_1+(h_2s-h_1)/eta_P;

h_6=XSteam('h_px',P_O,0);
s_6=XSteam('s_ph',P_O,h_6);

h_7s=XSteam('h_ps',P_b,s_6);
h_7=h_6+(h_7s-h_6)/eta_P;

h_9=XSteam('h_px',P_c1,0);
h_10=h_9;
T_9=XSteam('t_ph',P_c1,h_9);

T_8=T_9;
h_8=XSteam('h_pt',P_b,T_8);

h_11=XSteam('h_pt',P_b,T_max);
s_11=XSteam('s_ph',P_b,h_11);

h_14=XSteam('h_pt',P_R,T_max);
s_14=XSteam('s_ph',P_R,h_14);

for i=1:c

    h_4(i)=XSteam('h_px',P_c2(i),0);
    h_5(i)=h_4(i);
    T_4(i)=XSteam('t_ph',P_c2(i),h_4(i));
    T_3(i)=T_4(i);

    if P_c2(i)<P_O

        h_3(i)=XSteam('h_pt',P_O,T_3(i));

        h_12s(i)=XSteam('h_ps',P_c1,s_11);
        h_12(i)=h_11-(h_11-h_12s(i))*eta_T;
        s_12(i)=XSteam('s_ph',P_c1,h_12(i));

        h_13s(i)=XSteam('h_ps',P_R,s_12(i));
        h_13(i)=h_12(i)-(h_12(i)-h_13s(i))*eta_T;

        h_15s(i)=XSteam('h_ps',P_O,s_14);
        h_15(i)=h_14-(h_14-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_O,h_15(i));

        h_16s(i)=XSteam('h_ps',P_c2(i),s_15(i));
        h_16(i)=h_15(i)-(h_15(i)-h_16s(i))*eta_T;
        s_16(i)=XSteam('s_ph',P_c2(i),h_16(i));

        h_17s(i)=XSteam('h_ps',P_c,s_16(i));
        h_17(i)=h_16(i)-(h_16(i)-h_17s(i))*eta_T;

        m_1(i)=(h_8-h_7)/(h_12(i)-h_9);
        m_2(i)=(h_6-m_1(i)*h_10+h_3(i)*(m_1(i)-1))/(h_15(i)-h_3(i));
        m_3(i)=((1-m_1(i)-m_2(i))*(h_3(i)-h_2))/(h_16(i)-h_4(i));

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i))*h_17(i)+m_3(i)*h_5(i)-(1-m_1(i)-m_2(i))*h_1;
        Q_dot_add(i)=h_11-h_8+(h_14-h_13(i))*(1-m_1(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    elseif P_c2(i)<P_R

        h_3(i)=XSteam('h_pt',P_b,T_3(i));

        h_12s(i)=XSteam('h_ps',P_c1,s_11);
        h_12(i)=h_11-(h_11-h_12s(i))*eta_T;
        s_12(i)=XSteam('s_ph',P_c1,h_12(i));

        h_13s(i)=XSteam('h_ps',P_R,s_12(i));
        h_13(i)=h_12(i)-(h_12(i)-h_13s(i))*eta_T;

        h_16s(i)=XSteam('h_ps',P_c2(i),s_14);
        h_16(i)=h_14-(h_14-h_16s(i))*eta_T;
        s_16(i)=XSteam('s_ph',P_c2(i),h_16(i));

        h_15s(i)=XSteam('h_ps',P_O,s_16(i));
        h_15(i)=h_16(i)-(h_16(i)-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_O,h_15(i));

        h_17s(i)=XSteam('h_ps',P_c,s_15(i));
        h_17(i)=h_15(i)-(h_15(i)-h_17s(i))*eta_T;

        m_1(i)=(h_8-h_3(i))/(h_12(i)-h_9);
        m_3(i)=(h_3(i)-h_7+m_1(i)*(h_4(i)-h_10))/(h_16(i)-h_4(i));
        m_2(i)=(h_6-h_2+(m_1(i)+m_3(i))*(h_2-h_5(i)))/(h_15(i)-h_2);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i))*(h_17(i)-h_1);
        Q_dot_add(i)=h_11-h_8+(h_14-h_13(i))*(1-m_1(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));    

    elseif P_c2(i)<P_c1

        h_3(i)=XSteam('h_pt',P_b,T_3(i));

        h_12s(i)=XSteam('h_ps',P_c1,s_11);
        h_12(i)=h_11-(h_11-h_12s(i))*eta_T;
        s_12(i)=XSteam('s_ph',P_c1,h_12(i));

        h_16s(i)=XSteam('h_ps',P_c2(i),s_12(i));
        h_16(i)=h_12(i)-(h_12(i)-h_16s(i))*eta_T;
        s_16(i)=XSteam('s_ph',P_c2(i),h_16(i));

        h_13s(i)=XSteam('h_ps',P_R,s_16(i));
        h_13(i)=h_16(i)-(h_16(i)-h_13s(i))*eta_T;

        h_15s(i)=XSteam('h_ps',P_O,s_14);
        h_15(i)=h_14-(h_14-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_O,h_15(i));

        h_17s(i)=XSteam('h_ps',P_c,s_15(i));
        h_17(i)=h_15(i)-(h_15(i)-h_17s(i))*eta_T;

        m_1(i)=(h_8-h_3(i))/(h_12(i)-h_9);
        m_3(i)=(h_3(i)-h_7+m_1(i)*(h_4(i)-h_10))/(h_16(i)-h_4(i));
        m_2(i)=(h_6-h_2+(m_1(i)+m_3(i))*(h_2-h_5(i)))/(h_15(i)-h_2);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i))*(h_17(i)-h_1);
        Q_dot_add(i)=h_11-h_8+(h_14-h_13(i))*(1-m_1(i)-m_3(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));           

    else

        h_3(i)=XSteam('h_pt',P_b,T_3(i));

        h_16s(i)=XSteam('h_ps',P_c2(i),s_11);
        h_16(i)=h_11-(h_11-h_16s(i))*eta_T;
        s_16(i)=XSteam('s_ph',P_c2(i),h_16(i));

        h_12s(i)=XSteam('h_ps',P_c1,s_16(i));
        h_12(i)=h_16(i)-(h_16(i)-h_12s(i))*eta_T;
        s_12(i)=XSteam('s_ph',P_c1,h_12(i));

        h_13s(i)=XSteam('h_ps',P_R,s_12(i));
        h_13(i)=h_12(i)-(h_12(i)-h_13s(i))*eta_T;

        h_15s(i)=XSteam('h_ps',P_O,s_14);
        h_15(i)=h_14-(h_14-h_15s(i))*eta_T;
        s_15(i)=XSteam('s_ph',P_O,h_15(i));

        h_17s(i)=XSteam('h_ps',P_c,s_15(i));
        h_17(i)=h_15(i)-(h_15(i)-h_17s(i))*eta_T;

        m_3(i)=(h_3(i)-h_8)/(h_16(i)-h_4(i));
        m_1(i)=(h_8-h_7+m_3(i)*(h_9-h_5(i)))/(h_12(i)-h_9);
        m_2(i)=(h_6-h_2+(m_1(i)+m_3(i))*(h_2-h_10))/(h_15(i)-h_2);

        Q_dot_rej(i)=(1-m_1(i)-m_2(i)-m_3(i))*(h_17(i)-h_1);
        Q_dot_add(i)=h_11-h_3(i)+(h_14-h_13(i))*(1-m_1(i)-m_3(i));

        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i)); 


    end
    
end


[eta_max,n]=max(eta_th);
[W_max,m]=max(W_dot_net);


P_c2_eta=P_c2(n);
P_c2_W=P_c2(m);

disp('[Efficiency]')
fprintf('\nThe optimum CFWH_2 Pressure that maximizes the efficiency is : \n(%f)\n',P_c2_eta)
fprintf('\nThe maximum efficiency is : \n(%f)\n\n\n',eta_max)

disp('[Work Net]')
fprintf('\nThe optimum CFWH_2 Pressure that maximizes the net work is : \n(%f)\n',P_c2_W)
fprintf('\nThe maximum Work Net is : \n(%f)\n',W_max)


figure(1)
% yyaxis left
plot(P_c2,eta_th) 
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-2} VS Pressure')
% yyaxis right
% plot(P_c2,W_dot_net)
grid on

