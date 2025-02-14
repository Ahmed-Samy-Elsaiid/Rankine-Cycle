% NOTE THAT ::
% To run the code you need to add XSteam.m to the same folder
% that contain this m.file code


clc;clear,format long

P_c=0.1;
P_b=150;
P_R=30.54;
P_O=5.93;
T_max=500;
eta_T=0.85;
eta_P=0.95;

s=0.01;P_c1=(P_c+s):s:(P_b-s);c=length(P_c1);

h_1=XSteam('h_px',P_c,0);
s_1=XSteam('s_ph',P_c,h_1);

h_2s=XSteam('h_ps',P_O,s_1);
h_2=h_1+(h_2s-h_1)/eta_P;

h_3=XSteam('h_px',P_O,0);
s_3=XSteam('s_ph',P_O,h_3);

h_4s=XSteam('h_ps',P_b,s_3);
h_4=h_3+(h_4s-h_3)/eta_P;

h_8=XSteam('h_pt',P_b,T_max);
s_8=XSteam('s_ph',P_b,h_8);

h_10=XSteam('h_pt',P_R,T_max);
s_10=XSteam('s_ph',P_R,h_10);

for i=1:c

    h_6(i)=XSteam('h_px',P_c1(i),0);
    h_7(i)=h_6(i);
    T_6(i)=XSteam('t_ph',P_c1(i),h_6(i));

    T_5(i)=T_6(i);
    

    if P_c1(i)<P_O

        h_5(i)=XSteam('h_pt',P_O,T_5(i));

        h_9s(i)=XSteam('h_ps',P_R,s_8);
        h_9(i)=h_8-(h_8-h_9s(i))*eta_T;

        h_12s(i)=XSteam('h_ps',P_O,s_10);
        h_12(i)=h_10-(h_10-h_12s(i))*eta_T;
        s_12(i)=XSteam('s_ph',P_O,h_12(i));

        h_11s(i)=XSteam('h_ps',P_c1(i),s_12(i));
        h_11(i)=h_12(i)-(h_12(i)-h_11s(i))*eta_T;
        s_11(i)=XSteam('s_ph',P_c1(i),h_11(i));

        h_13s(i)=XSteam('h_ps',P_c,s_11(i));
        h_13(i)=h_11(i)-(h_11(i)-h_13s(i))*eta_T;

        m_2(i)=(h_3-h_5(i))/(h_12(i)-h_5(i));
        m_1(i)=((1-m_2(i))*(h_5(i)-h_2))/(h_11(i)-h_6(i));

        Q_dot_add(i)=h_8-h_4+h_10-h_9(i);
        Q_dot_rej(i)=(1-m_1(i)-m_2(i))*h_13(i)+m_1(i)*h_7(i)-h_1*(1-m_2(i));
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));

    elseif P_c1(i)<P_R

        h_5(i)=XSteam('h_pt',P_b,T_5(i));

        h_9s(i)=XSteam('h_ps',P_R,s_8);
        h_9(i)=h_8-(h_8-h_9s(i))*eta_T;

        h_11s(i)=XSteam('h_ps',P_c1(i),s_10);
        h_11(i)=h_10-(h_10-h_11s(i))*eta_T;
        s_11(i)=XSteam('s_ph',P_c1(i),h_11(i));

        h_12s(i)=XSteam('h_ps',P_O,s_11(i));
        h_12(i)=h_11(i)-(h_11(i)-h_12s(i))*eta_T;
        s_12(i)=XSteam('s_ph',P_O,h_12(i));

        h_13s(i)=XSteam('h_ps',P_c,s_12(i));
        h_13(i)=h_12(i)-(h_12(i)-h_13s(i))*eta_T;

        m_1(i)=(h_5(i)-h_4)/(h_11(i)-h_6(i));
        m_2(i)=(h_3-m_1(i)*h_7(i)+h_2*(m_1(i)-1))/(h_12(i)-h_2);

        Q_dot_add(i)=h_8-h_5(i)+h_10-h_9(i);
        Q_dot_rej(i)=(1-m_1(i)-m_2(i))*(h_13(i)-h_1);
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
        
    else 

        h_5(i)=XSteam('h_pt',P_b,T_5(i));

        h_11s(i)=XSteam('h_ps',P_c1(i),s_8);
        h_11(i)=h_8-(h_8-h_11s(i))*eta_T;
        s_11(i)=XSteam('s_ph',P_c1(i),h_11(i));

        h_9s(i)=XSteam('h_ps',P_R,s_11(i));
        h_9(i)=h_11(i)-(h_11(i)-h_9s(i))*eta_T;

        h_12s(i)=XSteam('h_ps',P_O,s_10);
        h_12(i)=h_10-(h_10-h_12s(i))*eta_T;
        s_12(i)=XSteam('s_ph',P_O,h_12(i));

        h_13s(i)=XSteam('h_ps',P_c,s_12(i));
        h_13(i)=h_12(i)-(h_12(i)-h_13s(i))*eta_T;

        m_1(i)=(h_5(i)-h_4)/(h_11(i)-h_6(i));
        m_2(i)=(h_3-m_1(i)*h_7(i)+h_2*(m_1(i)-1))/(h_12(i)-h_2);

        Q_dot_add(i)=h_8-h_5(i)+(h_10-h_9(i))*(1-m_1(i));
        Q_dot_rej(i)=(1-m_1(i)-m_2(i))*(h_13(i)-h_1);
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
    end
end


[eta_max,n_2]=max(eta_th);
[W_max,n_1]=max(W_dot_net);


P_c1_eta=P_c1(n_2);
P_c1_W=P_c1(n_1);

disp('[Efficiency]')
fprintf('\nThe optimum CFWH_1 Pressure that maximizes the efficiency is : \n(%f)\n',P_c1_eta)
fprintf('\nThe maximum efficiency is : \n(%f)\n\n\n',eta_max)

disp('[Work Net]')
fprintf('\nThe optimum CFWH_1 Pressure that maximizes the net work is : \n(%f)\n',P_c1_W)
fprintf('\nThe maximum Work Net is : \n(%f)\n',W_max)


figure(1)
% yyaxis left
plot(P_c1,eta_th)
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-1} VS Pressure')
% yyaxis right
% plot(P_c1,W_dot_net)
grid on