% NOTE THAT ::
% To run the code you need to add XSteam.m to the same folder
% that contain this m.file code


clc;clear,format long


P_c=0.1;
P_b=150;
P_R=30.54;
T_max=500;
eta_T=0.85;
eta_P=0.95;


s=0.01;P_O=(P_c+s):s:(P_b-s);c=length(P_O);


h_1=XSteam('h_px',P_c,0);
s_1=XSteam('s_ph',P_c,h_1);
h_5=XSteam('h_pt',P_b,T_max);
s_5=XSteam('s_ph',P_b,h_5);
h_8=XSteam('h_pt',P_R,T_max);
s_8=XSteam('s_ph',P_R,h_8);


for i=1:c

    h_2s(i)=XSteam('h_ps',P_O(i),s_1);
    h_2(i)=h_1+(h_2s(i)-h_1)/eta_P;
    h_3(i)=XSteam('h_px',P_O(i),0);
    s_3(i)=XSteam('s_ph',P_O(i),h_3(i));
    h_4s(i)=XSteam('h_ps',P_b,s_3(i));
    h_4(i)=h_3(i)+(h_4s(i)-h_3(i))/eta_P;
    if P_O(i)<P_R
        h_7s(i)=XSteam('h_ps',P_R,s_5);
        h_7(i)=h_5-(h_5-h_7s(i))*eta_T;
        h_6s(i)=XSteam('h_ps',P_O(i),s_8);
        h_6(i)=h_8-(h_8-h_6s(i))*eta_T;
        s_6(i)=XSteam('s_ph',P_O(i),h_6(i));
        h_9s(i)=XSteam('h_ps',P_c,s_6(i));
        h_9(i)=h_6(i)-(h_6(i)-h_9s(i))*eta_T;
        m_1(i)=(h_3(i)-h_2(i))/(h_6(i)-h_2(i));
        Q_dot_add(i)=h_5-h_4(i)+h_8-h_7(i);
        Q_dot_rej(i)=(1-m_1(i))*(h_9(i)-h_1);
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
    else 
        h_6s(i)=XSteam('h_ps',P_O(i),s_5);
        h_6(i)=h_5-(h_5-h_6s(i))*eta_T;
        s_6(i)=XSteam('s_ph',P_O(i),h_6(i));
        h_7s(i)=XSteam('h_ps',P_R,s_6(i));
        h_7(i)=h_6(i)-(h_6(i)-h_7s(i))*eta_T;
        h_9s(i)=XSteam('h_ps',P_c,s_8);
        h_9(i)=h_8-(h_8-h_9s(i))*eta_T;
        m_1(i)=(h_3(i)-h_2(i))/(h_6(i)-h_2(i));
        Q_dot_add(i)=h_5-h_4(i)+(1-m_1(i))*(h_8-h_7(i));
        Q_dot_rej(i)=(1-m_1(i))*(h_9(i)-h_1);
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
        eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
    end
end


[eta_max,n_2]=max(eta_th);
[W_max,n_1]=max(W_dot_net);


P_O_eta=P_O(n_2);
P_O_W=P_O(n_1);


disp('[Efficiency]')
fprintf('\nThe optimum OFWH Pressure that maximizes the efficiency is : \n(%f)\n',P_O_eta)
fprintf('\nThe maximum efficiency is : \n(%f)\n\n\n',eta_max)

disp('[Work Net]')
fprintf('\nThe optimum OFWH Pressure that maximizes the net work is : \n(%f)\n',P_O_W)
fprintf('\nThe maximum Work Net is : \n(%f)\n',W_max)


figure(1)
yyaxis left
plot(P_O,eta_th)
xlabel('Pressure')
ylabel('\eta_{th}')

yyaxis right
plot(P_O,W_dot_net)
xlabel('Pressure')
ylabel('W_{net}')

leg=legend('\eta_{th}','W_{net}');
set(leg,'Location','best')
grid on

