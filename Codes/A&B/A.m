% NOTE THAT ::
% To run the code you need to add XSteam.m to the same folder
% that contain this m.file code


clc;clear,format long


P_c=0.1;
P_b=150;
T_max=500;
eta_T=0.85;
eta_P=0.95;


s=0.01;P_R=(P_c+s):s:(P_b-s);c=length(P_R);


h_1=XSteam('h_px',P_c,0);
s_1=XSteam('s_ph',P_c,h_1);
h_2s=XSteam('h_ps',P_b,s_1);
h_2=h_1+(h_2s-h_1)/eta_P;
h_3=XSteam('h_pt',P_b,T_max);
s_3=XSteam('s_ph',P_b,h_3);


for i=1:c

    h_4s(i)=XSteam('h_ps',P_R(i),s_3);
    h_4(i)=h_3-(h_3-h_4s(i))*eta_T;
    h_5(i)=XSteam('h_pt',P_R(i),T_max);
    s_5(i)=XSteam('s_ph',P_R(i),h_5(i));
    h_6s(i)=XSteam('h_ps',P_c,s_5(i));
    h_6(i)=h_5(i)-(h_5(i)-h_6s(i))*eta_T;
    Q_dot_add(i)=h_3-h_2+h_5(i)-h_4(i);
    Q_dot_rej(i)=h_6(i)-h_1;
    eta_th(i)=1-(Q_dot_rej(i)/Q_dot_add(i));
    W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
end


[W_max,n_1]=max(W_dot_net);
[eta_max,n_2]=max(eta_th);


P_R_eta=P_R(n_2);
P_R_W=P_R(n_1);


disp('[Efficiency]')
fprintf('\nThe optimum Reheat Pressure that maximizes the efficiency is : \n(%f)\n',P_R_eta)
fprintf('\nThe maximum efficiency is : \n(%f)\n\n\n',eta_max)

disp('[Work Net]')
fprintf('\nThe optimum Reheat Pressure that maximizes the net work is : \n(%f)\n',P_R_W)
fprintf('\nThe maximum Work Net is : \n(%f)\n',W_max)


figure(1)
plot(P_R,eta_th)
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th} VS Pressure')
grid on


figure(2)
plot(P_R,W_dot_net)
xlabel('Pressure')
ylabel('W_{net}')
title('W_{net} VS Pressure')
grid on

