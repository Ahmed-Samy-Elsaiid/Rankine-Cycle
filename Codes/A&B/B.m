 clc;clear,clf,format long


P_c=0.1;
T_max=500;
eta_T=0.85;
eta_P=0.95;
h_1=XSteam('h_px',P_c,0);
s_1=XSteam('s_ph',P_c,h_1);


for k=1:901
    P_b(k)=119+0.1*k;
    h_2s(k)=XSteam('h_ps',P_b(k),s_1);
    h_2(k)=h_1+(h_2s(k)-h_1)/eta_P;
    h_3(k)=XSteam('h_pt',P_b(k),T_max);
    s_3(k)=XSteam('s_ph',P_b(k),h_3(k));
    for i=1:k+1199
        P_R(i)=0.1*i+0.05;
        h_4s(i)=XSteam('h_ps',P_R(i),s_3(k));
        h_4(i)=h_3(k)-(h_3(k)-h_4s(i))*eta_T;
        h_5(i)=XSteam('h_pt',P_R(i),T_max);
        s_5(i)=XSteam('s_ph',P_R(i),h_5(i));
        h_6s(i)=XSteam('h_ps',P_c,s_5(i));
        h_6(i)=h_5(i)-(h_5(i)-h_6s(i))*eta_T;
        Q_dot_add(i)=h_3(k)-h_2(k)+h_5(i)-h_4(i);
        Q_dot_rej(i)=h_6(i)-h_1;
        eta_th(i)=1-Q_dot_rej(i)/Q_dot_add(i);
        W_dot_net(i)=Q_dot_add(i)-Q_dot_rej(i);
    end
    [eta_max(k),n_2]=max(eta_th);
    [W_max(k),n_1]=max(W_dot_net);


    P_R_eta(k)=P_R(n_2);
    P_R_W(k)=P_R(n_1);


    ratio_eta(k)=P_R_eta(k)/P_b(k);
    ratio_W(k)=P_R_W(k)/P_b(k);
end

figure(1)
plot(P_b,ratio_eta)
grid on

figure(2)
plot(P_b,ratio_W)
grid on
