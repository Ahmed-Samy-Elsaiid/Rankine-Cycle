clc; clf;
p_c=0.1;
T_max=500;
efficiency_pump=.95;
efficiency_turbine=.85;
x_1=0;

h_1=XSteam('h_px',p_c,x_1);
s_1=XSteam('s_ph',p_c,h_1);
T_1=XSteam('T_ps',p_c,s_1);


for k=1:901
    p_b(k)=119+0.1*k;

    h_2s(k)=h_1+(0.1.*(p_b(k)-p_c));
    h_2(k)=((h_2s(k)-h_1)./efficiency_pump)+h_1;
    s_2s(k)=XSteam('s_ph',p_b(k),h_2s(k));
    s_2(k)=XSteam('s_ph',p_b(k),h_2(k));
    T_2s(k)=XSteam('T_ps',p_b(k),s_2s(k));
    T_2(k)=XSteam('T_ph',p_b(k),h_2(k));


    h_3(k)=XSteam('h_pT',p_b(k),T_max);
    s_3(k)=XSteam('s_ph',p_b(k),h_3(k));
    T_3(k)=XSteam('T_ph',p_b(k),h_3(k));

    for i=1:k+1199
        p_r(i)=0.1*i+0.05;

        s_4s(i)=s_3(k);
        h_4s(i)=XSteam('h_ps',p_r(i),s_4s(i));
        h_4(i)=h_3(k)-(efficiency_turbine.*(h_3(k)-h_4s(i)));
        s_4(i)=XSteam('s_ph',p_r(i),h_4(i));
        T_4(i)=XSteam('T_hs',h_4(i),s_4(i));
        T_4(i)=XSteam('T_hs',h_4s(i),s_4s(i));

        h_5(i)=XSteam('h_pT',p_r(i),T_max);
        s_5(i)=XSteam('s_pT',p_r(i),T_max);
        T_5(i)=XSteam('T_hs',h_5(i),s_5(i));



        s_6s(i)=s_5(i);
        h_6s(i)=XSteam('h_ps',p_c,s_6s(i));
        h_6(i)=h_5(i)-(efficiency_turbine*(h_5(i)-h_6s(i)));
        s_6(i)=XSteam('s_ph',p_c,h_6(i));
        T_6(i)=XSteam('T_hs',h_6(i),s_6(i));
        T_6s(i)=XSteam('T_hs',h_6s(i),s_6s(i));



        W_net(i)=((h_3(k)-h_4(i))+(h_5(i)-h_6(i)))-(h_2(k)-h_1);

        Q_in(i)=(h_3(k)-h_2(k))+(h_5(i)-h_4(i));

        Q_out(i)=(h_6(i)-h_1);

        efficiency(i)= 1-Q_out(i)/Q_in(i); 
    end
    [efficiency_opt(k),n]=max(efficiency);
    [W_net_opt(k),m]=max(W_net);

    p_r_opt1(k)=p_r(n);
    p_r_opt2(k)=p_r(m);

    ratio_eta(k)=p_r_opt1(k)/p_b(k);
    ratio_W(k)=p_r_opt2(k)/p_b(k);
end

figure(1)
plot(p_b,ratio_eta)
grid on
% hold on
figure(2)
plot(p_b,ratio_W)
grid on 


