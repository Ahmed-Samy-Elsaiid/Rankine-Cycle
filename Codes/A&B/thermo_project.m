clc;clear,format long


P_boiler=150;
T_max=500;
P_conenser=0.1;
P_reheat=30.54;

eta_T=0.85;
eta_P=0.95;

x=P_conenser:0.1:P_boiler-0.1;

%case1 1-OFWH ===========================================================

h1=XSteam('hL_p',P_conenser);
s1=XSteam('s_ph',P_conenser,h1);
h2=zeros(size(x));%pump
h3=zeros(size(x));
s3=zeros(size(x));
h4=zeros(size(x));%pump
h5=XSteam('h_pT',P_boiler,T_max); s5=XSteam('s_pT',P_boiler,T_max);
h6=zeros(size(x));%pump
h7=XSteam('h_pT',P_reheat,500); s7=XSteam('s_pT',P_reheat,500);
h8=zeros(size(x));%turbine
s8=zeros(size(x));
h9=zeros(size(x));%turbine


m1=zeros(size(x));

q_in=zeros(size(x));
q_out=zeros(size(x));

W_net=zeros(size(x));
eta_th=zeros(size(x));

for i=1:length(x)
    h2(i)=XSteam('h_ps',x(i),s1);
    h2(i)=((h2(i)-h1)/eta_P) + h1;
    h3(i)=XSteam('hL_p',x(i));
    s3(i)=XSteam('s_ph',x(i),h3(i));
    h4(i)=XSteam('h_ps',P_boiler,s3(i));
    h4(i)=((h4(i)-h3(i))/eta_P) + h3(i);
    
    if x(i)< P_reheat
        h6(i)=XSteam('h_ps',P_reheat,s5);%turbine
        h6(i)=h5-eta_T*(h5-h6(i));
        h8(i)=XSteam('h_ps',x(i),s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',x(i),h8(i));
        h9(i)=XSteam('h_ps',P_conenser,s8(i));%turbine
        h9(i)=h8(i)-eta_T*(h8(i)-h9(i));
        
        m1(i)=(h3(i)-h2(i))/(h8(i)-h2(i));
        
        q_out(i)=(h9(i)-h1)*(1-m1(i));
        q_in(i)=(h5-h4(i)+h7-h6(i));  
    else
        h8(i)=XSteam('h_ps',x(i),s5);
        h8(i)=h5-eta_T*(h5-h8(i));
        s8(i)=XSteam('s_ph',x(i),h8(i));
        h6(i)=XSteam('h_ps',P_reheat,s8(i));%turbine
        h6(i)=h8(i)-eta_T*(h8(i)-h6(i));
        h9(i)=XSteam('h_ps',P_conenser,s7);%turbine
        h9(i)=h7-eta_T*(h7-h9(i));
        
        m1(i)=(h3(i)-h2(i))/(h8(i)-h2(i));
        
        q_out(i)=(h9(i)-h1)*(1-m1(i));
        q_in(i)=(h5-h4(i))+(h7-h6(i))*(1-m1(i));
        
    end
    W_net(i)=q_in(i)-q_out(i);
    eta_th(i)=1-(q_out(i)/q_in(i));
end

[eta_th_max,i]=max(eta_th);
P_OFWH_opt=x(i);

figure(1)
plot(x,eta_th,x,(W_net/2500),':r')
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,OFWH} VS Pressure')
leg1=legend('\eta_{th}','W_{net}');
set(leg1,'Location','best')

h2=h2(i);
h3=h3(i);
s3=s3(i);
h4=h4(i);
h6=h6(i);
h8=h8(i);
s8=s8(i);
h9=h9(i);

%case2  1-CFWH & 1-OFWH =============================================

h10=zeros(size(x));
s10=zeros(size(x));
h11=zeros(size(x));
h12=zeros(size(x));
m2=zeros(size(x));
W_net_2=zeros(size(x));
eta_th_2=zeros(size(x));

for i=1:length(x)
    h11(i)=XSteam('hL_p',x(i));
    
    if x(i) < P_OFWH_opt
        h10(i)=XSteam('h_ps',x(i),s8);
        h10(i)=h8-eta_T*(h8-h10(i));
        s10(i)=XSteam('s_ph',x(i),h10(i));
        h9(i)=XSteam('h_ps',P_conenser,s10(i));%turbine
        h9(i)=h10(i)-eta_T*(h10(i)-h9(i));
        h12(i)=XSteam('h_pT',P_OFWH_opt,XSteam('Tsat_p',x(i)));
        
        m1(i)=(h3-h12(i))/(h8-h12(i));
        m2(i)=((1-m1(i))*(h12(i)-h2))/(h10(i)-h11(i));
        
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i))+(h11(i)-h1)*m2(i);
        q_in(i)=(h5-h4+h7-h6);
    elseif x(i) < P_reheat
        h10(i)=XSteam('h_ps',x(i),s7);
        h10(i)=h7-eta_T*(h7-h10(i));
        s10(i)=XSteam('s_ph',x(i),h10(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s10(i));
        h8(i)=h10(i)-eta_T*(h10(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h9(i)=XSteam('h_ps',P_conenser,s8(i));%turbine
        h9(i)=h8(i)-eta_T*(h8(i)-h9(i));
        h12(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
         
        m2(i)=(h12(i)-h4)/(h10(i)-h11(i));
        m1(i)=(m2(i)*(h2-h11(i))+h3-h2)/(h8(i)-h2);
        
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i));
        q_in(i)=(h5-h12(i)+h7-h6);
    else
        h10(i)=XSteam('h_ps',x(i),s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',x(i),h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h12(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h9(i)=XSteam('h_ps',P_conenser,s8(i));%turbine
        h9(i)=h8(i)-eta_T*(h8(i)-h9(i));
        
        m2(i)=(h12(i)-h4)/(h10(i)-h11(i));
        m1(i)=(m2(i)*(h2-h11(i))+h3-h2)/(h8(i)-h2);
         
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i));
        q_in(i)=(h5-h12(i)+(h7-h6(i))*(1-m2(i)));
    end
    W_net_2(i)=q_in(i)-q_out(i);
    eta_th_2(i)=1-q_out(i)/q_in(i);
end

[eta_th_2_max,i]=max(eta_th_2);
P_CFWH_1_opt=x(i);

figure(2)
plot(x,eta_th_2,x,(W_net_2/2500),':r')
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-1} VS Pressure')
leg2=legend('\eta_{th}','W_{net}');
set(leg2,'Location','best')

h6=h6(i);
h8=h8(i);
s8=s8(i);
h9=h9(i);
h10=h10(i);
s10=s10(i);
h11=h11(i);
h12=h12(i);

%case3  2-CFWH & 1-OFWH =============================================

h13=zeros(size(x));
s13=zeros(size(x));
h14=zeros(size(x));
h15=zeros(size(x));
m3=zeros(size(x));
W_net_3=zeros(size(x));
eta_th_3=zeros(size(x));

for i=1:length(x)
    h14(i)=XSteam('hL_p',x(i));
    
    if x(i) < P_OFWH_opt
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',x(i),s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',x(i),h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));%turbine
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h15(i)=XSteam('h_pT',P_OFWH_opt,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h4)/(h10(i)-h11);
        m1(i)=(m2(i)*(h15(i)-h11)+h3-h15(i))/(h8(i)-h15(i));
        m3(i)=((1-m1(i)-m2(i))*(h15(i)-h2))/(h13(i)-h14(i));
         
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i))+(h14(i)-h1)*m3(i);
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_reheat
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h13(i)=XSteam('h_ps',x(i),s7);
        h13(i)=h7-eta_T*(h7-h13(i));
        s13(i)=XSteam('s_ph',x(i),h13(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s13(i));
        h8(i)=h13(i)-eta_T*(h13(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h9(i)=XSteam('h_ps',P_conenser,s8(i));%turbine
        h9(i)=h8(i)-eta_T*(h8(i)-h9(i));
        h15(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h15(i))/(h10(i)-h11);
        m3(i)=(m2(i)*(h14(i)-h11)+h15(i)-h4)/(h13(i)-h14(i));
        m1(i)=((m2(i)+m3(i))*(h2-h14(i))+h3-h2)/(h8(i)-h2);
         
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i));
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_CFWH_1_opt
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h13(i)=XSteam('h_ps',x(i),s10(i));
        h13(i)=h10(i)-eta_T*(h10(i)-h13(i));
        s13(i)=XSteam('s_ph',x(i),h13(i));
        h6(i)=XSteam('h_ps',P_reheat,s13(i));%turbine
        h6(i)=h13(i)-eta_T*(h13(i)-h6(i));
        h15(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h9(i)=XSteam('h_ps',P_conenser,s8(i));%turbine
        h9(i)=h8(i)-eta_T*(h8(i)-h9(i));
        
        m2(i)=(h12-h15(i))/(h10(i)-h11);
        m3(i)=(m2(i)*(h14(i)-h11)+h15(i)-h4)/(h13(i)-h14(i));
        m1(i)=((m2(i)+m3(i))*(h2-h14(i))+h3-h2)/(h8(i)-h2);
         
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i));
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i)-m3(i));
    else
        h13(i)=XSteam('h_ps',x(i),s5);
        h13(i)=h5-eta_T*(h5-h13(i));
        s13(i)=XSteam('s_ph',x(i),h13(i));
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s13(i));
        h10(i)=h13(i)-eta_T*(h13(i)-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h15(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h9(i)=XSteam('h_ps',P_conenser,s8(i));%turbine
        h9(i)=h8(i)-eta_T*(h8(i)-h9(i));
        
        m3(i)=(h15(i)-h12)/(h13(i)-h14(i));
        m2(i)=(m3(i)*(h11-h14(i))+h12-h4)/(h10(i)-h11);
        m1(i)=((m2(i)+m3(i))*(h2-h11)+h3-h2)/(h8(i)-h2);
         
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i));
        q_in(i)=h5-h15(i)+(h7-h6(i))*(1-m2(i)-m3(i));
    end
    W_net_3(i)=q_in(i)-q_out(i);
    eta_th_3(i)=1-q_out(i)/q_in(i);
end

[eta_th_3_max,i]=max(eta_th_3);
P_CFWH_2_opt=x(i);

figure(3)
plot(x,eta_th_3,x,(W_net_3/2500),':r')
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-2} VS Pressure')
leg3=legend('\eta_{th}','W_{net}');
set(leg3,'Location','best')

h6=h6(i);
h8=h8(i);
s8=s8(i);
h9=h9(i);
h10=h10(i);
s10=s10(i);
h13=h13(i);
s13=s13(i);
h14=h14(i);
h15=h15(i);

%case4  3-CFWH & 1-OFWH =============================================

h16=zeros(size(x));
s16=zeros(size(x));
h17=zeros(size(x));
h18=zeros(size(x));
m4=zeros(size(x));
W_net_4=zeros(size(x));
eta_th_4=zeros(size(x));

for i=1:length(x)
    h17(i)=XSteam('hL_p',x(i));
    
    if x(i) < P_CFWH_2_opt
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h16(i)=XSteam('h_ps',x(i),s13(i));
        h16(i)=h13(i)-eta_T*(h13(i)-h16(i));
        s16(i)=XSteam('s_ph',x(i),h16(i));
        h9(i)=XSteam('h_ps',P_conenser,s16(i));%turbine
        h9(i)=h16(i)-eta_T*(h16(i)-h9(i));
        h18(i)=XSteam('h_pT',P_OFWH_opt,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h4)/(h10(i)-h11);
        m1(i)=(m2(i)*(h15-h11)+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i))*(h15-h18(i)))/(h13(i)-h14);
        m4(i)=((1-m1(i)-m2(i))*(h18(i)-h2)+m3(i)*(h17(i)-h14))/(h16(i)-h17(i));
         
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i))+(h17(i)-h1)*(m3(i)+m4(i));
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_OFWH_opt
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h16(i)=XSteam('h_ps',x(i),s8(i));
        h16(i)=h8(i)-eta_T*(h8(i)-h16(i));
        s16(i)=XSteam('s_ph',x(i),h16(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s16(i));
        h13(i)=h16(i)-eta_T*(h16(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));%turbine
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h18(i)=XSteam('h_pT',P_OFWH_opt,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h4)/(h10(i)-h11);
        m1(i)=(m2(i)*(h18(i)-h11)+h3-h18(i))/(h8(i)-h18(i));
        m4(i)=((1-m1(i)-m2(i))*(h18(i)-h15))/(h16(i)-h17(i));
        m3(i)=((1-m1(i)-m2(i))*(h15-h2)+m4(i)*(h14-h17(i)))/(h13(i)-h14);
          
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i))+(h14-h1)*(m3(i)+m4(i));
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_reheat
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h16(i)=XSteam('h_ps',x(i),s7);
        h16(i)=h7-eta_T*(h7-h16(i));
        s16(i)=XSteam('s_ph',x(i),h16(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s16(i));
        h8(i)=h16(i)-eta_T*(h16(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h18(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h18(i))/(h10(i)-h11);
        m4(i)=((h18(i)-h4)+m2(i)*(h17(i)-h11))/(h16(i)-h17(i));
        m1(i)=((m2(i)+m4(i))*(h15-h17(i))+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i)-m4(i))*(h15-h2))/(h13(i)-h14); 
          
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i))+(h14-h1)*m3(i);
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    else
        h8(i)=XSteam('h_ps',P_OFWH_opt,s7);
        h8(i)=h7-eta_T*(h7-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));%turbine
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h16(i)=XSteam('h_ps',x(i),s5);
        h16(i)=h5-eta_T*(h5-h16(i));
        s16(i)=XSteam('s_ph',x(i),h16(i));
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s16(i));
        h10(i)=h16(i)-eta_T*(h16(i)-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));%turbine
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h18(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
        
        m4(i)=(h18(i)-h12)/(h16(i)-h17(i));
        m2(i)=(m4(i)*(h11-h17(i))+h12-h4)/(h10(i)-h11);
        m1(i)=((m2(i)+m4(i))*(h15-h11)+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i)-m4(i))*(h15-h2))/(h13(i)-h14);
         
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i))+(h14-h1)*m3(i);
        q_in(i)=(h5-h18(i))+(h7-h6(i))*(1-m2(i)-m4(i));
    end
    W_net_4(i)=q_in(i)-q_out(i);
    eta_th_4(i)=1-q_out(i)/q_in(i);
end

[eta_th_4_max,i]=max(eta_th_4);
P_CFWH_3_opt=x(i);

figure(4)
plot(x,eta_th_4,x,(W_net_4/2500),':r')
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-3} VS Pressure')
leg4=legend('\eta_{th}','W_{net}');
set(leg4,'Location','best')

h6=h6(i);
h8=h8(i);
s8=s8(i);
h9=h9(i);
h10=h10(i);
s10=s10(i);
h13=h13(i);
s13=s13(i);
h16=h16(i);
s16=s16(i);
h17=h17(i);
h18=h18(i);

%case5  4-CFWH & 1-OFWH =============================================

h19=zeros(size(x));
s19=zeros(size(x));
h20=zeros(size(x));
h21=zeros(size(x));
m5=zeros(size(x));
W_net_5=zeros(size(x));
eta_th_5=zeros(size(x));

for i=1:length(x)
    h20(i)=XSteam('hL_p',x(i));
        
    if x(i) < P_CFWH_2_opt
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h16(i)=XSteam('h_ps',P_CFWH_3_opt,s7);
        h16(i)=h7-eta_T*(h7-h16(i));
        s16(i)=XSteam('s_ph',P_CFWH_3_opt,h16(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s16(i));
        h8(i)=h16(i)-eta_T*(h16(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h19(i)=XSteam('h_ps',x(i),s13(i));
        h19(i)=h13(i)-eta_T*(h13(i)-h19(i));
        s19(i)=XSteam('s_ph',x(i),h19(i));
        h9(i)=XSteam('h_ps',P_conenser,s19(i));
        h9(i)=h19(i)-eta_T*(h19(i)-h9(i));
        h21(i)=XSteam('h_pT',P_OFWH_opt,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h18)/(h10(i)-h11);
        m4(i)=((h18-h4)+m2(i)*(h17-h11))/(h16(i)-h17);
        m1(i)=((m2(i)+m4(i))*(h15-h17)+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i)-m4(i))*(h15-h21(i)))/(h13(i)-h14);
        m5(i)=((1-m1(i)-m2(i)-m4(i))*(h21(i)-h2)+m3(i)*(h20(i)-h14))/(h19(i)-h20(i));
        
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i)-m5(i))+(h20(i)-h1)*(m3(i)+m5(i));
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_OFWH_opt
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h16(i)=XSteam('h_ps',P_CFWH_3_opt,s7);
        h16(i)=h7-eta_T*(h7-h16(i));
        s16(i)=XSteam('s_ph',P_CFWH_3_opt,h16(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s16(i));
        h8(i)=h16(i)-eta_T*(h16(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h19(i)=XSteam('h_ps',x(i),s8(i));
        h19(i)=h8(i)-eta_T*(h8(i)-h19(i));
        s19(i)=XSteam('s_ph',x(i),h19(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s19(i));
        h13(i)=h19(i)-eta_T*(h19(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h21(i)=XSteam('h_pT',P_OFWH_opt,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h18)/(h10(i)-h11);
        m4(i)=((h18-h4)+m2(i)*(h17-h11))/(h16(i)-h17);
        m1(i)=((m2(i)+m4(i))*(h21(i)-h17)+h3-h21(i))/(h8(i)-h21(i));
        m5(i)=((1-m1(i)-m2(i)-m4(i))*(h21(i)-h15))/(h19(i)-h20(i));
        m3(i)=((1-m1(i)-m2(i)-m4(i))*(h15-h2)+m5(i)*(h14-h20(i)))/(h13(i)-h14);
            
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i)-m5(i))+(h14-h1)*(m3(i)+m5(i));
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_CFWH_3_opt
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h16(i)=XSteam('h_ps',P_CFWH_3_opt,s7);
        h16(i)=h7-eta_T*(h7-h16(i));
        s16(i)=XSteam('s_ph',P_CFWH_3_opt,h16(i));
        h19(i)=XSteam('h_ps',x(i),s16(i));
        h19(i)=h16(i)-eta_T*(h16(i)-h19(i));
        s19(i)=XSteam('s_ph',x(i),h19(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s19(i));
        h8(i)=h19(i)-eta_T*(h19(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h21(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
        
        m2(i)=(h12-h18)/(h10(i)-h11);
        m4(i)=((h18-h21(i))+m2(i)*(h17-h11))/(h16(i)-h17);
        m5(i)=((h21(i)-h4)+(m2(i)+m4(i))*(h20(i)-h17))/(h19(i)-h20(i));
        m1(i)=((m2(i)+m4(i)+m5(i))*(h15-h20(i))+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i)-m4(i)-m5(i))*(h15-h2))/(h13(i)-h14);
        
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i)-m5(i))+(h14-h1)*m3(i);
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_reheat
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h19(i)=XSteam('h_ps',x(i),s7);
        h19(i)=h7-eta_T*(h7-h19(i));
        s19(i)=XSteam('s_ph',x(i),h19(i));
        h16(i)=XSteam('h_ps',P_CFWH_3_opt,s19(i));
        h16(i)=h19(i)-eta_T*(h19(i)-h16(i));
        s16(i)=XSteam('s_ph',P_CFWH_3_opt,h16(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s16(i));
        h8(i)=h16(i)-eta_T*(h16(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h21(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
          
        m2(i)=(h12-h21(i))/(h10(i)-h11);
        m5(i)=((h21(i)-h18)+m2(i)*(h20(i)-h11))/(h19(i)-h20(i));
        m4(i)=((h18-h4)+(m2(i)+m5(i))*(h17-h20(i)))/(h16(i)-h17);
        m1(i)=((m2(i)+m4(i)+m5(i))*(h15-h20(i))+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i)-m4(i)-m5(i))*(h15-h2))/(h13(i)-h14);
           
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i)-m5(i))+(h14-h1)*m3(i);
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i));
    elseif x(i) < P_CFWH_1_opt
        h16(i)=XSteam('h_ps',P_CFWH_3_opt,s7);
        h16(i)=h7-eta_T*(h7-h16(i));
        s16(i)=XSteam('s_ph',P_CFWH_3_opt,h16(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s16(i));
        h8(i)=h16(i)-eta_T*(h16(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s5);
        h10(i)=h5-eta_T*(h5-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h19(i)=XSteam('h_ps',x(i),s10(i));
        h19(i)=h10(i)-eta_T*(h10(i)-h19(i));
        s19(i)=XSteam('s_ph',x(i),h19(i));
        h6(i)=XSteam('h_ps',P_reheat,s19(i));
        h6(i)=h19(i)-eta_T*(h19(i)-h6(i));
        h21(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
          
        m2(i)=(h12-h21(i))/(h10(i)-h11);
        m5(i)=((h21(i)-h18)+m2(i)*(h20(i)-h11))/(h19(i)-h20(i));
        m4(i)=((h18-h4)+(m2(i)+m5(i))*(h17-h20(i)))/(h16(i)-h17);
        m1(i)=((m2(i)+m4(i)+m5(i))*(h15-h17)+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i)-m4(i)-m5(i))*(h15-h2))/(h13(i)-h14);
            
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i)-m5(i))+(h14-h1)*m3(i);
        q_in(i)=h5-h12+(h7-h6(i))*(1-m2(i)-m5(i));
    elseif x(i) < P_boiler
        h16(i)=XSteam('h_ps',P_CFWH_3_opt,s7);
        h16(i)=h7-eta_T*(h7-h16(i));
        s16(i)=XSteam('s_ph',P_CFWH_3_opt,h16(i));
        h8(i)=XSteam('h_ps',P_OFWH_opt,s16(i));
        h8(i)=h16(i)-eta_T*(h16(i)-h8(i));
        s8(i)=XSteam('s_ph',P_OFWH_opt,h8(i));
        h13(i)=XSteam('h_ps',P_CFWH_2_opt,s8(i));
        h13(i)=h8(i)-eta_T*(h8(i)-h13(i));
        s13(i)=XSteam('s_ph',P_CFWH_2_opt,h13(i));
        h9(i)=XSteam('h_ps',P_conenser,s13(i));
        h9(i)=h13(i)-eta_T*(h13(i)-h9(i));
        h19(i)=XSteam('h_ps',x(i),s5);
        h19(i)=h5-eta_T*(h5-h19(i));
        s19(i)=XSteam('s_ph',x(i),h19(i));
        h10(i)=XSteam('h_ps',P_CFWH_1_opt,s19(i));
        h10(i)=h19(i)-eta_T*(h19(i)-h10(i));
        s10(i)=XSteam('s_ph',P_CFWH_1_opt,h10(i));
        h6(i)=XSteam('h_ps',P_reheat,s10(i));
        h6(i)=h10(i)-eta_T*(h10(i)-h6(i));
        h21(i)=XSteam('h_pT',P_boiler,XSteam('Tsat_p',x(i)));
         
        m5(i)=(h21(i)-h12)/(h19(i)-h20(i));
        m2(i)=((h12-h18)+m5(i)*(h11-h20(i)))/(h10(i)-h11);
        m4(i)=((h18-h4)+(m2(i)+m5(i))*(h17-h11))/(h16(i)-h17);
        m1(i)=((m2(i)+m4(i)+m5(i))*(h15-h17)+h3-h15)/(h8(i)-h15);
        m3(i)=((1-m1(i)-m2(i)-m4(i)-m5(i))*(h15-h2))/(h13(i)-h14);
          
        q_out(i)=(h9(i)-h1)*(1-m1(i)-m2(i)-m3(i)-m4(i)-m5(i))+(h14-h1)*m3(i);
        q_in(i)=(h5-h21(i))+(h7-h6(i))*(1-m2(i)-m5(i));   
    end
    W_net_5(i)=q_in(i)-q_out(i);
    eta_th_5(i)=1-q_out(i)/q_in(i);
end

[eta_th_5_max,i]=max(eta_th_5);
P_CFWH_4_opt=x(i);

figure(5)
plot(x,eta_th_5,x,(W_net_5/2500),':r')
xlabel('Pressure')
ylabel('\eta_{th}')
title('\eta_{th,CFWH-4} VS Pressure')
leg5=legend('\eta_{th}','W_{net}');
set(leg5,'Location','best')

h6=h6(i);
h8=h8(i);
s8=s8(i);
h9=h9(i);
h10=h10(i);
s10=s10(i);
h13=h13(i);
s13=s13(i);
h16=h16(i);
s16=s16(i);
h19=h19(i);
s19=s19(i);  
h20=h21(i);
h21=h21(i);