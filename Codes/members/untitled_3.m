clc;clear

%givens
Pc=0.1;%bar
Pb=150;%bar
Pr=30.54;%bar
Tmax=500;%c

%Entropy
s5=XSteam('s_pT',Pb,Tmax);
s1=XSteam('sL_p',Pc);

%Enthalpy
h1=XSteam('hL_p',Pc);
h5=XSteam('h_pT',Pb,Tmax);

step=0.05;
c=(Pb-Pc-step)/step;

for i=1:c
    Po(i)=Pc+step*i;

if Po(i)<Pr

    h2s=XSteam('h_ps',Po(i),s1);
    h2=((h2s-h1)/0.95)+h1;

    h3=XSteam('hL_p',Po(i));
    s3=XSteam('sL_p',Po(i));
    h4s=XSteam('h_ps',Pb,s3);
    h4=((h4s-h3)/0.95)+h3;
    
    h6s=XSteam('h_ps',Pr,s5);
    h6=-(((h5-h6s)*0.85)-h5);

    h7=XSteam('h_pT',Pr,Tmax);
    s7=XSteam('s_pT',Pr,Tmax);

    h8s=XSteam('h_ps',Po(i),s7);
    h8=-((0.85*(h7-h8s))-h7);
    s8=XSteam('s_ph',Po(i),h8);

    h9s=XSteam('h_ps',Pc,s8);
    h9=-(0.85*(h8-h9s)-h8);

    x=(h3-h2)/(h8-h2);
    qrej=(1-x)*(h9-h1);
    qin=(h5-h4)+(h7-h8);
    Wnet(i)=qin-qrej;
    thermal_efficiency(i)=((1-qrej)/qin)*100;
else
    h2s=XSteam('h_ps',Po(i),s1);
    h2=((h2s-h1)/0.95)+h1;

    h3=XSteam('hL_p',Po(i));
    s3=XSteam('sL_p',Po(i));
    h4s=XSteam('h_ps',Pb,s3);
    h4=((h4s-h3)/0.95)+h3;
    
    h6s=XSteam('h_ps',Po(i),s5);
    h6=-(((h5-h6s)*0.85)-h5);
    s6=XSteam('s_ph',Po(i),h6);

    h7s=XSteam('h_ps',Pr,s6);
    h7=-(0.85*(h6-h7s)-h6);

    h8=XSteam('h_pT',Pr,Tmax);
    s8=XSteam('s_pT',Pr,Tmax);

    h9s=XSteam('h_ps',Pc,s8);
    h9=-(0.85*(h8-h9s)-h8);

    x=(h3-h2)/(h6-h2);
    qrej=(h9-h1);
    qin=(h5-h4)+(1-x)*(h8-h7);
    Wnet(i)=qin-qrej;
    thermal_efficiency(i)=((1-qrej)/qin)*100;

end
end
figure(1)
yyaxis left
plot(Po,thermal_efficiency)
yyaxis right
plot(Po,Wnet)
