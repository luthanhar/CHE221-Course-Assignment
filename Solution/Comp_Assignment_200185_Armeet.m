syms P1sat(T) P2sat(T) z12(T) z21(T) v1(T,x) v2(T,x) F(T) G(T) J(T) d(T) f(T)
P1sat(T)=10^(7.47429-1314.188/(T+186.5-273.15));
P2sat(T)=10^(7.07638-1571.005/(T+209.728-273.15));
z12(T)=(138.93/92.35)*exp(-389.1597/(1.9872*T)); 
z21(T)=(92.35/138.93)*exp(-691.3689/(1.9872*T));
v1(T,x)=(1/(x+z12(T)*(1-x)))*exp((1-x)*(z12(T)/(x+z12(T)*(1-x))-z21(T)/(z21(T)*x+(1-x))));
v2(T,x)=(1/(1-x+z21(T)*x))*exp(-x*(z12(T)/(x+z12(T)*(1-x))-z21(T)/(z21(T)*x+(1-x))));
Eq(T,x)=760/((x*( P1sat(T)/P2sat(T)*v1(T,x)-v2(T,x))+v2(T,x))); 

mi=100;  %maximum iteration
mae=0.001; %tolerance
%given experimental data
Te=[164.9,148.2,142.2,134.1,120.8,114.1,110.1,106.9,104.5,103,101.8,100.7,100.2,99.9,99.5];
xe=[0,0.04,0.06,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.94,0.96,1];
ye=[0,0.3660,0.470,0.5960,0.758,0.825,0.8610,0.89,0.912,0.927,0.942,0.962,0.975,0.982,1];


y=zeros(1,101); 
T=zeros(1,101);
n=1;

for x=0:0.01:1
Temp=x*99.5988+ 164.71*(1-x)+273; %for getting a good estimate of initial assumption of Equilibrium
                                   %Temperature of

[T(1,n), P2s]=fpi(Eq,x,mi,mae,Temp); %calling function to return Eqlbm Temp and P2sat(T eqlbm)
y(1,n)=1-P2s*(1-x)*double(v2(T(1,n),x))/760; %vapor phase mole fraction of component 1
n=n+1;
end



T=T-273.15;

figure(1)
title("Isobaric T vs X,Y Diagram at P=760mm Hg");
plot(0:0.01:1,T(1,:),"g-");
hold on;
plot(y,T(1,:),'r-');
plot(xe,Te,'k.');
plot(ye,Te,'b.');
xlabel('x,y');
ylabel('Temperature(Celsius)');
legend('T-x','T-y','Experimental data of x','Experimental data of y');
hold off;

figure(2)
title("X vs Y Plot");
plot(0:0.01:1,y(1,:),'b-');
hold on 
plot(xe,ye,'k.');
legend('Plot from Wilson Model','Experimental data');
xlabel('x');
ylabel('y');


%values of T and y from Wilson model
ya=[y(1),y(5),y(7),y(11:10:91),y(95),y(97),y(101)];
Ta=[T(1),T(5),T(7),T(11:10:91),T(95),T(97),T(101)];

%displaying relative error from given experimental data
fprintf( "At P=760mm Hg, Pure component T2sat= %f Celsius ,T1sat= %f Celsius\n",164.7161,99.5955 );
disp("The relative error(in %) from experimental data for y");
fprintf( " %f\n",(abs(ya-ye)./ye).*100);
disp("The relative error(in %) from experimental data for Equilibrium Temperature");
fprintf( "%f\n",(abs(Ta-Te)./Te).*100);


%% Code for fixed point iteration method to get equilibrium temperature
function [T P3]=fpi(Eq,x,mi,mae,Temp)
T=Temp;
P3=0;
ae=100;
for i=1:mi
    if(ae<mae)
    return
    end
    P2=Eq(T,x);
    P3=double(P2);
    T=1571.005/(7.07638-log10(P3))-209.728+273.15;
    ae=abs(T-Temp)/T;
    Temp=T;
end
end
