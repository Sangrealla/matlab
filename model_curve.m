clear
hold on
colormap([0 0 0])
t=-120:1/10:120;
y1=sqrt(3)/3*t;
y2=-sqrt(3)/3*t;
plot(t,y1,'k',t,y2,'k')
plot([0 0],[0 120],'k')
plot([0 120],[0 0],'k')
syms x y
a1=100+sqrt(6)/3*y;
a2=100-sqrt(6)/6*y-sqrt(2)/2*x;
a3=100-sqrt(6)/6*y+sqrt(2)/2*x;
%摩尔库伦准则
ezplot(a1-3*a3,[-300*sqrt(2)/7,0,100*sqrt(6)/7,40*sqrt(6)])  %a1>a2>a3
ezplot(a1-3*a2,[0,300*sqrt(2)/7,100*sqrt(6)/7,40*sqrt(6)])   %a1>a3>a2
ezplot(a2-3*a3,[-60*sqrt(2),-300*sqrt(2)/7,-20*sqrt(6),100*sqrt(6)/7])  %a2>a1>a3
ezplot(a2-3*a1,[-60*sqrt(2),0,-40*sqrt(6),-20*sqrt(6)])  %a2>a3>a1
ezplot(a3-3*a1,[0,60*sqrt(2),-40*sqrt(6),-20*sqrt(6)])  %a3>a2>a1
ezplot(a3-3*a2,[300*sqrt(2)/7,60*sqrt(2),-20*sqrt(6),100*sqrt(6)/7])  %a3>a1>a2
%SMP准则
ezplot(1/a1+1/a2+1/a3-7/180,[-150,150,-150,150])
%Lade准则
ezplot(1/(a1*a2*a3)-1/648000,[-150,150,-150,150])
title('P=100，\phi=30，C=0时的莫尔库伦准则,Lade准则,SMP准则在偏平面的破坏曲线')
axis square;
grid on;
legend('M-C','Lade','SMP');

