%% 
%   1_随机事件；
clc,clear,n=10000;
a=rand(1,n);
n1=sum(a<0.3);n2=sum(a>=0.3&a<0.8);n3=sum(a>=0.8);
f=[n1,n2,n3]/n;
%% 
%   2_Buffon投针；
% 计算pi的值，p=(2*l)/(a*pi)=m/n;  pi=2*n*l/(a*m)
clc,clear
a=45;L=36;n=5000;
x=unifrnd(0,a/2,1,n);%产生[0,a/2]上均匀分布随机数；
phi=unifrnd(0,pi,1,n);%产生[0,pi]上均匀分布随机数；
m=sum(x<=L*sin(phi)/2);
pis=2*n*L/(a*m);
%% 
%   3_条件概率
clc,clear,p1=0.99;p2=0.05;k=0:3;
PAH=0.99.^(3-k).*0.55.^k;
PHi=[];
for i=0:3
    PHi=[PHi,nchoosek(96,3-i)*nchoosek(4,i)/nchoosek(100,3)];
end
PA=dot(PAH,PHi);
%% 
%   4_画均值lambda=3的泊松分布
clc,clear,lambda=3;
x0=0:20;
subplot(121),plot(x0,poisspdf(x0,lambda),'*-'),title('分布律图');
subplot(122),fplot(@(x)poisscdf(x,lambda),[0,20]),title('分布函数图');
%% 
%   5_画正态分布，及上a分位点
clc,clear 
alpha=[0.001 0.005 0.01 0.025 0.05 0.10];
za=norminv(1-alpha);
fplot(@(x)normpdf(x),[-4,4]);
x0=za(end):0.01:4;
y0=normpdf(x0);
xx0=[x0,4,za(end)];
yy0=[y0,0,0];
hold on ,fill(xx0,yy0,'g');
text(1.9,0.05,'\leftarrow\alpha=0.01');
text(za(end),0.02,['z_a=',num2str(za(end))]);
text(1.2,0.2,'\it X~N(0,1)');
%% 
%   6_test1. X~N(2,9);求P{1<X<5}；确定c使P{-c<X<2c}=0.8
clc,clear
p=normcdf(5,2,3)-normcdf(1,2,3);
fc=@(c)normcdf(2*c,2,3)-normcdf(-c,2,3)-0.8;
c1=fzero(fc,[0,12]);
c2=fsolve(fc,rand);
% test2.X~N(98.6,2),Y=5/9(X-32),求Y的密度函数
syms x
y1=1/sqrt(4*sum(pi))*exp(-(x-98.6)^2/4);
y2=5/9*(x-32);
y3=finverse(y2);dy3=diff(y3);
y4=compose(y1,y3)*dy3;
y4=simplify(y4);
%% 
%   7_多为随机变量
clc,clear
fxy=@(x,y)1/8*(6-x-y).*(x>0&x<2).*(y>2&y<4).*(x+y<=4);
p=integral2(fxy,0,2,2,4);%二重积分
%% 
%   8_求联合分布律
clc,clear
X=[-3 0 1 3 5];
Y=[-2 0 1 2];
[X,Y]=meshgrid(X,Y);
P=[0.036 0.0198 0.0297 0.0209 0.0180
    0.0372 0.0558 0.0837 0.0589 0.0744
    0.0516 0.0774 0.1161 0.0817 0.1032
    0.0264 0.0270 0.0405 0.0285 0.0132];
PP=P.*(X.^2-3*Y>=1);
sump=sum(sum(PP));
%% 
%   9_二维正态分布N(1,4,2,4,0.5),求E(X)、E(Y)
clc,clear,syms f(x,y)
f=@(x,y)mvnpdf([x,y],[1,2],[4,2;2,4]);
EX=int(int(f*x,y,-inf,inf),x,-inf,inf);
disp(['EX=',num2str(double(EX))]);
EY=int(int(f*y,x,-inf,inf),y,-inf,inf);
disp(['EY=',num2str(double(EY))]);
%% 
%   10_大数定律和中心极限定理
%   1、大数定律指出，n→inf，随机变量的算数平均值依概率收敛到数学期望a；
%   2、中心极限定理指出，X~B(n，p)，xi近似服从N（np，np（1-p））
clc,clear,format long
n=90000;p=1/3;
mu=n*p;sigma=sqrt(n*p*(1-p));
p1=binocdf(30500,n,p)-binocdf(29499,n,p);
p2=normcdf(30500,mu,sigma)-normcdf(29500,mu,sigma);
delta_p=p1-p2;format
%% 
%   11_参数估计
clc,clear
data=[2 24 32 24 10 6 3 1 1];
bin=200:400:3400;
subplot(121),bar(bin,data/100);
mu1=2*log(dot(data,bin)/100)-0.5*log(dot(data,bin.^2)/100);
var1=log(dot(data,bin.^2)/100)-2*log(dot(data,bin)/100);
mu2=dot(data,log(bin))/100;
var2=dot((log(bin)-mu2).^2,data)/100;
subplot(122),x2=0:3600;y2=lognpdf(x2,mu1,sqrt(var1));
plot(x2,y2);
%% 
%   12_求mu的置信水平为0.95的置信区间
x0=[129 134 114 120 116 133 142 138 148 129 133 141 142];
alpha=0.05;mu=mean(x0);sigma=std(x0);n=length(x0);
muci1=[mu-sigma/sqrt(n)*tinv(1-alpha/2,n-1),mu+sigma/sqrt(n)*tinv(1-alpha/2,n-1)];
[mu2,sigma2,muci2]=normfit(x0,alpha);
[ms3,msci3]=mle(x0,'distribution','norm','alpha',alpha);
%% 
%   13_假设检验
