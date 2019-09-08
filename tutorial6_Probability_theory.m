%% 
%   1_����¼���
clc,clear,n=10000;
a=rand(1,n);
n1=sum(a<0.3),n2=sum(a>=0.3&a<0.8),n3=sum(a>=0.8);
f=[n1,n2,n3]/n
%% 
%   2_BuffonͶ�룻
% ����pi��ֵ��p=(2*l)/(a*pi)=m/n;  pi=2*n*l/(a*m)
clc,clear
a=45;L=36;n=5000;
x=unifrnd(0,a/2,1,n);%����[0,a/2]�Ͼ��ȷֲ��������
phi=unifrnd(0,pi,1,n);%����[0,pi]�Ͼ��ȷֲ��������
m=sum(x<=L*sin(phi)/2);
pis=2*n*L/(a*m)
%% 
%   3_��������
clc,clear,p1=0.99;p2=0.05;k=0:3;
PAH=0.99.^(3-k).*0.55.^k;
PHi=[];
for i=0:3
    PHi=[PHi,nchoosek(96,3-i)*nchoosek(4,i)/nchoosek(100,3)];
end
PA=dot(PAH,PHi);
%% 
%   4_����ֵlambda=3�Ĳ��ɷֲ�
clc,clear,lambda=3;
x0=0:20;
subplot(121),plot(x0,poisspdf(x0,lambda),'*-'),title('�ֲ���ͼ');
subplot(122),fplot(@(x)poisscdf(x,lambda),[0,20]),title('�ֲ�����ͼ');
%% 
%   5_����̬�ֲ�������a��λ��
clc,clear 
alpha=[0.001 0.005 0.01 0.025 0.05 0.10];
za=norminv(1-alpha);
fplot(@(x)normpdf(x),[-4,4]);
x0=[za(end):0.01:4];
y0=normpdf(x0);
xx0=[x0,4,za(end)];
yy0=[y0,0,0];
hold on ,fill(xx0,yy0,'g');
text(1.9,0.05,'\leftarrow\alpha=0.01');
text(za(end),0.02,['z_a=',num2str(za(end))]);
text(1.2,0.2,'\it X~N(0,1)');
%% 
%   6_test1. X~N(2,9);��P{1<X<5}��ȷ��cʹP{-c<X<2c}=0.8
clc,clear
p=normcdf(5,2,3)-normcdf(1,2,3);
fc=@(c)normcdf(2*c,2,3)-normcdf(-c,2,3)-0.8;
c1=fzero(fc,[0,12]);
c2=fsolve(fc,rand);
% test2.X~N(98.6,2),Y=5/9(X-32),��Y���ܶȺ���
syms x
y1=1/sqrt(4*sum(pi))*exp(-(x-98.6)^2/4);
y2=5/9*(x-32);
y3=finverse(y2);dy3=diff(y3);
y4=compose(y1,y3)*dy3;
y4=simplify(y4);
%% 
%   7_��Ϊ�������
clc,clear
fxy=@(x,y)1/8*(6-x-y).*(x>0&x<2).*(y>2&y<4).*(x+y<=4);
p=dblquad(fxy,0,2,2,4)%���ػ���
%% 
%   8_�����Ϸֲ���
clc,clear
X=[-3 0 1 3 5];
Y=[-2 0 1 2];
[X,Y]=meshgrid(X,Y);
P=[0.036 0.0198 0.0297 0.0209 0.0180
    0.0372 0.0558 0.0837 0.0589 0.0744
    0.0516 0.0774 0.1161 0.0817 0.1032
    0.0264 0.0270 0.0405 0.0285 0.0132];
PP=P.*(X.^2-3*Y>=1)
sump=sum(sum(PP));
%% 
%   9_��ά��̬�ֲ�N(1,4,2,4,0.5),��E(X)��E(Y)
clc,clear,syms f(x,y)
f=@(x,y)mvnpdf([x,y],[1,2],[4,2;2,4]);
EX=int(int(f*x,y,-inf,inf),x,-inf,inf);
disp(['EX=',num2str(double(EX))]);
EY=int(int(f*y,x,-inf,inf),y,-inf,inf);
disp(['EY=',num2str(double(EY))]);

