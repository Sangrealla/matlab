%% 
%1 求逆矩阵inv(a)
clc,clear
a = [4 2 1; 1 1 0; 1 2 0];
a=sym(a);
b=inv(a)
%% 
%2 求伴随矩阵*1
clc,clear
a=[3 1 -1 2; -5 1 3 -4; 2 0 1 -1; 1 5 3 -3];
a=sym(a);
n=length(a); b=zeros(n);
for i =1:n
    for j=1:n
        hij=a;hij(i,:)=[];hij(:,j)=[];
        b(j,i)=(-1)^(i+j)*det(hij);
    end
end
disp('方阵a的伴随矩阵b如下所示：'),b
%% 
%3 求伴随矩阵*2(Hamilton-Cayley)
clc,clear
a=[3 1 -1 2; -5 1 3 -4; 2 0 1 -1; 1 5 3 -3];
n=length(a);p1= poly(a);%计算a的特征多项式(tE-A)
p2=p1(1:end-1)
b=(-1)^(n-1)*polyvalm(p2,a)
%% 
%4矩阵变换
clc,clear
A=sym('a%d%d',3)
A=rot90(A)
%% 
% 5对称变换
clc,clear
axis equal
syms x0 y0 X Y t
eq1=(Y+y0)/2==3*(X+x0)/2+5;
eq2=3*(Y-y0)==-(X-x0);
[X Y]=solve(eq1,eq2,X,Y)
X=subs(X,{x0 y0},{1+cos(t),sin(t)});
Y=subs(Y,{x0 y0},{1+cos(t),sin(t)});
t0=0:0.02:2*pi;
x=1+cos(t0);
y=sin(t0);
plot(x,y,'k'),hold on,fplot(@ (x)3*x+5, [-5 2.5],'k')%对称轴y=3x+5
X=double(subs(X,t0));
Y=double(subs(Y,t0));
plot(X,Y,'k');%画镜像圆
text(1,2,'\downarrow原来的圆','FontSize',12)
text(-4,0,'\uparrow镜像圆','FontSize',12)
grid on
%% 
% 6三维透视
clc,clear
a=[3 1 5;5 1 5;5 0 5;3 0 5;3 1 4;5 1 4;5 0 4;3 0 4]';
b=[a;ones(1,length(a))]
b=sym(b)
P=eye(4);
P(3,3)=0;
P(4,3)=-1/10;
P=sym(P)
c=P*b
cc=c([1:3],:)./repmat(c(4,:),3,1)
%% 
%7 图像空间变换
clc,clear
a=imread('peppers.png');
tf1=affine2d([cosd(30),-sind(30),0;sind(30),cosd(30),0;0 0 1;]);
ta1=imwarp(a,tf1);
tf2=affine2d([5 0 0;0 10.5 0;0 0 1]);
ta2=imwarp(a,tf2);
subplot(131),imshow(a);
subplot(132),imshow(ta1);
subplot(133),imshow(ta2);