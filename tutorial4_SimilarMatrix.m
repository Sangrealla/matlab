%% 
%   1_Schmidt正交化
clc,clear
a=[1 -1 4; 2 3 -1; -1 1 0];
a=sym(a);
[m, n]=size(a);
e(:,1)=a(:,1)/norm(a(:,1));
for j=2:n
    bj=a(:,j);
    for i=1:j-1
        bj=bj-(e(:,i)'*bj)*e(:,i);
    end
    e(:,j)=bj/norm(bj);
end
e;
%% 
% [v,d] = eig(A);  返回v为特征向量，d为对角线为A的特征值的对角阵；
% p=poly(A); 返回A的特征多项式f(A)=A^n+a1*A^n-1+...+an*E=0
% trace(A);  返回A的迹，
%   2_计算特征值和特征向量
clc,clear 
a=[-1 1 0; -4 3 0; 1 0 2];
p=charpoly(a);
r=roots(p); %特征多项式的根（特征值）
[vec,val]=eig(a);
%% 
%   3_(p^-1)*A*P=B,求P
clc,clear,syms x y 
A=[1 -2 -4; -2 x -2; -4 -2 1]; B = diag([5 -4 y]);
eq1=det(A)==det(B);
eq2=trace(A)==trace(B);
[x0,y0]=solve(eq1,eq2);
A=subs(A,{x y},{x0 y0});
[v1,d1] = eig(double(A));   %数值矩阵的特征值（向量）
[v2,d2] = eig(A);
v2 = myschmidth(v2);
%% 
%   4_矩阵多项式的值
clc,clear
f1=@(x,k)(factorial(10)/factorial(10-k)*x.^(10-k)-6*factorial(9)/factorial(9-k)...
    *x.^(9-k)+5*factorial(8)/factorial(9-k)*x.^(8-k))*(k<=8)+...
    (factorial(10)*x-6*factorial(9))*(k==9)+factorial(10)*(k==10);
f2=@(x)x^10-6*x^9+5*x^8;%定义多项式的匿名函数，效率最高
A=[2 1 2;1 2 2;2 2 1];
tic,B1 = funm(A,f1);toc;
tic,B2 = A^10-6*A^9+5*A^8;toc;
tic,B3 = f2(A);toc;
% fumn(A,fun)计算矩阵函数，一般要求函数具有收敛半径r=+\infty 的泰勒级数，...
% 函数fun的返回值是函数的各阶导数,对于本例中的函数\phi（x）=x^10-6x^9+5x^8
%% 
%   5_求最大特征值及特征向量
%[V,D] = eigs(A,k,sigma)
% sigma='lm',求前k个；sigma='sm',求后k个
clc,clear
a=[3 1 2;1 6 -2;2 -2 1];
[v1,d1]=eigs(a,1);       %求最大特征值及特征向量
[v2,d2]=eigs(a,1,'sm');%求最小特征值及特征向量
%另一种方法
[v,d]=eig(a);
vv1=v*[0 0 1]';
vv2=v*[1 0 0]';
%% 
% 层次分析法（the Analytic Hierarchy Process,AHP）
clc,clear 
a=[1 1 1 4 1 1/2;...
    1 1 2 4 1 1/2;...
    1 1/2 1 5 3 1/2;...
    1/4 1/4 1/5 1 1/3 1/3;...
    1 1 1/3 3 1 1;...
    2 2 2 3 1 1;];
[v1,d1] = eigs(a,1);%求A最大特征值及特征向量
B1=v1/sum(v1);%归一化
a1=[1 1/4 1/2; 4 1 3; 2 1/3 1];
a2=[1 1/4 1/5; 4 1 1/2; 5 2 1];
a3=[1 3 1/3; 1/3 1 1; 3 1 1];
a4=[1 1/4 1/2 ; 4 1 3; 2 1/3 1];
a5=[1 1/4 1/5; 4 1 1/2; 5 2 1];
a6=[1 3 1/3; 1/3 1 1; 3 1 1];
lambda = []; B2=[];
for i = 1:6
    str=['[v,d]=eigs(a',int2str(i),',',int2str(1),');v=v/sum(v);']
    eval(str)
    lambda = [lambda,d];B2=[B2,v]
end
lambda,B2;
B3=B2*B1;
%% 
% 马尔科夫链（Markov Chain）

