%% 
clc,clear
% 1逆序数(t=t1+t2+...+tn)
a=32514;
b=num2str(a);
for i = 1:length(b);
    c(i)=str2num(b(i));
end
s=0;
for i =2:length(b)
    s=s+sum(c(1:i-1)-c(i)>0);
end
s
%% 
clc,clear
%2计算行列式的值
A1=[2 1 4 1; 3 -1 2 1; 1 2 3 2; 5 0 6 2];
val1=det(A1);
B=sym(A1),val2=det(B);
syms a b c d 
A2=[a b c; b c a; c a b]
val=det(A2)
%% 
% 3计算（-2,-2)(0,3)(4,-1)(6,4)确定平行四边形面积
clc,clear
a=[-2 -2; 0 3; 4 -1; 6 4];
b=a-repmat(a(1,:),size(a,1),1)%把其中一个顶点平移到坐标原点
s=abs(det(b([2,3],:)))
%%
%4计算(x1^2/a^2)+(x2^2/b^2)=1的面积
clc,clear
syms a b u1 u2 x1 x2
u=[u1,u2]';x=[x1,x2]'
a=[a 0; 0 b]
s=det(a)*pi
%% 
clc,clear
% 5解线性方程组
a=[2 1 -5 1; 1 -3 0 -6; 0 2 -1 2; 1 4 -7 6];
b=[8 9 -5 0]';
a1=a; a1(:,1)=b;
a2=a; a2(:,2)=b; 
a3=a; a3(:,3)=b; 
a4=a; a4(:,4)=b;
for i = 1:4
    str=['x',int2str(i),'=det(a',int2str(i),')/det(a)']
    %上面是构造xi=det(ai)/det(a)的字符串
    eval(str)   %执行字符串对应的命令
end
%% 
% 6计算特征值lambda
clc,clear,syms t
a=[5-t 2 2; 2 6-t 0; 2 0 4-t];
D=det(a),DD=factor(D)
s=solve(D)