%% 
%1���ȱ任���������
clc,clear
a=[1 2 3;2 1 2;1 3 4];
ae=[a,eye(3)];
b=rref(ae);%�ѣ�A l E)��Ϊ�������
b(:,1:3)=[]%ɾȥBǰ����
c=inv(sym(a));%ֱ���������
%% 
%2��������
clc,clear
format bank
a=[1 -2 2 -1;2 -4 8 0;-2 4 -2 3;3 -6 0 -6];
a
b=rref(a);
b
r=rank(a);
r
%% 
%3�����η�����
clc,clear,format bank %��������ʾ��ʽ
a=[1 -2 3 -1;3 -1 5 -3;2 1 2 -2];
b=[1 2 3]';
r1=rank(a), r2=rank([a b]);
x=pinv(a)*b%��С���˽�  
%A������ʱ����ʹ��x=A\b
%% 
%4��С�������
clc,clear
format short
t=1:0.5:8;
x0=[1:8]';
y0=[8 12 7 14 15 16 18 21]';
xs=[x0,ones(8,1)];
ab1=xs\y0;
ab2=polyfit(x0,y0,1)%1�ζ���ʽ
cde=polyfit(x0,y0,2)%2�˶���ʽ
xs3=[log(x0),1./x0];
fg=xs3\y0
y1=1.7738*t+5.8929
y2=0.1131*t.^2++0.7560*t+7.5893
y3=8.5557*log(t)+7.6742/4
figure
scatter(x0,y0,'ro'),hold on;
plot(t,y1,'r-'),hold on;
plot(t,y2,'g--'),hold on;
plot(t,y3,'b-.'),hold on;
%% 
%5ϡ�����Ч��
clc,clear 
b=ones(1,999);
a1=4*eye(1000)+diag(b,1)+diag(b,-1);
c=[1:1000]';
a2=sparse(a1);
tic,x1=a1\c,toc;
tic,x2=a2\c,toc;
cha=sum(abs(x1-x2));
