%% 
%   1_矩阵范数
clc,clear 
a=[5 2 -2;-1 4 3;2 6 5];
L1 = norm(a,1);     %1范数/列和范数
Linf = norm(a,inf);     %无穷范数/行和范数
Lf = norm(a,'fro');     %Frobenius范数/(sigma(xi))^2
L2 = sqrt(eigs(a'*a,1));     %2范数 /A'*A的最大特征值开平方
L22 = max(svd(a));      %矩阵最大奇异值
%% 
%   2_矩阵奇异值分解
%最优化问题、特征值问题、最小二乘问题、广义逆矩阵问题及统计问题
clc,clear
a=[1 0 1;0 1 1;0 0 0];
a=sym(a);
sigma=svd(a);
[p,d,q]=svd(a);
%% 
%   3_图像压缩
clc,clear
X=imread('2_media_screenshots_Dva01.jpg');
if (size(X,3)~=1),X=rgb2gray(X);
end
[U,V,D]=svd(double(X));
plot(diag(D),'b-','LineWidth',1.5);
title('图像矩阵的奇异值')
ylabel('奇异值');
[m,n]=size(X);
Rank=rank(double(X));
figure,subplot(1,2,1),imshow(X);
Image_Rank=['图像矩阵的秩=',int2str(Rank)];
title(Image_Rank,'color','k');
it=0;
for K=100:100:Rank
    R=U(:,1:K)*D(1:K,1:K)*V(:,1:K)';
    T=uint8(R);subplot(1,2,2),imshow(T);
    SVD_number=['选取的奇异值个数=',int2str(K)];
    title(SVD_number,'color','b')
    src_elements=m*n;compress_elements=K*(m+n+1);
    compress_ratio = (1-compress_elements/src_elements)*100;
    it=it+1;
    CR(it)=compress_ratio;
    fprint f ('Rank=% d:K=% d个:compress_ratio=%.2f\\n',Rank,K,compress_ratio);
    pause(2);
end
figure,plot((1:10:10*it),CR,'ob-','LineWidth',1.5);
title('奇异值个数与压缩比率的关系');
xlabel('奇异值个数');
ylabel('压缩比率');
%% 
%   4_对应分析
% 10个省市农村家庭人均消费支出结构选取8项指标，食品支出(x1)，衣着支出(x2)，...
% 居住支出(x5)，家庭设备及服务支出(x4)，交通和通信支出(x5)，文教娱乐及服务用品支出(x6)，...
% 医疗保健支出(x7)，其他商品及服务支出(x8)。
%i从（1，10），j从（1，8）；
% 概率矩阵pij=aij/a..;
% 变换矩阵bij=(pij-p.j*pi.)/(sqrt(p.j*pi.));
% 协方差矩阵R=(B^T)*B
%% 
%   5_满秩分解  A=BC
clc,clear
a=[1 4 -1 2 3;2 0 0 0 -4;1 2 -4 -6 -10;2 6 3 12 17];
[r,ind]=rref(a);
b=a(:,ind);c=r(1:length(ind),:);
%% 
%   6_广义逆矩阵
% 满足Moore-Penrose方程1、AGA=A;2、GAG=G；3、(GA)^H=GA;4、(AG)^H=AG;
% A=BC,G=C^H*(C*C^H)^-1(B^H*B)^-1*B^H
clc,clear
a=[0 2i 1i 0 1+1i 1;0 0 0 -3 -2 -1-1i;0 2 1 1 1-1i 1];
aa=sym(a)
[r,ind]=rref(a),rr=rref(aa);
b=aa(:,ind),c=rr(1:length(ind),:);%满秩分解
ap1=c'*inv(c*c')*inv(b'*b)*b;%公式求解伪逆
ap2=pinv(aa);%工具箱求解伪逆
%% 
