%% 
%   1_������
clc,clear 
a=[5 2 -2;-1 4 3;2 6 5];
L1 = norm(a,1);     %1����/�кͷ���
Linf = norm(a,inf);     %�����/�кͷ���
Lf = norm(a,'fro');     %Frobenius����/(sigma(xi))^2
L2 = sqrt(eigs(a'*a,1));     %2���� /A'*A���������ֵ��ƽ��
L22 = max(svd(a));      %�����������ֵ
%% 
%   2_��������ֵ�ֽ�
%���Ż����⡢����ֵ���⡢��С�������⡢������������⼰ͳ������
clc,clear
a=[1 0 1;0 1 1;0 0 0];
a=sym(a);
sigma=svd(a);
[p,d,q]=svd(a);
%% 
%   3_ͼ��ѹ��
clc,clear
X=imread('2_media_screenshots_Dva01.jpg');
if (size(X,3)~=1),X=rgb2gray(X);
end
[U,V,D]=svd(double(X));
plot(diag(D),'b-','LineWidth',1.5);
title('ͼ����������ֵ')
ylabel('����ֵ');
[m,n]=size(X);
Rank=rank(double(X));
figure,subplot(1,2,1),imshow(X);
Image_Rank=['ͼ��������=',int2str(Rank)];
title(Image_Rank,'color','k');
it=0;
for K=100:100:Rank
    R=U(:,1:K)*D(1:K,1:K)*V(:,1:K)';
    T=uint8(R);subplot(1,2,2),imshow(T);
    SVD_number=['ѡȡ������ֵ����=',int2str(K)];
    title(SVD_number,'color','b')
    src_elements=m*n;compress_elements=K*(m+n+1);
    compress_ratio = (1-compress_elements/src_elements)*100;
    it=it+1;
    CR(it)=compress_ratio;
    fprint f ('Rank=% d:K=% d��:compress_ratio=%.2f\\n',Rank,K,compress_ratio);
    pause(2);
end
figure,plot((1:10:10*it),CR,'ob-','LineWidth',1.5);
title('����ֵ������ѹ�����ʵĹ�ϵ');
xlabel('����ֵ����');
ylabel('ѹ������');
%% 
%   4_��Ӧ����
% 10��ʡ��ũ���ͥ�˾�����֧���ṹѡȡ8��ָ�꣬ʳƷ֧��(x1)������֧��(x2)��...
% ��ס֧��(x5)����ͥ�豸������֧��(x4)����ͨ��ͨ��֧��(x5)���Ľ����ּ�������Ʒ֧��(x6)��...
% ҽ�Ʊ���֧��(x7)��������Ʒ������֧��(x8)��
%i�ӣ�1��10����j�ӣ�1��8����
% ���ʾ���pij=aij/a..;
% �任����bij=(pij-p.j*pi.)/(sqrt(p.j*pi.));
% Э�������R=(B^T)*B
%% 
%   5_���ȷֽ�  A=BC
clc,clear
a=[1 4 -1 2 3;2 0 0 0 -4;1 2 -4 -6 -10;2 6 3 12 17];
[r,ind]=rref(a);
b=a(:,ind);c=r(1:length(ind),:);
%% 
%   6_���������
% ����Moore-Penrose����1��AGA=A;2��GAG=G��3��(GA)^H=GA;4��(AG)^H=AG;
% A=BC,G=C^H*(C*C^H)^-1(B^H*B)^-1*B^H
clc,clear
a=[0 2i 1i 0 1+1i 1;0 0 0 -3 -2 -1-1i;0 2 1 1 1-1i 1];
aa=sym(a)
[r,ind]=rref(a),rr=rref(aa);
b=aa(:,ind),c=rr(1:length(ind),:);%���ȷֽ�
ap1=c'*inv(c*c')*inv(b'*b)*b;%��ʽ���α��
ap2=pinv(aa);%���������α��
%% 
