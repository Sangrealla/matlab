%%%%�������س���

clear;
%%���ٻ� ������������
%�ٽ�״̬����he��ʼӦ��״̬
%Mc=1.3;Me=1.3;laimuda=0.162;kaipa=0.0243; 
%miu=0.3;e0=0.974;
%st3=300;st1=st3;jmax=-500; AAA=1;beta1=1;
%����
Mc=1.3;Me=0.65;laimuda=0.162;kaipa=0.0243; 
miu=0.3;e0=0.974;
st3=50;st1=st3;jmax=500; beta1=1;


%%%%Li ������
%Mc=0.772; Me=0.614;laimuda=0.173;kaipa=0.034;miu=0.30;e0=1.021;st1=500;st3=500;jmax=350;     
%AAA=0.75;beta1=4;

%%%Weald clay ������
%Mc=0.95;Me=0.95;laimuda=0.093;kaipa=0.025;miu=0.3;e0=0.632;st1=207;st3=207;max=119;               
%cy=0; r=2.; R=0.50;   beta=0;  L=2.;numb=0;     uu=2.6;UU=120;

%Mc=1.1392; Me=0.8256;laimuda=0.0508;kaipa=0.0108;miu=0.15;e0=0.62;
%st1=393;st3=393;max=362.5; cy=0; beta=-7.2; r=2; 
%R=0.50; HL=4.9;%numb=0; uu=3.4;UU=60;


v0=1+e0;                    %�����������ѹ�̽�����
X0=beta1*v0/(laimuda-kaipa);%��ʱ����

dj=0.1;
djmax=0:dj:jmax;            %��jmax����dj�������������仮�ֳ������ı����ʽ
dqj=dj+zeros(1,jmax/dj+1);  %�����׶�  ��ЧӦ����ʽ ƫӦ������

dpzj=dqj/3;                 %�����׶�  ��Ӧ����ʽ ƽ��Ӧ������
qj=st1-st3+djmax;
pzj=(st1+djmax+2*st3)/3;    %�����׶�  ��Ӧ�� ƽ��Ӧ��ȫ��

pj(1)=pzj(1);               %�����׶�  ��ЧӦ��  ƽ��Ӧ����ֵ
pc(1)=pj(1);                %��ʼ�̽�ѹ��

%�ڼ����棨Loading surface������̽�Ĳ���pc
uwj(1)=0;       %����ѭ������ʱ���ۼƿ�ѹֵ,��ֵ����0
stn1j(1)=0;     %����Ӧ���ֵ
stn3j(1)=0;     %����Ӧ���ֵ

pvstnj(1)=0;    %�������Ӧ��
pqstnj(1)=0;    %���Լ���Ӧ��
pqstnabsj(1)=0; %���ڼ������Լ���Ӧ�������ľ���ֵ

evstnj(1)=0;    %�������Ӧ��
eqstnj(1)=0;    %���Լ���Ӧ��
tqstnj(1)=0;    %ѭ�����ؿ�ʼ�׶εļ���Ӧ��ȫ����ʼֵ

for j=1:jmax/dj+1
       
        M=Mc;       
        
        p1(j)=pj(j);    
        q1(j)=qj(j);
                       
        %pFpp(j)=2*p1(j)-pc(j);           %����߽���F�Ա߽����ϵ�p���ƫ������
        pFpp(j)=M^2*(2*p1(j)-pc(j));
        pFpq(j)=2*q1(j);                              %����߽���F�Ա߽����ϵ�q���ƫ������
        pFppc(j)=-M^2*p1(j);
       
        %pFpvstn(j)=-X0*pc(j).*p1(j);                 %����߽���F���������Ӧ���ƫ����            
        %kp(j)=-(pFpp(j).*pFpvstn(j));                %����߽����ϵ�����ģ��,
        
        kp(j)=-X0*(pFpp(j).*pFppc(j).*pc(j));         %����߽����ϵ�����ģ��,          
        K(j)=((1+e0)/kaipa)*p1(j);                    %���嵯�����Ӧ��ģ��
        G(j)=(1.5*(1-2*miu)/(1+miu))*K(j);            %���嵯�Լ���ģ�� 
       
        dpj(j)=-(pFpp(j)*pFpq(j)*dqj(j))./(pFpp(j).^2+kp(j)./K(j));   
        duwj(j)=dpzj(j)-dpj(j);  
        
       
        L(j)= (pFpp(j)*dpj(j)+pFpq(j)*dqj(j))/kp(j);           %%����������ӣ����Լ������ӽ���M��������
       
        dpc(j)=X0*pc(j)*L(j)*pFpp(j);                %����Ӳ������ 
               
        dpvstnj(j)=L(j).*pFpp(j);                     %�������Ӧ������
        dpqstnj(j)=L(j).*pFpq(j);                     %���Լ���Ӧ������
        pvstnj(j+1)=pvstnj(j)+dpvstnj(j);             %ȷ���������Ӧ������
        pqstnj(j+1)=pqstnj(j)+dpqstnj(j);             %ȷ�����Լ���Ӧ������
        pqstnabsj(j+1)=pqstnabsj(j)+abs(dpqstnj(j));  %���Լ���Ӧ�䳤��
           
        devstnj(j)=dpj(j)./K(j);             %���嵯�����Ӧ������
        deqstnj(j)=dqj(j)./(3*G(j));         %���嵯�Լ���Ӧ������       
        evstnj(j+1)=evstnj(j)+devstnj(j);    %ȷ���������Ӧ������
        eqstnj(j+1)=eqstnj(j)+deqstnj(j);    %ȷ�����Լ���Ӧ������
    
        dtvstnj(j)=dpvstnj(j)+devstnj(j);    %ȷ���ܵ����Ӧ������,������ڲ���ˮ�����£���ֵӦ��Ϊ���Ӧ������Ӧ�õ�����
        dtqstnj(j)=deqstnj(j)+dpqstnj(j);    %ȷ���ܵļ���Ӧ�����������ڵ��Լ���Ӧ�������������Լ���Ӧ������
        tqstnj(j+1)=tqstnj(j)+dtqstnj(j);    %�����ܵļ���Ӧ�����������������Ӧ��ȫ��
    
        dstn1j(j)=dtqstnj(j);               %��������Ӧ������
        dstn3j(j)=-dstn1j(j)/2;             %��������Ӧ������
        stn1j(j+1)=stn1j(j)+dstn1j(j);       %��������Ӧ��
        stn3j(j+1)=stn3j(j)+dstn3j(j);       %���㻷����Ӧ��
    
        uwj(j+1)=uwj(j)+duwj(j);          %�����϶ˮѹ�����ۼ�ֵ
        pj(j+1)=pj(j)+dpj(j);             %������ЧӦ��·���е�pֵ
    
        pc(j+1)=pc(j)+dpc(j);
       

    if stn1j(j)>0.15;
         break
     end 
  
         if kp(j)<0.0000;
         break
     end 
   
end
for i=1:j
   qjj(i)=qj(i);
   pjj(i)=pj(i); 
   stn1jj(i)=stn1j(i);
   uwjj(i)=uwj(i);   
end
    qjj=qjj';
    pjj=pjj';
    stn1jj=100*stn1jj';    
    uwjj=(uwjj)';   
  
    %����Ӧ��·��������Ӧ�䡢��ѹ-ƫӦ��֮��Ĺ�ϵ
    %  plot(stn1jj,uwjj)  
    %plot(stn1jj,qjj,'r');
     plot(pjj,qjj,'r');
       
    hold on    
    grid on
    
    