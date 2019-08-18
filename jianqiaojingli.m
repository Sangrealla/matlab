%%%%单调加载程序

clear;
%%杨召焕 静三轴试验结果
%临界状态参数he初始应力状态
%Mc=1.3;Me=1.3;laimuda=0.162;kaipa=0.0243; 
%miu=0.3;e0=0.974;
%st3=300;st1=st3;jmax=-500; AAA=1;beta1=1;
%杨腾
Mc=1.3;Me=0.65;laimuda=0.162;kaipa=0.0243; 
miu=0.3;e0=0.974;
st3=50;st1=st3;jmax=500; beta1=1;


%%%%Li 静三轴
%Mc=0.772; Me=0.614;laimuda=0.173;kaipa=0.034;miu=0.30;e0=1.021;st1=500;st3=500;jmax=350;     
%AAA=0.75;beta1=4;

%%%Weald clay 静三轴
%Mc=0.95;Me=0.95;laimuda=0.093;kaipa=0.025;miu=0.3;e0=0.632;st1=207;st3=207;max=119;               
%cy=0; r=2.; R=0.50;   beta=0;  L=2.;numb=0;     uu=2.6;UU=120;

%Mc=1.1392; Me=0.8256;laimuda=0.0508;kaipa=0.0108;miu=0.15;e0=0.62;
%st1=393;st3=393;max=362.5; cy=0; beta=-7.2; r=2; 
%R=0.50; HL=4.9;%numb=0; uu=3.4;UU=60;


v0=1+e0;                    %重塑土三轴等压固结后比容
X0=beta1*v0/(laimuda-kaipa);%临时变量

dj=0.1;
djmax=0:dj:jmax;            %将jmax按照dj的增量步，将其划分成增量的表达形式
dqj=dj+zeros(1,jmax/dj+1);  %静力阶段  有效应力形式 偏应力增量

dpzj=dqj/3;                 %静力阶段  总应力形式 平均应力增量
qj=st1-st3+djmax;
pzj=(st1+djmax+2*st3)/3;    %静力阶段  总应力 平均应力全量

pj(1)=pzj(1);               %静力阶段  有效应力  平均应力初值
pc(1)=pj(1);                %初始固结压力

%在加载面（Loading surface）等向固结的参量pc
uwj(1)=0;       %定义循环加载时的累计孔压值,该值等于0
stn1j(1)=0;     %轴向应变初值
stn3j(1)=0;     %剪切应变初值

pvstnj(1)=0;    %塑性体积应变
pqstnj(1)=0;    %塑性剪切应变
pqstnabsj(1)=0; %用于计算塑性剪切应变增量的绝对值

evstnj(1)=0;    %弹性体积应变
eqstnj(1)=0;    %弹性剪切应变
tqstnj(1)=0;    %循环加载开始阶段的剪切应变全量初始值

for j=1:jmax/dj+1
       
        M=Mc;       
        
        p1(j)=pj(j);    
        q1(j)=qj(j);
                       
        %pFpp(j)=2*p1(j)-pc(j);           %定义边界面F对边界面上的p点的偏导数；
        pFpp(j)=M^2*(2*p1(j)-pc(j));
        pFpq(j)=2*q1(j);                              %定义边界面F对边界面上的q点的偏导数；
        pFppc(j)=-M^2*p1(j);
       
        %pFpvstn(j)=-X0*pc(j).*p1(j);                 %定义边界面F对塑性体积应变的偏导数            
        %kp(j)=-(pFpp(j).*pFpvstn(j));                %定义边界面上的塑性模量,
        
        kp(j)=-X0*(pFpp(j).*pFppc(j).*pc(j));         %定义边界面上的塑性模量,          
        K(j)=((1+e0)/kaipa)*p1(j);                    %定义弹性体积应变模量
        G(j)=(1.5*(1-2*miu)/(1+miu))*K(j);            %定义弹性剪切模量 
       
        dpj(j)=-(pFpp(j)*pFpq(j)*dqj(j))./(pFpp(j).^2+kp(j)./K(j));   
        duwj(j)=dpzj(j)-dpj(j);  
        
       
        L(j)= (pFpp(j)*dpj(j)+pFpq(j)*dqj(j))/kp(j);           %%计算加载因子，并对加载因子进行M括号运算
       
        dpc(j)=X0*pc(j)*L(j)*pFpp(j);                %等向硬化规则 
               
        dpvstnj(j)=L(j).*pFpp(j);                     %塑性体积应变增量
        dpqstnj(j)=L(j).*pFpq(j);                     %塑性剪切应变增量
        pvstnj(j+1)=pvstnj(j)+dpvstnj(j);             %确定塑性体积应变总量
        pqstnj(j+1)=pqstnj(j)+dpqstnj(j);             %确定塑性剪切应变总量
        pqstnabsj(j+1)=pqstnabsj(j)+abs(dpqstnj(j));  %塑性剪切应变长度
           
        devstnj(j)=dpj(j)./K(j);             %定义弹性体积应变增量
        deqstnj(j)=dqj(j)./(3*G(j));         %定义弹性剪切应变增量       
        evstnj(j+1)=evstnj(j)+devstnj(j);    %确定弹性体积应变总量
        eqstnj(j+1)=eqstnj(j)+deqstnj(j);    %确定弹性剪切应变总量
    
        dtvstnj(j)=dpvstnj(j)+devstnj(j);    %确定总的体积应变增量,如果是在不排水条件下，该值应该为体积应变增量应该等于零
        dtqstnj(j)=deqstnj(j)+dpqstnj(j);    %确定总的剪切应变增量，等于弹性剪切应变增量加上塑性剪切应变增量
        tqstnj(j+1)=tqstnj(j)+dtqstnj(j);    %利用总的剪切应变增量，计算出剪切应变全量
    
        dstn1j(j)=dtqstnj(j);               %计算轴向应变增量
        dstn3j(j)=-dstn1j(j)/2;             %计算周向应变增量
        stn1j(j+1)=stn1j(j)+dstn1j(j);       %计算轴向应变
        stn3j(j+1)=stn3j(j)+dstn3j(j);       %计算环向向应变
    
        uwj(j+1)=uwj(j)+duwj(j);          %定义空隙水压力的累计值
        pj(j+1)=pj(j)+dpj(j);             %定义有效应力路径中的p值
    
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
  
    %绘制应力路径，轴向应变、孔压-偏应力之间的关系
    %  plot(stn1jj,uwjj)  
    %plot(stn1jj,qjj,'r');
     plot(pjj,qjj,'r');
       
    hold on    
    grid on
    
    