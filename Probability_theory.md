# 简单笔记
-------------

对数正态分布（$\mu%$,$\sigma^2$）矩估计值  
$\mu=2ln(\frac{1}{n}\sum\limits_{i=1}^n x_i)-\frac{1}{2}ln(\frac{1}{n}\sum\limits_{i=1}^n (x_i)^2)$  
$\sigma^2=ln(\frac{1}{n}\sum\limits_{i=1}^n (x_i)^2)-2ln(\frac{1}{n}\sum\limits_{i=1}^n x_i)$  
**************
## 随机数产生函数


|函数名称  | 函数说明|调用格式|
|:---------: | :--------:|:--------|
|betarnd | $\beta$分布的随机数 |R=betarnd(A,B,m,n)|
|binornd | 二项分布随机数|R=binornd(N,P,m,n)|
