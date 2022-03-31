function result=SchafferF6(x0)
%Schaffer 函数
%输入x,给出相应的y值,在x=(0,0,…,0) 处有全局最小点.
[row,col]=size(x0);
if row>1
    error('输入的参数错误');
end
x=x0(1,1);
y=x0(1,2);
temp=x^2+y^2;
result=0.5+((sin(sqrt(temp)))^2-0.5)/(1+0.001*temp)^2;