clear all
clc
load 'dataset/test.txt';
load 'dataset/train.txt';
new=[train;test];
y1=new(:,10);
x=new(:,2:9);
y1=y1-mean(y1);
for j=1:8
    x(:,j)=(x(:,j)-mean(x(:,j)))/std(x(:,j));
end
new=[x y1];
train=new(1:67,:);
test=new(68:97,:);
A=train(:,1:8);
bs=train(:,9);
val_in=test(:,1:8);
A=[A;val_in];
val_tar=test(:,9);
bs=[bs;val_tar];
save A.mat A
save bs.mat bs
for j=1:1
    j
b=bs(:,j);
save b.mat b
x0=zeros(8,1);
tic;
[T,X] = ode15s('right1',[0 150], x0);
Time=toc;
X=X';
[row col]=size(X);
y=X(:,col);
Id=find(abs(y)<=0.0001);
[row col]=size(Id);
save y.mat y
save row.mat row
end
error1=mean((val_in*y-val_tar).^2); 