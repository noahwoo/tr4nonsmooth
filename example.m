
function example

b = gen_img(32);
R = 0.2*randn(32); 
b = b + R;
imshow(b);

% we use identity matrix as A
% b is the image pixels with random noise

% generate the image
function A = gen_img(n) 
A = zeros(n,n);
% white 
for r=(n/2-5):(n/2-1)
  A = A + gen_circle(n, n/2, n/2, r, 1);
end

% a single cicle
function A = gen_circle(n, x0, y0, r, c)
A = zeros(n,n);
t = 0:0.1:2*pi;
len = length(t);
x = r*sin(t)+x0; 
y = r*cos(t)+y0;
x = round(x); y = round(y); 
for k=1:len
    A(x(k),y(k)) = c;
end