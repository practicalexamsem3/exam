**newton forward interpolation**
x=input('enter value of x:');
y=input('enter value of y:');
X=input('value to approx');
n=length(x);
d=zeros(n,n);
d(:,1)=y;
disp(d);
for j=2:n
    for i=1:(n-j+1)
        d(i,j)=d(i+1,j-1)-d(i,j-1);
    end
end
h=x(2)-x(1);
u=(X-x(1))/h;
A=y(1);
G=u;
for k=1:n-1
    A=A+G*d(1,k+1);
    G=G*(u-k)/k+1;
end
disp(A);
**langrange interpolation**
x=[1 ,3 ,5 ,8 ,10, 2];
y=[200 300 400 500 600 700];
X=2.5;
sum=0;
for i=1:6
    prodct=y(i);
    for j=1:6
        if i~=j
        prodct=prodct*((X-x(j)) / (x(i)-x(j)));
        end
    end
    sum=sum+prodct;
end
disp(sum);
**sarrus method**
A=[1,2,3;4,5,6;7,8,9];
x=0;
for i=1:3
    k=1;
    for j=1:3
        v=i+j-1;
        if(v>3)
            v=v-3;
        end
        k=k*A(j,v);
    end
    x=x+k;
end
disp(x);
y=0;
for i=1:3
    k=1;
    for j=1:3
        v=i+j-1;
        if(v>3)
            v=v-3;
        end
        k=k*A(j,4-v);
    end
    y=y+k;
end
disp(y);
res=x-y;
disp(res);

**upper triangular matrix**
a=[1,2,3;4,5,6;7,8,9];
a=[1 2 3;4 5 6;7 8 9];
for i=1:3
    for j=1:3
        if(i>j)
            p=-1*a(i,j)/a(j,i);
            for k=1:3
                a(i,k)=p*a(j,k)+a(i,k);
            end
        end
    end
end
disp(a);
**guass elimination**
A=[1,1,2;1,-3,2;2,-1,2];
B=[14;10;15];
X=[1;1;1];
Ag=[A,B];

for i=1:3
    for j=1:4
        if i>j
            p=-1*Ag(i,j)/Ag(j,j);
            for k=1:4
                Ag(i,k)=p*Ag(j,k)+Ag(i,k);
            end
        end
    end
end
disp(Ag);

X(3)=Ag(3,4)/Ag(3,3);
for i=2:-1:1
    sum=0;
    for j=i+1:3
        sum=sum+Ag(i,j)*X(j);
    end
    x(i)=(Ag(i,4)-sum)/Ag(i,i);
end


disp(p);
disp(X);
