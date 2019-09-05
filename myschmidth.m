function P = myschmidth(a);
a=sym(a);
[m,n]=size(a);
e(:,1)=a(:,1)/norm(a(:,1));
for j = 2:n;
    bj=a(:,j);
    for i = 1:j-1
        bj=bj-(e(:,1)'*bj)*e(:,i);
    end
    e(:,j)=bj/norm(bj);
end
P = simplify(e)


