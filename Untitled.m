clear all
[x,y,z]=meshgrid(-10:0.5:10,-10:0.5:10,-10:0.5:10);

f=(y.^2-z.^2-1);
[d,m]=isosurface(x,y,z,f,1);


grid on;view(3);axis equal;
xlabel('X'),ylabel('Y'),zlabel('Z');;