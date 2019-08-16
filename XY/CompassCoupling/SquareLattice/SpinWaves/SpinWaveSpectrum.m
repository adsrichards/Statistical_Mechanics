%%
clear
%%
N = 100;

kx = pi*(-N:N)/N;
ky = pi*(-N:N)/N;

f2 = figure(2); clf; box on;
V = zeros(2*N+1,2*N+1,5);

% 1
th = 0;
A = zeros(2*N+1,2*N+1);
for ix=-N:N
    A(ix+N+1,:) = -cos(ky).' - (cos(kx(ix+N+1)).' - cos(ky).')*(cos(th)^2);
end
V(:,:,1) = A;
% 2
th = pi/6;
A = zeros(2*N+1,2*N+1);
for ix=-N:N
    A(ix+N+1,:) = -cos(ky).' - (cos(kx(ix+N+1)).' - cos(ky).')*(cos(th)^2);
end
V(:,:,2) = A;
% 3
th = pi/4;
A = zeros(2*N+1,2*N+1);
for ix=-N:N
    A(ix+N+1,:) = -cos(ky).' - (cos(kx(ix+N+1)).' - cos(ky).')*(cos(th)^2);
end
V(:,:,3) = A;
% 4
th = pi/3;
A = zeros(2*N+1,2*N+1);
for ix=-N:N
    A(ix+N+1,:) = -cos(ky).' - (cos(kx(ix+N+1)).' - cos(ky).')*(cos(th)^2);
end
V(:,:,4) = A;
% 5
th = pi/2;
A = zeros(2*N+1,2*N+1);
for ix=-N:N
    A(ix+N+1,:) = -cos(ky).' - (cos(kx(ix+N+1)).' - cos(ky).')*(cos(th)^2);
end
V(:,:,5) = A;

Vp = permute(V,[1,3,2]);

xslice = [1,2,3,4,5];
yslice = [];
zslice = [];
s = slice(1:5,kx,ky,Vp,xslice,yslice,zslice);

axis([1 5 -pi pi -pi pi]);

xticks([1,2,3,4,5])
xticklabels(["0","\pi/6","\pi/4","\pi/3","\pi/2"]);
xlabel("$\theta$","Interpreter","latex")

yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels(["-\pi","-\pi/2","0","\pi/2","\pi"]);
ylabel("$k_x$","Interpreter","latex")

zticks([-pi,-pi/2,0,pi/2,pi])
zticklabels(["-\pi","-\pi/2","0","\pi/2","\pi"]);
zlabel("$k_y$","Interpreter","latex")

shading interp
set(gca,"fontsize",24)

view(-13.5,20)
colorbar
colormap(bone);
%%
f2=figure(2);
set(f2,"papersize",[10 4]);
print(f2,"./Plots/dispersion.pdf","-dpdf")