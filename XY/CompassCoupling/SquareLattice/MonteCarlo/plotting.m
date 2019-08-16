
blue    = [1  , 0  , 98 ]/256;
teal    = [0  , 106, 110]/256;
green   = [2  , 120, 0  ]/256;
orange  = [127, 74 , 0  ]/256;
red     = [130, 0  , 0  ]/256;

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%

L = 4;
N = L*L;
data = load("data/L=4.txt");

runs = 80;
TN = 64;

Tmax = 0.5;
Tmin = 0.01;
T = Tmax - (Tmax-Tmin)*(0:(TN-1))/(TN-1);

data2 = zeros(TN,6); 
for Ti = 1:TN
    for i = 1:(runs*TN)
        if abs(data(i,1) - T(Ti)) < 0.001 
            data2(Ti,:) = data2(Ti,:) + data(i,:);
        end
    end
end
data2 = data2/(runs);


for Ti = 1:TN
    count = 0;
    for i = 1:length(data(:,1))
        if abs(data(i,1) - T(Ti)) < 0.001 
            data2(Ti,:) = data2(Ti,:) + data(i,:);
            count = count + 1;
        end
    end
    data2(Ti,:) = data2(Ti,:)/count;
end


data2(:,3) = N*(data2(:,3)-data2(:,2).^2)./data2(:,1);
data2(:,5) = N*(data2(:,5)-data2(:,4).^2)./(data2(:,1).^2);

disp(data2);

figure(1); clf; hold on; box on;
plot(data2(:,1),data2(:,2),"o","color",blue,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$m$","Interpreter","latex");
set(gca,"fontsize",24);

figure(2); clf; hold on; box on;
plot(data2(:,1),data2(:,3),"o","color",blue,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$\chi$","Interpreter","latex");
set(gca,"fontsize",24);

figure(3); clf; hold on; box on;
plot(data2(:,1),data2(:,4),"o","color",blue,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$u$","Interpreter","latex");
set(gca,"fontsize",24);

figure(4); clf; hold on; box on;
plot(data2(:,1),data2(:,5),"o","color",blue,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$c_V$","Interpreter","latex");
set(gca,"fontsize",24);

figure(5); clf; hold on; box on;
plot(data2(:,1),data2(:,6),"o","color",blue,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$m_L$","Interpreter","latex");
set(gca,"fontsize",24);


%%

runs = 40;

L = 8;
N = L*L;
data = load("data/L=8.txt");

data2 = zeros(TN,6);
for Ti = 1:TN
    for i = 1:runs*TN
        if abs(data(i,1) - T(Ti)) < 0.001 
            data2(Ti,:) = data2(Ti,:) + data(i,:);
        end
    end
end
data2 = data2/runs;

data2(:,3) = N*(data2(:,3)-data2(:,2).^2)./data2(:,1);
data2(:,5) = N*(data2(:,5)-data2(:,4).^2)./(data2(:,1).^2);

disp(data2);

figure(1); hold on; box on;
plot(data2(:,1),data2(:,2),"o","color",teal,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$q$","Interpreter","latex");
set(gca,"fontsize",24);

figure(2); hold on; box on;
plot(data2(:,1),data2(:,3),"o","color",teal,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$\chi$","Interpreter","latex");
set(gca,"fontsize",24);

figure(3); hold on; box on;
plot(data2(:,1),data2(:,4),"o","color",teal,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$u$","Interpreter","latex");
set(gca,"fontsize",24);

figure(4); hold on; box on;
plot(data2(:,1),data2(:,5),"o","color",teal,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$c_V$","Interpreter","latex");
set(gca,"fontsize",24);

figure(5); hold on; box on;
plot(data2(:,1),data2(:,6),"o","color",teal,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$m_L$","Interpreter","latex");
set(gca,"fontsize",24);

%%

runs = 80;

L = 16;
N = L*L;
data = load("data/L=16.txt");

data2 = zeros(TN,6);
for Ti = 1:TN
    for i = 1:runs*TN
        if abs(data(i,1) - T(Ti)) < 0.001 
            data2(Ti,:) = data2(Ti,:) + data(i,:);
        end
    end
end
data2 = data2/runs;

data2(:,3) = N*(data2(:,3)-data2(:,2).^2)./data2(:,1);
data2(:,5) = N*(data2(:,5)-data2(:,4).^2)./(data2(:,1).^2);

disp(data2);

figure(1); hold on; box on;
plot(data2(:,1),data2(:,2),"o","color",green,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$q$","Interpreter","latex");
set(gca,"fontsize",24);

figure(2); hold on; box on;
plot(data2(:,1),data2(:,3),"o","color",green,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$\chi$","Interpreter","latex");
set(gca,"fontsize",24);

figure(3); hold on; box on;
plot(data2(:,1),data2(:,4),"o","color",green,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$u$","Interpreter","latex");
set(gca,"fontsize",24);

figure(4); hold on; box on;
plot(data2(:,1),data2(:,5),"o","color",green,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$c_V$","Interpreter","latex");
set(gca,"fontsize",24);
axis([0.1 0.5 0.5 1.])

%%

runs = 10;

L = 32;
N = L*L;
data = load("data/L=32.txt");

data2 = zeros(TN,6);
for Ti = 1:TN
    for i = 1:(runs*TN)
        if abs(data(i,1) - T(Ti)) < 0.001 
            data2(Ti,:) = data2(Ti,:) + data(i,:);
        end
    end
end
data2 = data2/runs;

data2(:,3) = N*(data2(:,3)-data2(:,2).^2)./data2(:,1);
data2(:,5) = N*(data2(:,5)-data2(:,4).^2)./(data2(:,1).^2);

disp(data2);

figure(1); hold on; box on;
plot(data2(:,1),data2(:,2),"o","color",orange,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$q$","Interpreter","latex");
set(gca,"fontsize",24);

figure(2); hold on; box on;
plot(data2(:,1),data2(:,3),"o","color",orange,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$\chi$","Interpreter","latex");
set(gca,"fontsize",24);

figure(3); hold on; box on;
plot(data2(:,1),data2(:,4),"o","color",orange,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$u$","Interpreter","latex");
set(gca,"fontsize",24);

figure(4); hold on; box on;
plot(data2(:,1),data2(:,5),"o","color",orange,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$c_V$","Interpreter","latex");
set(gca,"fontsize",24);


%%

L = 64;
N = L*L;
data = load("data/L=64.txt");

data2 = zeros(TN,6);

for Ti = 1:TN
    count = 0;
    for i = 1:length(data(:,1))
        if abs(data(i,1) - T(Ti)) < 0.001 
            data2(Ti,:) = data2(Ti,:) + data(i,:);
            count = count + 1;
        end
    end
    data2 = data2/(count);
end

data2(:,3) = N*(data2(:,3)-data2(:,2).^2)./data2(:,1);
data2(:,5) = N*(data2(:,5)-data2(:,4).^2)./(data2(:,1).^2);

disp(data2);

figure(1); hold on; box on;
plot(data2(:,1),data2(:,2),"o","color",red,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$q$","Interpreter","latex");
set(gca,"fontsize",24);

figure(2); hold on; box on;
plot(data2(:,1),data2(:,3),"o","color",red,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$\chi$","Interpreter","latex");
set(gca,"fontsize",24);

figure(3); hold on; box on;
plot(data2(:,1),data2(:,4),"o","color",red,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$u$","Interpreter","latex");
set(gca,"fontsize",24);
legend("L=4","L=8","L=16","L=32","L=64","location","northwest")

figure(4); hold on; box on;
plot(data2(:,1),data2(:,5),"o","color",red,"LineWidth",1.5);
xlabel("$T$","Interpreter","latex");
ylabel("$c_V$","Interpreter","latex");
set(gca,"fontsize",24);
%%

f1 = figure(1);
set(f1,"papersize",[4.5 4])
print(f1,"./Plots/q.pdf","-dpdf")

f2 = figure(2);
set(f2,"papersize",[4.5 4])
print(f2,"./Plots/chi.pdf","-dpdf")

f3 = figure(3);
set(f3,"papersize",[4.5 4])
print(f3,"./Plots/u.pdf","-dpdf")

f4 = figure(4);
set(f4,"papersize",[4.5 4])
print(f4,"./Plots/cv.pdf","-dpdf")