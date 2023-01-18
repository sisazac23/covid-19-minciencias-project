clc
clear all
format long
% 
q=[-0.6,-0.70,-0.80,-0.90,-0.99];
Nq=20;
for i=1:length(q)
for j=1:Nq
    eq(i,j)=(1+(1-q(i))*(j-1))^(1/(1-q(i)));
end

plot(0:1:Nq-1,eq(i,:),'LineWidth',2)
hold on

end
legend('q = -0.6','q = -0.7','q = -0.80','q = -0.9','q = -0.99')

% ______________________________________________

% Nq=20;
% dt=1/365;
% q=0.9;
% m=100;
% for j=1:Nq
%     eq(j)=(1+(1-q)*m*(j-1)*dt)^(1/(1-q));
% end
% 
% plot(0:1:Nq-1,eq)
% 
% hold on

%_______________________________________


%  y=[1;1;1;3;3;9;13;16;24;45;57;75;102;128;158;210;...
%    240;306;378];  
% % 
% N=length(y);
% % 
%  plot(0:1:N-1,y)
% 
% Nq=N;
% dt=1/365;
% q=0.999999;
% f=1.5;
% m=-94;
% % 
% for j=1:Nq
%      eq(j)=f*(1+(1-q)*m*(j-1)*dt)^(1/(1-q));
%  end
% 
% plot(0:1:N-1,y,0:1:Nq-1,eq)
% 
% hold on

%______________________________________

% 
% y=[1;1;1;3;3;9;13;16;24;45;57;75;102;128;158;210;...
%     240;306;378];  
% 
% N=length(y);
% dt=1/365;
% d=0:dt:dt*(N-1);
% q=0.9;
% fun= @(x)(1+(1-q)*x(1)*d).^(1/(1-q))-y;
% x0 = [100];
% [x,resnorm] = lsqnonlin(fun,x0);
% x,resnorm;
% f=2.64;
% for j=1:N
%     eq(j)=f*(1+(1-q)*x(1,1)*(j-1)*dt)^(1/(1-q));
% end
% 
% plot(0:1:N-1,y,0:1:N-1,eq)
% 
% hold on


%____________________________________________

y=[1;1;1;3;3;9;13;16;24;45;57;75;102;128;158;210;...
     240;306;378;470]; 
%Uk data
cases=[2
2
2
2
2
2
2
3
3
4
4
8
8
9
9
9
9
9
9
9
9
9
9
9
13
13
13
13
16
18
23
36
40
51
85
115
163
206
273
321
373
456
590
707
1140
1391
1543
1950
2630
3277
3983
5018
5683
6650
8077
9529
11658
]

% Spain data
% cases=[1
% 1
% 1
% 1
% 1
% 1
% 1
% 1
% 1
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 2
% 3
% 7
% 12
% 25
% 34
% 66
% 83
% 114
% 151
% 200
% 261
% 374
% 430
% 589
% 1204
% 1639
% 2140
% 3004
% 4231
% 5753
% 7753
% 9191
% 11178
% 13716
% 17147
% 19980
% 24926
% 28572
% 33089
% 39673
% 47610
%  ]
% Usa
% y=[0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 0
% 1
% 1
% 1
% 1
% 2
% 2
% 5
% 5
% 5
% 5
% 6
% 7
% 8
% 11
% 11
% 11
% 12
% 12
% 12
% 12
% 12
% 13
% 13
% 14
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 16
% 35
% 35
% 35
% 53
% 53
% 59
% 60
% 66
% 69
% 89
% 103
% 125
% 159
% 233
% 338
% 433
% 554
% 754
% 1025
% 1312
% 1663
% 2174
% 2951
% 3774
% 4661
% 6427
% 9415
% 14250
% 19624
% 26747
% 35206
% 46442
% 55231
% 69194
% 85991
% ];
day=52;
data=[];

for i=day-1:length(cases)-1
    d=cat(1,cases(1:i,1),(zeros(1,30))');
    d=d(1:length(cases),:)';
data=cat(1,data,d);
end
[a b]=size(data);
dataTable=[];
for i=1:a
y=data(i,:);
y=y(y~=0)
y=y'

N=length(y);
q=0.8:0.01:0.99;
nq=length(q);
dt=1/365;
d=0:dt:dt*(N-1); 

for h=1:nq
    fun= @(x)(1+(1-q(h))*x(1)*d).^(1/(1-q(h)))-y;
    x0 = [100];
    [x(h),resnorm] = lsqnonlin(fun,x0);
    x(h),resnorm;
end

%plot(x)

f=2:0.01:5;
nf=length(f);

for i=1:nf
    for h=1:nq
        for j=1:N
            ysq(i,h,j)=f(i)*(1+(1-q(h))*x(1,h)*(j-1)*dt)^(1/(1-q(h)));
        end
    end
end

for i=1:nf
    for h=1:nq
        for j=1:N
            E(i,h)=mse(ysq(i,h,j),y(j));
        end
    end
end

mE=min(min(E));

[k1 k2] = find(E==mE);


for j=1:N+2
   ypsq(j)=f(k1)*(1+(1-q(k2))*x(1,k2)*(j-1)*dt)^(1/(1-q(k2)));   
end

plot(0:1:N-1,y,0:1:N+1,ypsq,'-or')


%%Datos para tabla
real=cases(length(y)+1,:);
prediccion=ypsq(:,length(y)+1);
e=(prediccion-real)/real;
qvalue=q(k2);
fvalue=f(k2);
datos=[real prediccion e mE qvalue fvalue];
dataTable=cat(1,dataTable,datos);
end




% % % Error=(y(N)-ypqs(N))/y(N);
% % z=[ypsq(N+1),a(k1),q(k2),x(1,k2)]%,Error]
% plot(0:1:N-1,y,0:1:N,ypsq,'-or')