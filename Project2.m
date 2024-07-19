clc 
clearvars

global pmax dab sigmaB alpha epsilon Pw hab hae
pmax = 5;
dab = 5;
sigmaB = 10^(-10/3);
alpha = 4;
epsilon = 0.1;
Pw = 0.7;
noiseab = (randn(1,1)+1i*randn(1,1))/sqrt(2);
noiseae = (randn(1,1)+1i*randn(1,1))/sqrt(2);
hab = abs(noiseab);
hae = abs(noiseae);

a  = linspace(0.00001,20,10);
b = zeros(size(a));

% do thi PLS EXhaustive search
for i = 1:length(a)
   sigmaE = 10^(-10/3);
   b(i)=PLS_OPT_ExhaustiveSearch(a(i),sigmaE);
end
plot(a,b,'*-','LineWidth', 2, 'MarkerSize', 10);
xlabel('dae[m]');
ylabel('Secrecy rate [bps/HZ]');
xlim([0,max(a)]);
ylim([0,max(b)]);
hold on;
%do thi PLS_2
for i = 1:length(a)
   sigmaE = 10^(-10/3);
   b(i)=PLS_OPT_LP(a(i),sigmaE);
end
plot(a,b,'x-','LineWidth', 2, 'MarkerSize', 10);
xlabel('dae[m]');
ylabel('Secrecy rate [bps/HZ]');
xticks(0:5:20)
legend('Exhaustive search','Use of total power');

%do thi convert
figure;
%1
for i = 1:length(a)
   sigmaE = 10^(-10/3);
   b(i)=Convert_Exhaustive_Search(a(i),sigmaE);
end
array1 = unique(b)
xticks(0:5:20)
plot(a,b,'x-','LineWidth', 2, 'MarkerSize', 10);
hold on
%2
for i = 1:length(a)
   sigmaE = 10^(-10/3);
   b(i)=Convert_Total_Power(a(i),sigmaE);
end
plot(a,b,'*-','LineWidth', 2, 'MarkerSize', 10);
xlabel('dae[m]');
ylabel('Convert rate [bps/HZ]');
array2 = unique(b)
xticks(0:5:20)
yticks([array1,array2])
legend('Exhaustive search','Use of total power')


% tong hop
figure;
%PLS 10
for i = 1:length(a)
   sigmaE = 10^(-10/3);
   b(i)=PLS_OPT_LP(a(i),sigmaE);
end
plot(a,b,'*-','LineWidth', 2, 'MarkerSize', 10);
hold on;

%PLS 30
for i = 1:length(a)
   sigmaE = 10^(-30/3);
   b(i)=PLS_OPT_LP(a(i),sigmaE);
end
plot(a,b,'x-','LineWidth', 2, 'MarkerSize', 10);
hold on;

%Convert 10
for i = 1:length(a)
   sigmaE = 10^(-10/3);
   b(i)=Convert_Exhaustive_Search(a(i),sigmaE);
end
plot(a,b,'o-','LineWidth', 2, 'MarkerSize', 10);
array1 = unique(b);
hold on;

%Convert 30
for i = 1:length(a)
   sigmaE = 10^(-30/3);
   b(i)=Convert_Total_Power(a(i),sigmaE);
end
plot(a,b,'+-','LineWidth', 2, 'MarkerSize', 10);
hold on;
xlabel('dae[m]');
ylabel('Rate [bps/HZ]');
array2 = unique(b)
xticks(0:5:20)
yticks(sort([array1,array2,0:0.2:2]))
legend('\sigma_e = -10dB, Secrecy rate','\sigma_e=-30dB, Secrecy rate','\sigma_e = -10dB, Convert rate','\sigma_e = -30dB, Convert rate')


%nhieu tang
figure;

a  = linspace(-30,-10,10);
%pls nhieu tang
for i = 1:length(a)
   sigmaE = 10^(a(i)/3);
   b(i)=PLS_OPT_LP(10,sigmaE);
end
plot(a,b,'*-','LineWidth', 2, 'MarkerSize', 10);
hold on;

%Convert nhieu tang
for i = 1:length(a)
   sigmaE = 10^(a(i)/3);
   b(i)=Convert_Total_Power(10,sigmaE);
end
plot(a,b,'+-','LineWidth', 2, 'MarkerSize', 10);
hold on;
xlabel('dae[m]');
ylabel('Rate [bps/HZ]');
legend('PLS dae = 10m','Convert Communication dae = 10m')


%PLS Exhaustive Search
function result = PLS_OPT_ExhaustiveSearch(dae,sigmaE)
    global pmax dab sigmaB alpha epsilon Pw hab hae
    maxval = 0;
    for sum_val = 0:0.01:5
        for pab = 0:0.01:sum_val
            pj = sum_val-pab;
            sinre = pab*hae^2/(dae^alpha*sigmaE+pj*hae^2);
            sinrb = pab*(hab^2)/(dab^alpha*sigmaB+pj*hab^2);
            Rsec = log10(1+sinrb)-log10(1+sinre);
            if Rsec>=0
                maxval = max(Rsec,maxval);
            end
        end
    end
    result = maxval;
end

%PLS 
function result = PLS_OPT_LP(dae,sigmaE)
global pmax dab sigmaB alpha epsilon Pw hab hae
    yb = pmax*hab^2/dab^alpha/sigmaB;
    ye = pmax*hae^2/dae^alpha/sigmaE;
    cvx_begin quiet
         variables x 
         maximize((1+yb)/(1+ye)*(x*(1-ye/yb)+ye/yb))
         subject to
         x>=1/(1+yb);
         x<=1;
    cvx_end
    if cvx_optval>=0
        result = log10(cvx_optval);
    else 
        result = 0;
    end
end

%Convert_Communication
function result = Convert_Exhaustive_Search(dae,sigmaE)
    global pmax dab sigmaB alpha epsilon Pw hab hae
    maxval = 0;
    for sum_val = 0:0.01:5
        for pab = 0:0.01:sum_val
           pj = sum_val-pab;
           expr1 = log(pj/(pj+pab))/pab;
           FA_MD = 1-exp(expr1*pj)+exp(expr1*(pj+pab));
      
           if FA_MD>=1-epsilon
                expr2 = pab*hab^2/(dab^alpha*sigmaB+pj*hab^2);
                maxval = max(maxval,Pw*log10(1+expr2));
           end
        end
    end
    result = maxval;
end

function result = Convert_Total_Power(dae,sigmaE)
global pmax dab sigmaB alpha epsilon Pw hab hae
    yb = pmax*hab^2/dab^alpha/sigmaB;
    ye = pmax*hae^2/dae^alpha/sigmaE;
     cvx_begin quiet
         variables x 
         maximize(x)
         subject to
            x>=0;
            x<=1;
           rel_entr(x*pmax,1)+rel_entr((1-x)*pmax,1)-pmax*log(pmax)-x*pmax*log(epsilon)<=0;
     cvx_end
     expr2 = yb/(1/x+(1/x-1)*yb);
     result =  Pw*log10(1+expr2)
end

