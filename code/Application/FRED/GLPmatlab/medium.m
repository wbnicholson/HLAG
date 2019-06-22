% Application script Macro Medium-Large data set
clear;

% Load data
format long g;
load Ymedium_trainh1.dat;
load Ymedium_testh1.dat;

p = 4;
h = 1; %forecast horizon
check = size(Ymedium_trainh1);
rolsize = check(1)/76;
k = check(2);
MSFEs_GLP_medium_h1=zeros(76,k);

        
for j=1:76
        disp(sprintf('now running j= %d', j));
        Yraw=Ymedium_trainh1(((j-1)*rolsize+1):(rolsize*j), 1:end);
        Ytest=Ymedium_testh1(j,1:end);

        % GLP paper
        YrawGLP=Yraw;
        res = bvarGLP(YrawGLP, p, h); 
        YhatGLP = res.postmax.forecast(end,1:end);
        MSFEs_GLP_medium_h1(j,1:end)=(YhatGLP-Ytest).^2;

end


