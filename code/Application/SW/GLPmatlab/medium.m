% Application script Macro Medium-Large data set
clear;

% Load data
format long g;
load Ymediumtrainh1.dat;
load Ymediumtesth1.dat;

p = 4;
h = 1; %forecast horizon
check = size(Ymediumtrainh1);
rolsize = check(1)/61;
k = check(2);
MSFEs_GLP_medium_h1=zeros(61,k);

        
for j=1:61
        disp(sprintf('now running j= %d', j));
        Yraw=Ymediumtrainh1(((j-1)*rolsize+1):(rolsize*j), 1:end);
        Ytest=Ymediumtesth1(j,1:end);

        % GLP paper
        YrawGLP=Yraw;
        res = bvarGLP(YrawGLP, p, h); 
        YhatGLP = res.postmax.forecast(end,1:end);
        MSFEs_GLP_medium_h1(j,1:end)=(YhatGLP-Ytest).^2;

end
