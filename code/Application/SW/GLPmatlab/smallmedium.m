% Application script Macro Medium-Large data set
clear;

% Load data
format long g;
load Ysmallmediumtrainh1.dat;
load Ysmallmediumtesth1.dat;

p = 4;
h = 1; %forecast horizon
check = size(Ysmallmediumtrainh1);
rolsize = check(1)/61;
k = check(2);
MSFEs_GLP_smallmedium_h1=zeros(61,k);

        
for j=1:61
                disp(sprintf('now running j= %d', j));
        Yraw=Ysmallmediumtrainh1(((j-1)*rolsize+1):(rolsize*j), 1:end);
        Ytest=Ysmallmediumtesth1(j,1:end);

        % GLP paper
        YrawGLP=Yraw;
        res = bvarGLP(YrawGLP, p, h); 
        YhatGLP = res.postmax.forecast(end,1:end);
        MSFEs_GLP_smallmedium_h1(j,1:end)=(YhatGLP-Ytest).^2;
end

