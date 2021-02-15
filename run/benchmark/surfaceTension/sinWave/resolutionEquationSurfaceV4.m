%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the ration A(t)/Ao and 2 comparison analitic / numeric
%%
%% author  DENIS Flavien // TAS // DLR BREMEN
%% version V3
%% update  26.08/2016 
%% input   Read text-file extract from the OpenFoam simulation
%% input   Please check the value of the "Parameters User" (1st section)
%% output  Plot the evolution of the ratio A(t)/A0
%% output  Plot the error Anum/Aana
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

%oldFolder = cd()

%%%%%%%%%%%%%%%%% Intput variables

rho1=1         %Density liquid 1
rho2=1        %Density liquid 2
lambda=0.003   %Wavelength
H0=3*10^-5;     %Initial height
nu=0.001       %Kinematic viscosity
sigma=1       % Surface tension

autoSave=false; % Auto-save of th figure in the folder of the text-files


%%%%%%%%%%%%%%%%%%%%%%%

% adresse= cellstr(char('results32.txt', 'results64.txt', 'results128.txt'));
% adresse=transpose(cellstr(uigetfile('*.txt','MultiSelect','on')))
[adresse,path]=(uigetfile('*.txt','MultiSelect','on'))

adresse=transpose(cellstr(adresse))

cd(path)

T=0:10^-7:0.0004;

%Constantes
waveNb=(2*pi())/lambda %wavenumber
omega0=abs(sqrt((sigma*waveNb^3)/(rho1+rho2)))   %fundamental pulsation
epsilon=(nu*waveNb^2)/omega0
beta=(rho1*rho2)/(rho1+rho2)^2
TPrime=omega0*T;


%roots

p=[1 -4*beta*(epsilon*omega0)^0.5 2*(1-6*beta)*epsilon*omega0 4*(1-3*beta)*(epsilon*omega0)^(3/2) (1-4*beta)*(epsilon*omega0)^2+omega0^2];
r=roots(p);
%plot([-100:0.01:100],polyval(p,[-100:0.01:100]))


Z(1)=(r(2)-r(1))*(r(3)-r(1))*(r(4)-r(1));
Z(2)=(r(1)-r(2))*(r(3)-r(2))*(r(4)-r(2));
Z(3)=(r(1)-r(3))*(r(2)-r(3))*(r(4)-r(3));
Z(4)=(r(1)-r(4))*(r(2)-r(4))*(r(3)-r(4));

fun1 = @(t) exp(-t.^2);
%fun2 = @(x) 1-(2/sqrt(pi())).*integral(fun1,0,x);
fun2 = @(x) 1-(2/sqrt(pi())).*quadgk(fun1,0,x);

plotStyle = {'m','c','r','g','k','m','c','r','g','k','m','c','r','g','k'};


%a(t)

index=1;

for t=T

tPrime=omega0*t;
A1=(4*(1-4*beta)*epsilon^2)/(8*(1-4*beta)*epsilon^2+1);
A1=A1.*erfc(sqrt(epsilon*t));

sum=0;
for i=1:4
    sum=sum+(r(i)*omega0^2)/(Z(i)*(r(i)^2-epsilon*omega0)).*exp(((r(i)^2-epsilon*omega0).*tPrime)./(omega0)).*fun2((r(i)*tPrime.^0.5)/(omega0^0.5));

end

result(index)=sum+A1;

index=index+1;

clear sum A1;
end



figure(1);
plot(TPrime,abs(result))

xlabel('Dimensionless Time')
ylabel('|H/Ho|')
ylim([0 1])

hold on


clear result
%%%%%%%Comparison


 for i=1:length(adresse)
    
    %dataNum=textread(char(adresse(i)), '%d %d %d %d %d');
  % dataNum = load -ascii char(adresse(i))
    dataNum=dlmread(char(adresse(i)));
    figure(1)
    plot(omega0.*dataNum(:,1),dataNum(:,5)./H0,plotStyle{i});
    hold on
    
    index=1;
    for t=transpose(dataNum(:,1))
        tPrime=omega0*t;
        A1=(4*(1-4*beta)*epsilon^2)/(8*(1-4*beta)*epsilon^2+1);
        A1=A1.*erfc(sqrt(epsilon*t));

        sum=0;
        for k=1:4
            sum=sum+(r(k)*omega0^2)/(Z(k)*(r(k)^2-epsilon*omega0)).*exp(((r(k)^2-epsilon*omega0).*tPrime)./(omega0)).*fun2((r(k)*tPrime.^0.5)/(omega0^0.5));

        end

        result(index)=sum+A1;

        index=index+1;

        clear sum A1;
    end
    
    
    error=(dataNum(:,5)./H0)-abs(transpose(result));
    
    
    figure(2);
    grid on
    plot(omega0.*dataNum(:,1),error,plotStyle{i});
    hold on
    
%     condition=dataNum(:,1)>(25/omega0);
%     % Remove rows
%     dataNum(condition,:)=[];
%     error(condition)=[];
    
    error=error.^2;
%     norme(i)= sqrt((trapz(dataNum(:,1),error))/max(dataNum(:,1)))
    norme(i)= sqrt((trapz(dataNum(:,1).*omega0,error))/(max(dataNum(:,1)*omega0)-min(dataNum(:,1)*omega0)))
    
%     sum(error^2)
    
    
    clear dataNum
 end
 
 %%%%%%%%%%%%%%%%%%%
 
 figure(2);
 xlabel('Dimensionless Time')
 ylabel('Error a(t)-a_a_n_a_l_i_t_i_c')
 legend(adresse);
 
 
 figure(1)
 %adresse=[cellstr(char('Analytic : Prosperetti')); adresse];
 %legend(adresse);
 
 
 if autoSave==true
    prompt = 'Name of the figure files ? ';
    nameFile = input(prompt,'s')

    figure(1)
    savefig(strcat(nameFile,'1'));
    figure(2)
    savefig(strcat(nameFile,'2'));

 end
 

%cd(oldFolder)


