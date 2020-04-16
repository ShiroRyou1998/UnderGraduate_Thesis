%GEP programe main v0.9.9
%FOR TEST DATA
%this program is aimed to ge regression function via GEP algorithm
%attention:this program apply basic GEP and GEP-RNC const method
%main.m include 3 bodies, datainput model, processing model, conduct model
%decode, select, mutate, IS, cross, 5 parts covered

%update:
%1.body established, some bug corrected
%2.select,mutate,IS,RIS,3 cross have added
%3.plot function added

%created by shiro_ryou in 2020/4/14/10/28
%last edited by shiro_ryou in 2020/4/15/15/23

clc
clear
tic

%initialize the program
popSize=30;%even#,100 is recomand
iterationMax=500;%500
indexFinish=0.01;
countFinish=0;

%input gene setting info
geneHead=7; 
chromNum=3;%each chrom has 3 genes using '+' to connect
F=['/' '+' '-' '*' 'q' 's'];
Fnary=[2 2 2 2 1 1];
T=['a?'];%if variable changes,plz change GEPfitness.m
C=[3.14 2.718 1.414 1.732 2.236 ...
    2 3 5 7];%num=9, equal to genetail

geneTail=(max(Fnary)-1)*geneHead+1;%genecost = genetail

pmutate=0.035;pis=0.1;pris=0.1;
pcrosss=0.4;pcrossd=0.2;pcrossg=0.1;
isLength=3;

%input data
load('D:\FunctionSoftware\matlab2016\bin\Matlabprogram\graduate_thesis\2.GEPalgorithm\testData_1.mat')
sourceData=testData;
%the firt raw is y, while others are xi
%plz change setting info in the GEPfitness.m when variable number changes

%body
%initialize parameter
maxfitness=0;
bestfitness=0;
bestchrom='shiro_ryou';
bestindividual='shiro_ryou is the best';

%initialize plot function
maxfitPlot=zeros(1,iterationMax);
maxvarPlot=zeros(1,iterationMax);
popfitPlot=zeros(1,iterationMax);

%creat newpop
newpop=GEPnewpop(F,T,C,geneHead,geneTail,popSize,chromNum);

%do the loop
for i=1:iterationMax
    
    %caculate fitness
    [fitnessList,varList,maxfitness,maxMathexp,maxchrom]=...
        GEPfitness(newpop,geneHead,geneTail,chromNum,F,T,Fnary,C,sourceData);  

    %continue?
    if abs(bestfitness-maxfitness)<indexFinish
        countFinish=countFinish+1;
        if countFinish==30
            %plot info
            popfitPlot(i)=maxfitness;
            maxfitPlot(i)=bestfitness;
            maxvarPlot(i)=bestvariance;
            disp('GEP reach converge, so end loop ahead of schedule')
            break;
        end
    else
        countFinish=0;
    end
    
    popfitPlot(i)=maxfitness;%plot info
    
    %update the bestindividual
    if maxfitness>bestfitness
        bestfitness=maxfitness;
        bestvariance=1000/bestfitness-1;
        bestindividual=maxMathexp;
        bestchrom=maxchrom;
    end
    
    %update iretation plot
    maxfitPlot(i)=bestfitness;
    maxvarPlot(i)=bestvariance;
    
    %select individuals that join GEPopration
    newpop=GEPselect(fitnessList,newpop);
    
    %mutate
    newpop=GEPmutate(newpop,geneHead,geneTail,chromNum,F,T,C,pmutate);
    
    %IS & RIS
    newpop=GEPdis(newpop,geneHead,geneTail,chromNum,isLength,F,pis,pris);
    
    %cross (single double & gene)
    newpop=GEPtcross(newpop,chromNum,pcrosss,pcrossd,pcrossg);

end

%result

%plot
subplot(2,1,1);
plot(1:i,maxfitPlot(1:i),'r');
title('Subplot 1: bestfitness')
xlabel('iteration times')
ylabel('bestfitness')
grid on

subplot(2,1,2);
plot(1:i,maxvarPlot(1:i),'g');
title('Subplot 2: bestvariance')
xlabel('iteration times')
ylabel('bestvariance')
grid on
%{ 
plot2
figure(2)
plot(1:i,popfitPlot(1:i),'b',1:i,bestfitness*ones(1,i),'r');
title('Plot 1: maxfitness of each generation')
xlabel('genaration')
ylabel('fitness')
grid on
%}

toc