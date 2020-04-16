%GEPfitness function
%fitness function : MSE

%aimed to get fitness of a pop in 3 steps
%step1:decode the chrom
%step2:caculate the variance
%step3:caculate the fitness

%input:pop, geneinfo, chorminfo, F T C and their features, sourcedata
%output:fitness list,variance list,maxfiness and its mathexpression

%warning:if soursedata or variable numbers changed, this function need too!

function [fitnessList,varList,maxfitness,maxMathexp,maxchrom]=...
    GEPfitness(pop,geneHead,geneTail,chromNum,F,T,Fnary,C,sourceData)

%setting -- connection symbol
connectSym='+';
%sourcedata processing %if var changers,ND changing
[dataNum,~]=size(sourceData); %if u want to use RMSE et al.
P=sourceData(:,1);
A=sourceData(:,2);
%{
B=sourceData(:,3);
G=sourceData(:,4);%variable c has conflict with C, hence using G
D=sourceData(:,5);
E=sourceData(:,6);
%}

%body
[popSize,chromSize]=size(pop);
geneSize=chromSize/chromNum;

%initialize
fitnessList=zeros(popSize,1);
varList=ones(popSize,1);
maxfitness=-114514;
maxMathexp='an error occurred!';
maxchrom='shiro_ryou is best';

%do the loop

for i=1:popSize
    popTemp=pop(i,:);
    chromExp=[];
    
    %divide chrom into genes and transfer to mathexpression
    for j=1:chromNum
        subGene=popTemp((1+(j-1)*geneSize):j*geneSize);
        mathexp=GEPdecode(subGene,geneHead,geneTail,geneSize,F,Fnary,T,C);
        chromExp=[chromExp connectSym mathexp];
    end
    chromExp(1)=[];
    
    %calculate var and fitness
    voidChromExp=['0*(a)+' chromExp];%if var changers,ND changing
    regressFun=inline(vectorize(voidChromExp));
    regressP=regressFun(A);%if var changers,ND changing

    %variance=sum( (regressP-P).^2 );%MSE
    %variance=sum( abs(regressP-P) )/dataNum;%MAE
    variance=sqrt(sum( (regressP-P).^2 )/dataNum);%RMSE
    fitness=1000/(1+variance);
    
    if isnan(fitness)
       % disp('NaN occurred');
        fitness=0;
        variance=Inf;
    end
   
    fitnessList(i)=fitness;
    varList(i)=variance;
    
    if fitness>maxfitness
        maxfitness=fitness;
        maxMathexp=chromExp;
        maxchrom=popTemp;
    end
    
    
end

end