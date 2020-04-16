%GEPnewpop function
%aimed to creat new population

%input:function set and terminate set,gene info,population info
%output:newpop

function newpop=GEPnewpop(F,T,C,geneHead,geneTail,popSize,chromNum)

newpop=[];

FT=[F T];
Tlength=length(T);
FTlength=length(FT);
Clength=length(C);

for i=1:popSize
    chromBody=[];
    
    for j=1:chromNum
        geneBody=[];
        for k=1:geneHead
            geneBody=[geneBody FT(randperm(FTlength,1))];
        end
        
        for k=1:geneTail
            geneBody=[geneBody T(randperm(Tlength,1))];
        end
        
        for k=1:geneTail
            geneBody=[geneBody int2str(randperm(Clength,1))];
        end
        
        chromBody=[chromBody geneBody];
    end
    
    newpop=[newpop;chromBody];
end

end