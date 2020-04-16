%GEPdecode function
%this function is to decode a K-expression into a math expresstion
%uesing DoublyLinkedList method so skip the processing of ET
%version 1.0.0
%created by shiro_ryou in 2020/4/12/19/42
%last edited by shiro_ryou in 2020/4/14/15/35

%update:
%1. add constant set C,using DC field
%2. correct bug who use (sj) wrongly in the unary part
%3. correct bug who doesnot permit T element in the head
%4. replace sqrt with abs for occuring the i

%input:genehead,function set F,terminate set T, constant set C
%      a K-expresstion in the form of 'string'
%output:a math expresstion in the form of 'string'

function mathexp=GEPdecode(gene,geneHead,geneTail,geneSize,F,Fnary,T,C)
%input

kexpression.element=gene;
kexpression.position=[];
kexpression.dcfield=gene(geneHead+geneTail+1:geneSize);

%body
i=1;j=1;

if strfind(F,kexpression.element(1))
    
    %initialize the loop
    listexp={};
    kexpression.position(1)=[1];
    listexp=[listexp kexpression.element(i)];
    
    %do the loop
    while (i<=geneHead)
        
        if i>length(kexpression.position)
            break;
        end
        
        pos=kexpression.position(i);
        si=kexpression.element(i);
        
        if strfind(T,si)
            i=i+1;
            continue;
        end
        
        if Fnary(strfind(F,si))==1 %F is unary
            j=j+1;
            sj=kexpression.element(j);
            
           % if strfind(F,sj)
                listlength=length(listexp);
                listexp={listexp{1:pos} '(' sj ')' listexp{pos+1:listlength}};
                kexpression.position=(kexpression.position>pos)*3+...
                    kexpression.position;
                kexpression.position(j)=pos+2;
           % else
           %{
                listlength=length(listexp);
                listexp={listexp{1:pos} sj listexp{pos+1:listlength}};
                kexpression.position=(kexpression.position>pos)*1+...
                    kexpression.position;
                kexpression.position(j)=pos+1;
           %}
           % end
        else
            %F is binary
            j=j+1;
            sj=kexpression.element(j);
            
            if strfind(F,sj)
                listlength=length(listexp);
                listexp={listexp{1:pos-1} '(' sj ')' listexp{pos:listlength}};
                kexpression.position=(kexpression.position>=pos)*3+...
                    kexpression.position;
                kexpression.position(j)=pos+1;
            else
                listlength=length(listexp);
                listexp={listexp{1:pos-1} sj listexp{pos:listlength}};
                kexpression.position=(kexpression.position>=pos)*1+...
                    kexpression.position;
                kexpression.position(j)=pos;
            end
            
            pos=kexpression.position(i);
            j=j+1;
            sj=kexpression.element(j);
            
            if strfind(F,sj)
                listlength=length(listexp);
                listexp={listexp{1:pos} '(' sj ')' listexp{pos+1:listlength}};
                kexpression.position=floor(kexpression.position>pos)*3+...
                    kexpression.position;
                kexpression.position(j)=pos+2;
            else
                listlength=length(listexp);
                listexp={listexp{1:pos} sj listexp{pos+1:listlength}};
                kexpression.position=floor(kexpression.position>pos)*1+...
                    kexpression.position;
                kexpression.position(j)=pos+1;
            end
            
        end
        
        i=i+1;
        
    end
    
    %trans listexp into mathexp
    
    listlength=length(listexp);
    mathexp=[];
    
    ck=1;%constant
    for k=1:listlength
        if listexp{k}=='?' %constant
            mathexp=[mathexp num2str(C(str2double(kexpression.dcfield(ck))))];
            ck=ck+1;
        else if listexp{k}=='q' %abs
                mathexp=[mathexp 'abs'];
            else if listexp{k}=='s' %sin
                    mathexp=[mathexp 'sin'];
                else
                    mathexp=[mathexp listexp{k}];
                end
            end
        end

    end
    
else if kexpression.element(1)=='?'%const
        mathexp=num2str(C(1));
    else
        mathexp=kexpression.element(1);
    end

end

end