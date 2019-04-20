clc
% ----------  Khoi tao cac gia tri ban dau 
M = 5; %channel
N = 4; %user secondary
K = 4; %user primary
dmax = N-1;  %the maximum transmission power of secondary users
dmin = 1;   %the minimun transmission power of secondary users
DPR =2; 
pop_size =20;
dim = M*N;
bound = 30;
Cmax = 5;
% ------------ Khoi tao ma tran float ngau nhien -------------
A = rand(pop_size,dim)*2*bound -bound;
fprintf('Ma tran A float\n\n')
%disp(A)
% -------------Chuyen ve ma tran nhi phan --------------------

fprintf('Ma tran A binary\n\n')
Ab = feval('convertFloattoBin',A,dim,pop_size);
%disp(Ab)
% -------------Chuyen ve ma tran nhi phan --------------------
a = 1;
b = 10;
for ku = 1:M            % Primary User
%      xk(ku) = a + (b-a)*rand(1,1) ;     % Hoanh do cua k
%      yk(ku) = a + (b-a)*rand(1,1);       % Tung do cua k
    xk = csvread('xk.txt');
    yk = csvread('yk.txt'); 
    ck(ku) = randi([1 M]);              % Channel 

end

for un = 1:N            % Seconday User
    xu(un) = a + (b-a)*rand(1,1);       % Hoanh do cua n
    yu(un) = a + (b-a)*rand(1,1);       % Tung do cua n
    sun = csvread('sun.txt');
end

for n=1:N % n dong
    for m=1:M
        benhat=max(sun);       
        flag=0;
        for k=1:K % k cot
            Dist(n,k)=0;
            if ck(k)== m
                flag=1;
                Dist(n,k) = sqrt((xu(n)-xk(k))^2 + (yu(n)-yk(k))^2)-DPR;             
                if (benhat > Dist(n,k))
                    benhat = Dist(n,k);
                end
                 %disp(benhat)
            end      
        end
        if flag
            DSE(n,m)=min(dmax,benhat);
        else
            DSE(n,m)=min(dmax,0);
        end
        
%--------------------- Dieu kien DSEnm -----------------------  
        if DSE(n,m) > dmin
            B(n,m)=DSE(n,m).^2;
            L(n,m)= 1;
        else
            B(n,m)=0;
            L(n,m)=0;
        end
    end   
end

for n=1:N-1
    for i=n+1:N
        for m=1:M            
            DISTni(n,i)= sqrt((xu(n)-xu(i))^2  + (yu(n)-yu(i))^2);
            if  DSE(n,m)+DSE(i,m) >= DISTni(n,i)
                C(n,i,m)= 1;
                C(i,n,m)=C(n,i,m);
            else
                C(n,i,m)= 0;
                C(i,n,m)=C(n,i,m);                
            end
        end
    end
end
    %%%%%%%%%%%%%--- DUA MA TRAN QUA 3 RANG BUOC ----%%%%%%%%%%%%%
% AL = feval('matxl',Ab,L,dim,pop_size);  
% AC = feval('matxc',AL,C,M,N,dim,pop_size);
% ACmax = feval('matxcmax',AC,Cmax,M,N,dim,pop_size);
% if true
%     format = '%g';
%     fileL = fopen('M_L.txt','r');
%     AL = fscanf(fileL,format);
%     fileB = fscanf('M_B.txt','r');
%     AC = fscanf(fileB,format);
%     fileC = fscanf('M_C.txt','r');
%     fscanf(fileC,format);
%     fclose(fileL);
%     fclose(fileB);
%     fclose(fileC);
% end
L = csvread('M_L.txt'); 
B = csvread('M_B.txt');
ObjC = csvread('M_C.txt');
C = permute(reshape(ObjC.', 20,20,[]),[2 1 3]);
AL = feval('matxl',Ab,L,dim,pop_size);
AC = feval('matxc',AL,C,M,N,dim,pop_size);
ACmax = feval('matxcmax',x,Cmax,M,N,dim,pop_size);
%----------------Dua ma tran ve lai float
fprintf('Ma tran sau khi xu ly rang buoc\n\n');
matA = feval('convertBintoFloat',A,ACmax,dim,pop_size);
%disp(matA)

    %%%%%%%%%%%%%--- VORTEX SEARCH ----%%%%%%%%%%%%%
upperlimit = max(max(matA));
lowerlimit = min(min(matA));
maxitr = 10000;
%itr = maxitr;
u = 0.5 * (upperlimit + lowerlimit) * ones(1,dim);
x = 0.1;
ginv = (1/x)*gammaincinv(x,1);
r = ginv * ((upperlimit - lowerlimit) / 2);

fsbest = inf;
% for count = 1:30
%     fprintf('So lan chay thu %d sau 10000 vong lap:\n',count);
for itr = 1:maxitr
    C = r *r* matA;
    Cs = feval('Gauss',C,u,dim);
    %Cs = bsxfun(@plus, C, u); 
    Cs(Cs < lowerlimit) =  rand*(upperlimit - lowerlimit) + lowerlimit;
    Cs(Cs > upperlimit) =  rand*(upperlimit - lowerlimit) + lowerlimit; 
    f1 = feval('MSUMR',Cs,M,B,pop_size);
    fsmin = max(f1) ;
    MinFitInd = find(f1 == fsmin);
    if numel(MinFitInd) > 1
        MinFitInd = MinFitInd(1); 
    end
    itrBest = Cs(MinFitInd,1:dim);
    if fsmin < fsbest         
        sbest = itrBest; 
        fsbest = fsmin;
    end    
   fprintf('Iter = %d Fsbest = %g\n',itr, fsbest); 
    u = sbest;
    a = maxitr -itr/maxitr;
    ginv = (1/x)*gammaincinv(x,a)  ;
    r = ginv * ((upperlimit - lowerlimit) / 2);
end       
% end

function y = convertFloattoBin(x,dim,pop_size)
    for i=1:pop_size
        for j=1:dim
            s=1/(1+exp(-x(i,j))); %S1 transfer function
            if rand < s
                y(i,j)=0;
            else
                y(i,j)=1;
            end
        end
    end
end

% Function Convert Binary to Float

function k = convertBintoFloat(x,y,dim,pop_size)
    for i = 1:pop_size
        for j = 1:dim
            if y(i,j) == 1
                y(i,j) = x(i,j);
            else
                y(i,j) = 0;
            end
        end
        k = y;
    end
end

% Function Constraint with L

function xl = matxl(xl,L,dim,pop_size)
    l=reshape(L.',1,dim); %chuyen mtran L(0,1) sang vector
for i=1:pop_size
    for t=1:dim
        if l(t)==0
            xl(i,t)=0;
        end
    end
end
end

% Function Constraint with C  

function xc=matxc(xc,C,M,N,dim,pop_size)
for i=1:pop_size
    mat=vec2mat(xc(i,:),M);
    for n=1:N
        for k=n
            for m=1:M
                if (C(n,k,m)== 1 && mat(n,m) == 1 && mat(k,m) == 1)
                    v=rand();
                    if v<0.5
                        mat(n,m)= 0;
                    else
                        mat(k,m)=0;
                    end
                end
            end
        end
    end 
    xc(i,:)=reshape(mat',1,dim);
end

end

% Function Constraint with C max 
function x=matxcmax(x,Cmax,M,N,dim,pop_size)
for y=1:pop_size
    mat=vec2mat(x(y,:),M);
    for n=1:N
        z=sum(mat,2);
        if z(n) > Cmax
            r=zeros(1,z(n)-Cmax);
            for i=1:z(n)-Cmax
                flag=1;
                while flag
                    a=randi(z(n),1,1);
                    flag=0;
                    for j=1:z(n)-Cmax
                        if r(j)== a
                            flag=1;
                        else
                            flag=0;
                        end
                    end                    
                end
                r(i)=a;
            end           
            for t=1:numel(r)
               if (mat(n,r(t))== 1)
                mat(n,r(t))=0;
               end
            end
        end    
    end 
    x(y,:)=reshape(mat',1,dim);
end
end

%Function MAX SUM REWARD for ham so fitness 
function f = MSUMR(y,M,B,pop_size)
% b=reshape(B',1,dim);
for i=1:pop_size
    mat = vec2mat(y(i,:),M);
    sumr=sum(mat.*B,2);%Sum Reward
    sr(i)=max(sumr); % max sum reward
end
f=1000-reshape(sr,pop_size,1);
end
function f = Gauss(C,u,dim)
   p1 = exp(-0.5.*(C-u).^2);
   p2 = sqrt((2*pi).^dim .* abs(C));
   f = p1./p2;
end




