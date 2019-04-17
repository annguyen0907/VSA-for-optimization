clear all
close all
clc
%ackley: upperlimit = 32, lowerlimit = -32, dimension = 30
%khoi tao ma tran 5x5 [-32,32]
%----------------------------------------------
maxitr = 5*5;
itr = maxitr;
A = randi([-32,32],5);      disp(A)
upperlimit = max(max(A));   %fprintf('upperlimit A = %d\n',upperlimit);
lowerlimit = min(min(A));   %fprintf('lowerlimit A = %d\n',lowerlimit);
dim = 30;
u = 0.5 * (upperlimit + lowerlimit) * ones(1,dim);
%fprintf('Tam ban dau la:'); disp(u);
x = 0.1;
ginv = (1/x)*gammaincinv(x,1);
r = ginv*(upperlimit - lowerlimit)*0.5;
%fprintf('Ban kinh ban dau la %f\n',r);
t = 1;
fsbest = inf;
disp('----------------------------');
while (itr)
    C = r*randn(itr,dim);
    Cs = bsxfun(@plus,C,u);   
    Cs(Cs < lowerlimit) =  rand *(upperlimit - lowerlimit) + lowerlimit;
    Cs(Cs > upperlimit) =  rand *(upperlimit - lowerlimit) + lowerlimit;
    disp(Cs)
    ObjVal = feval('ackley',Cs);
    disp(ObjVal)
    fsmin = min(ObjVal); 
    MinFitInd = find(ObjVal == fsmin);
    if numel(MinFitInd) > 1
        MinFitInd = MinFitInd(1); 
    end
    itrBest = Cs(MinFitInd,1:dim);
    if fsmin < fsbest         
        sbest = itrBest; 
        fsbest = fsmin;
    end
%OUTPUT FSBEST
    fprintf('Iter=%d Fsbest=%g\n',t, fsbest); 
    u = sbest;
    t = t+1;
    itr = itr-1;   
    a = itr / maxitr;
    ginv = (1/x)*gammaincinv(x,a);  
    r = ginv * ((upperlimit - lowerlimit) / 2);
    %fprintf('ban kinh giam di con %f\n',r);
end
    % có tam luc sau va ban kinh luc sau => tim upper,lower / giai pt bac 2
    %fprintf('ban kinh toi uu nhat: %f\n',r);
    %fprintf('tam toi uu nhat:\n');  disp(u);
    
    % t dang thac mac cai nay la fsbest/sbest ?
    %co upperlimit', lowerlimit' => chay ra ma tran moi voi up,low nho hon
    % demo giai he phuong trinh tuyen tinh
        %x(1)/2 + x(2)/2 = u / ones(1,dim);
        %x(1)/2 - x(2)/2 = r;
     % x(1) la upperlimit moi x(2) la lowerlimit moi => Ma tran moi   
       % k = [lowerlimit; upperlimit];
       % F = [k(1) + k(2) - 2*u/ones(1,dim); k(1) - k(2) - 2*r];
       % [k,fval,exitflag] = fsolve(F,k);
       n = [1 1; 1 -1];
       m = [2*u/ones(1,dim); (2*r)/((1/x)*gammaincinv(x,1))];
       k = n\m;
       %disp(k);
       %disp(k(1));
       %disp(k(2));
      
       if (k(1)<0)
        up = feval('round',-k(1));
        low = feval('round',k(2));
       end
       
       if k(1)>0
        up = feval('round',k(1));
       low = feval('round',-k(2));
       end
       %end
       %if k(1)>0
       %up = feval('round',k(1));
       %low = feval('round',-k(1));
       %end
       
function [ObjVal] = ackley(Cs)
    dim = size(Cs,2);
    a = 20; b = 0.2; c = 2*pi;
    f1 = 0; 
    f2 = 0;
    for i = 1:dim
        f1 = f1 + Cs(:,i).^2;
        f2 = f2 + cos(c.*Cs(:,i));
    end
    % Ham tinh toi uu hoa fsbest theo ackley
    ObjVal = -a*exp(-b * sqrt(1/dim*f1)) - exp(1/dim*f2) + a + exp(1);
end








