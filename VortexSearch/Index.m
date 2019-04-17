%

% ------------- Intitial Parameters ------------------
M = 3;          % Number of channels
K = 3;          % Number of primary users
N = 3;          % Number of secondary users
dmax = 4;       % The maximum transmission power of secondary users
dmin = 1;       % The minimun transmission power of secondary users
DPR = 2;        % Digital Process Reporter
pop_size = 3;   % Number of individuals
dim = M*N;      % Dimensions
Cmax = 5;       % Radio Interface Limit
bound = 30; 

    %%% ---------------- Matrices processing ---------------------%%%
%   Khoi tao ma tran A:
A = rand(pop_size,dim) * 2 * bound - bound;
A_binary = feval('converFloatToBinary',A,pop_size,dim);

%   Initialization position for primary & seconday users
a = 1;
b = 10;
for ku = 1:M                            % Primary User
     xk(ku) = a + (b-a)*rand(1,1) ;     % x_k
     yk(ku) = a + (b-a)*rand(1,1);      % y_k
     ck(ku) = randi([1 M]);             % Channel 

end

for un = 1:N                            % Seconday User
    xu(un) = a + (b-a)*rand(1,1);       % x_n
    yu(un) = a + (b-a)*rand(1,1);       % y_n
end
%   Distance N , K
for n=1:N 
    for m=1:M
        benhat=max(xu);       
        flag=0;
        for k=1:K 
            Dist(n,k)=0;
            if ck(k)== m
                flag=1;
                Dist(n,k) = sqrt((xu(n)-xk(k))^2 + (yu(n)-yk(k))^2)-DPR;             
                if (benhat > Dist(n,k))
                    benhat = Dist(n,k);
                end
            end      
        end
        if flag
            DSE(n,m)=min(dmax,benhat);
        else
            DSE(n,m)=min(dmax,0);
        end
      
%--------------------- Condition: DSEnm -----------------------  
        if DSE(n,m) > dmin
            B(n,m)=DSE(n,m).^2;
            L(n,m)= 1;
        else
            B(n,m)=0;
            L(n,m)=0;
        end
    end   
end

% Initialization matrices:  L B C
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

%   Matrix Processing A with condition of L, C, Cmax
AL = feval('matL',A_binary,L,dim,pop_size);
AC = feval('matC',AL,C,M,N,dim,pop_size);
ACmax = feval('matCmax',AC,Cmax,M,N,dim,pop_size);

%   Convert binary after condition to float
for i = 1:pop_size
    for j = 1:dim
        if ACmax(i,j) == 1
            A_float(i,j) = A(i,j);
        else
            A_float(i,j) = 0;
        end
    end
end
disp(A_float)
     %%% ---------------- Vortex Search ---------------------%%%

% ------------- Intitial Parameters ------------------
upperlimit = 30;
lowerlimit = -30;
maxitr = 20;
itr = maxitr;                       % iteration
x = 0.1;                            % x = 0.1 for gammaincinv(x,a) function                                    
ginv = (1/x)*gammaincinv(x,1);      % initially a = 1 
fsbest = inf;                       % fitness of the global min
% ---- Processing :
u = 0.5 * (upperlimit + lowerlimit) *  ones(1,dim);    %  initial center of the circle 
r = 0.5 * (upperlimit - lowerlimit) * ginv;            % initial radius 
count =1;
while itr
    obj = feval('AriseIndividuals',A_float,pop_size,dim,lowerlimit,upperlimit);
    C = r * r * obj;
	Cs = feval('gaussian_distribution',C,u,dim);
    % Limit the variables
    Cs(Cs < lowerlimit) =  rand*(upperlimit - lowerlimit) + lowerlimit;
    Cs(Cs > upperlimit) =  rand*(upperlimit - lowerlimit) + lowerlimit;    
    plot(obj,Cs,'.')
    grid on
    title('Bell Curve')
    xlabel('Randomly produced numbers')
    ylabel('Gauss Distribution')
    % Evaluate the candidate solutions
    objval = feval('MaxSumReward',Cs,M,B,pop_size);
    
    fsmin = max(objval); % maximun reward fitness value;
    fitness = find(objval == fsmin);   % find the min. fitness index
    
    if numel(fitness) > 1
        fitness = fitness(1); % if more than one solution keep one of them
    end
    sbest = Cs(fitness,1:dim);
    
    if fsmin < fsbest
        fsbest = fsmin;
        s = sbest;
    end
    fprintf('itr = %d Fitness = %g \n',count,fsbest);
    u = s;
    count = count + 1;
    itr = itr-1;   
    a = itr / maxitr;
    ginv = (1/x)*gammaincinv(x,a);  
    r = ginv * ((upperlimit - lowerlimit) / 2);
end


















