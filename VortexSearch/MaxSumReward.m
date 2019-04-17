function f = MaxSumReward(y,M,B,pop_size)
    for i=1:pop_size
        mat = vec2mat(y(i,:),M);
        sumr=sum(mat.*B,2);         % Sum Reward
        sr(i)=max(sumr);            % Max sum reward
    end
    f=1000-reshape(sr,pop_size,1);
end
