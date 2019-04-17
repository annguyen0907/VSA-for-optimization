function f2b = converFloatToBinary(x,pop_size,dim)
    for i = 1:pop_size
        for j = 1:dim
            s = 1/(1+exp(-x(i,j)));
            if rand < s
                f2b(i,j) = 0;
            else
                f2b(i,j) = 1;
            end
        end
    end
end

