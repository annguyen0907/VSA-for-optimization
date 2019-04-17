function c = AriseIndividuals(x,pop_size,dim,a,b)
    for i = 1:pop_size
        for j = 1:dim
            if x(i,j) == 0
                x(i,j) = a + (b-a) * rand();
            end
        end
        c = x;
    end
end

