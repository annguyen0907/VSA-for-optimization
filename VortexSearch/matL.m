function xl = matL(xl,L,dim,pop_size)
    l=reshape(L',1,dim);
    for i=1:pop_size
        for j=1:dim
            if l(j)==0
                xl(i,j)=0;
            end
        end
    end
end

