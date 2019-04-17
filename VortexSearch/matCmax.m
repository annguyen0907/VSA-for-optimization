function xcmax = matCmax(x,Cmax,M,N,dim,pop_size)
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
                    if mat(n,r(t))== 1
                        mat(n,r(t))=0;
                    end
                end
            end    
        end 
        xcmax(y,:)=reshape(mat',1,dim);
    end
end


