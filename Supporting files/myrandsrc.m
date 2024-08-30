function X = myrandsrc(M, N, A)
% Draws a MxN matrix of random numbers given the possible values and
% probability distribution in A
% e.g. myrandsrc(1,1,[1 2 3; 0.1 0.4 1.0]);  
    while 1
        % the number of elements to chose from
        sz = size(A,2);
        % generate some random numbers
        r = rand(M*N,1)*ones(1,sz);    
        % determine the correct elements
        r = sum(A(2,:) < r,2)+1;    
        % select the correct elements
        X = A(1,r);    
        % Reshape into M x N matrix
        X = reshape(X,M,N);
        % If shape reasonabily large, try to approximate underlying distrubtion
        % as good as possible
        if M*N > 50   
            counts = histcounts(X,'BinMethod','integers');
            counts = counts(counts~=0);
            probs = cumsum(counts/(N*M));
            try
                if pdist2(probs,A(2,:),'minkowski') < 0.02
                    break;
                end
            end
        else
            break;
        end   
    end
end

