function X = topo_embed(x, M, opt)
% x: input before transform/embedding
% M: matrix, adjacency/laplacian/other
% X: input after embedding
switch opt.embed
        case 'none'
            % % X * matrix
            X = x*M;
        case 'eig'
            % % X * eigenvector
            [eigenvectors, eigenvalues] = eig(M);
%             X = x*eigenvectors;
            % % [TODO]: adjustable trunc 
            if opt.trunc
                eigenvectors(:,1:13) = [];
                X = x*eigenvectors;
            else
                X = x*eigenvectors;
            end
        case 'eig_v2'
            % % X * eigenvector
            [eigenvectors, eigenvalues] = eig(M);
            X = x*(M*eigenvectors);
        case 'svd'
            [U, S, V] = svd(M);
%             trace(S(1:23,1:23)) / trace(S)
            % % [TODO]: adjustable trunc 
            if opt.trunc
                V(27:end,:) = [];
                X = x*V';
            else
                X = x*V';
            end
        otherwise
            error('method not found')
end
end




















