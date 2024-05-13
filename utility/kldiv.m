function KL = kldiv(pVect1,pVect2,type)
% % Input
% pVect1: n*1 or 1*n vector, probability of vector 1
% pVect2: n*1 or 1*n vector, probability of vector 2
% type (type of KL divergence):
% 'dir': direct KL
% 'js': Jensen-Shannon divergence
% 'sym': symmetric KL, Jeffreys divergence
% % Output
% KL: KL divergence value, the smaller, the better

% % Check probabilities sum to 1:
if (abs(sum(pVect1) - 1) > .000001) || (abs(sum(pVect2) - 1) > .000001)
    error('Probablities don''t sum to 1.')
end
% % calculate KL
% avoid zeros
pVect1(pVect1==0) = 1e-6;
pVect2(pVect2==0) = 1e-6;
switch type
    case 'dir'
        KL = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));

    case 'js'
        logQvect = log2((pVect2+pVect1)/2);
        KL = .5 * (sum(pVect1.*(log2(pVect1)-logQvect)) + ...
            sum(pVect2.*(log2(pVect2)-logQvect)));

    case 'sym'
        KL1 = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
        KL2 = sum(pVect2 .* (log2(pVect2)-log2(pVect1)));
        KL = (KL1+KL2)/2;
        
    otherwise
        error(['KLDiv type not correct'])
end








