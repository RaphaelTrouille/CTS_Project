function H = CM_entropy(dat)

% H = CM_entropy(dat)
%
% input:
%   dat:  n_sig x n_obs data matrix. dat should contain integers from 1 to 
%         Nbin, Nbin being the number of different bins.

M = max(dat,[],2);
Phi = dat(1,:);
for n = 2:size(dat,1)
    Phi = Phi + (dat(n,:)-1)*prod(M(1:n-1));
end

Phi = sort(Phi);
N = [0 find(diff(Phi)) length(Phi)];
N = diff(N);
p = N/length(Phi);

H = -p*log2(p');
