function G = fcn_distance_dependent_threshold(A,dist,hemiid,frac)
% fcn_distance_dependent_threshold      generate group matrix
%
%   G = fcn_distance_dependent_threshold(A,dist,hemiid,frac) generates a
%       group-representative structural connectivity matrix by preserving
%       within-/between-hemisphere connection length distributions.
%
%   Inputs:
%               A,      [node x node x subject] structural connectivity
%                       matrices.
%               dist,   [node x node] distance matrix
%               hemiid, indicator matrix for left (1) and right (2)
%                       hemispheres
%               frac,   scales connection density (1 = connection density
%                       of average subject).
%
%   Outputs:
%               G,      group matrix (binary)
%
%   Richard Betzel, Indiana University, 2018
%
if ~exist('frac','var')
    frac = 1;
    
end
[n,~,nsub] = size(A); % number nodes (n) and subjects (nsub)
C = sum(A > 0,3);     % consistency
W = sum(A,3)./C;      % average weight
W(isnan(W)) = 0;      % remove nans
Grp = zeros(n,n,2);   % for storing inter/intra hemispheric connections (we do these separately)
for j = 1:2
    if j == 1         % make inter- or intra-hemispheric edge mask
        d = +(hemiid == 1)*(hemiid' == 2);
        d = d | d';
    else
        d = +(hemiid == 1)*(hemiid' == 1) | +(hemiid == 2)*(hemiid' == 2);
        d = d | d';
    end
    D = nonzeros(bsxfun(@times,(A > 0),dist.*triu(d))); % mask connections
    M = round(frac*length(D)/nsub);                     % fraction of connections to retain
    dist_hemi = dist.*d;                                
    [x,y] = ecdf(D);
    x = round(x.*M);    % x normally goes 0 to 1, we scale so it goes 0 to M
    G = zeros(n);       % temporary connectivity matrix
    for i = 1:M         % loop over edges
        ind = (x >= (i - 1)) & (x < i); % find all possible edges that fall within a distance bin
        if sum(ind)                     % total possible
            yy = y(ind);                % find range of distances
            mask = dist_hemi >= min(yy) & dist_hemi <= max(yy); % make mask
            [u,v] = find(triu(mask,1)); % find elements
            indx = (v - 1)*n + u;       % vectorize
            c = C(indx);                % get consistency
            w = W(indx);                % get weights
            zz = sum(c == max(c));      % is there a unique max?
            if zz == 1                  % if so, use the max
                [~,indmax] = max(c);
                G(indx(indmax)) = 1;
            else                        % if not, find all elements with max value and break tie with weight
                aa = find(c == max(c));
                ww = w(aa);
                [~,indsort] = sort(ww,'descend');
                G(indx(aa(indsort(1)))) = 1;
            end
        end
    end
    Grp(:,:,j) = G;
end
G = sum(Grp,3); G = G + G';