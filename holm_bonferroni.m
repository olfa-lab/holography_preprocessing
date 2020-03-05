function sig = holm_bonferroni(pvals, thresh)
    sig = [];
    
    [~, isort] = sort(pvals);
    n = numel(pvals);
    psort_adj = cummax(min((n + 1 - (1:n) ).*pvals(isort),1));
    p_adj(isort) = psort_adj;
    signif(isort) = psort_adj <= thresh;
    adj_thresh(isort) = thresh ./ (n + 1 - (1:n) );
    
    sig.p = pvals;
    sig.alpha = thresh;
    sig.p_adj = p_adj;
    sig.signif = signif;
    sig.adj_thresh = adj_thresh;
    sig.rank_seq = isort;
    sig.first_failed = find(isort == find(psort_adj > thresh, 1, 'first'),1,'first');
    
    
    
end