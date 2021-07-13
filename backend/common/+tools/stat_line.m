

function s = stat_line(d1,d2,dtext)
% function [s] = tools.stat_line(x1, x2, discription_text)
% compute some summary statistics for a pair of populations. The mean and 
%  standard deviation of each population are displayed, along with a small 
%  battery of comparison statistics (ttest2, ranksum, kstest). The struct
%  contains additional statistics (median,sem, paired t-test). Calling
%  stat_line with just a struct s reproduces the display (which is only 
%  generated if nargout == 0).  
%
% Version 0.2 Calvin Eiber 16-Apr-2020 (added paired t-test)

if nargin > 1 && islogical(d2)  
  d1a = d1(d2);  
  d2 = d1(~d2);
  d1 = d1a;
end

if nargin > 1

  d1 = reshape(d1(~isnan(d1)), [], 1);
  d2 = reshape(d2(~isnan(d2)), [], 1); 
  
  pair_ = @(f) [f(d1); f(d2)];
  
  s = struct;
  if nargin > 2, s.info = dtext; end
  s.n = pair_(@numel);

  s.avg = pair_(@nanmean);
  s.std = pair_(@nanstd);
  s.median = pair_(@nanmedian);
 
  s.sem = s.std ./ sqrt(s.n); 

  [~,s.p_ttest] = ttest2(d1,d2);
  [s.p_wilcox] = ranksum(d1,d2);
  [~,s.p_kstest] = kstest2(d1,d2);
  
  if s.n(1) == s.n(2)
    
    [s.corr_r,s.corr_r_pval] = corr(d1,d2); 
    [s.corr_tau,s.corr_tau_pval] = corr(d1,d2,'Type','Kendall'); 
    
    
    [~,s.p_paired_ttest] = ttest(d1,d2);  
    s.p_paired_wilcox = signrank(d1,d2);
  end
else s = d1; 
end
txt = sprintf('(%0.2f ± %0.2f) vs (%0.2f ± %0.2f) : p (%0.4f | %0.4f | %0.4f) [TT/WX/KS]\n', ...
                s.avg(1),s.std(1),s.avg(2),s.std(2),s.p_ttest, s.p_wilcox,s.p_kstest);
if isfield(s,'info'), txt = sprintf('%s:\t%s',s.info,txt); end
if nargout == 0, fprintf('%s',txt), clear, end
if nargout == 1 && nargin == 1, s = txt; end

