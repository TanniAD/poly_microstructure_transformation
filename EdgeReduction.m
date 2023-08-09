function [k,d, d_astm] = EdgeReduction(B)
    s = slice(B,[],[],length(B)/2);  
    c = get(s,'CData');
%         [u,v] = count_unique(c); % 2D: u = uniques, v = number of uniques
    b = c;
    b(2:end-1,2:end-1) = 0;
    m = b(:);
    for k = 1:length(m)
        c(c==m(k)) = NaN;
    end
    k = c;
    [p,q] = count_unique(k);
    d = nthroot(4*q/pi,2);
    figure()
    surf(k,'edgecolor', 'none'); view(2)
    [g,h] = count_unique(m);       % Discarded grain count
  for j = 1:length(unique(m))
      e(j) = 0.5*nthroot(4*h(j)/pi,2); % Half of edge grain dia
  end
  d_astm = vertcat(d,e');
end