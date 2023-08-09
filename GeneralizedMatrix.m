function [Final] = GeneralizedMatrix(bin_number, weight_shift)
w = weight_shift;
nbin = bin_number-1;
a=zeros(nbin,nbin);                   % Actual Matrix
for i=1:nbin
    for j=1:nbin
        if i<j
               a(i,j)=sqrt((j-w).^2-(i-1).^2)-sqrt((j-w).^2-i.^2);  % Coeff aij
        elseif i == j
            a(i,j) = sqrt((i-w).^2-(i-1).^2);       %   Coeff aii
        else
            a(i,j) = 0;
        end
    end
end
 Final = inv(a);       
end

