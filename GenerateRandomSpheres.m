function P = GenerateRandomSpheres(nWant, Width, Radius)
% INPUT:
%   nWant:  Number of spheres
%   Width:  Dimension of 3d box as [1 x 3] double vector
%   Radius: [nWant x 1] vector, radii of spheres
% OUTPUT:
%   P:      [nWant x 3] matrix, centers
P      = zeros(nWant, 3);
R2     = (2 * Radius(:)) .^ 2;  % Squared to avoid expensive SQRT
iLoop  = 1;                     % Security break to avoid infinite loop
index  = 1;
while index <= nWant && iLoop < 1e6
  newP = rand(1, 3) .* (Width - 2 * Radius(index)) + Radius(index);
  % Auto-expanding, 
  Dist2 = sum((P(1:index-1, :) - newP) .^ 2, 2);
  if all(Dist2 > R2(1:index-1) + R2(index))
    % Success: The new point does not touch existing sheres:
    P(index, :) = newP;
    index       = index + 1; 
  end
  iLoop = iLoop + 1;
end

% Error if too few values have been found:
if index < nWant
  error('Cannot find wanted number of points in %d iterations.', iLoop)
end
end