%
% d0calc(system,h,k,l,a,b,c); 
% calculates d-spacing for a given lattice parameter/hkl/system - cubic-tetra-ortho-hexag
%

function F = d0calc(system,h,k,l,a,b,c)
system=lower(system);
if system=='cubic'; F=a./(h.^2+k.^2+l.^2).^0.5; end
if system=='tetra'; F=1/sqrt((h^2+k^2)/a^2+l^2/c^2); end
if system=='ortho'; F=1/sqrt(h^2/a^2+k^2/b^2+l^2/c^2); end
if system=='hexag'; F=1./sqrt(4/3/a^2*(h.^2+h.*k+k.^2)+l.^2./c^2); end
