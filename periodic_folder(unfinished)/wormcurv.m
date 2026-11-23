function [r, d, d2] = wormcurv(t,varargin) 

a = 1;
b = 1;
if nargin > 1 && ~isempty(varargin{1})
    a = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    b = varargin{2};
end


r =  [a*cos(t(:).'); b*sin(t(:).') + 0.3*sin(1.4.*cos(t(:).'))];
d =  [-a*sin(t(:).'); b*cos(t(:).') + 0.3*(1.4)*cos(1.4.*cos(t(:).')).*(-sin(t(:).'))];
d2 = [-a*cos(t(:).');-b*sin(t(:).') + ...
    (-0.42)*(cos(1.4*(cos(t(:).')).*cos(t(:).')  ))-... 
    (0.42*1.4*(sin(1.4*cos(t(:).')).*(sin(t(:).').^2) ))];


end