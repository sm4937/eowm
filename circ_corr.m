function rho = circ_corr(a, b)

if size(a,2) > size(a,1)
	a = a';
end

if size(b,2) > size(b,1)
	b = b';
end

if length(a)~=length(b)
  error('Input dimensions do not match.')
end


a_diff = a - a';
b_diff = b - b';

ind = ones(size(a_diff));
ind = tril(ind, -1);

a_diff = a_diff(ind==1);
b_diff = b_diff(ind==1);

num = sin(a_diff(:))'*sin(b_diff(:));
den = sqrt((sin(a_diff(:))'*sin(a_diff(:)))*(sin(b_diff(:))'*sin(b_diff(:))));

rho = num/den;