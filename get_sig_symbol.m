function [symbol] = get_sig_symbol(p)
%get_sig_symbol 
%  Get proper sign for whether test is significant or not, based on p

symbol = 'x';
if p < 0.0001
    symbol = '***';
elseif p < 0.001
    symbol = '**';
elseif p < 0.05
    symbol = '*';
elseif p < 0.1
    symbol = '^';
end

end

