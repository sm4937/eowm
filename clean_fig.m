function [fig,ax] = clean_fig()
%clean_fig Clean figures up to one template 
%   Make figure white, remove in-ticks, increase font size to acceptable
%   size

ax = gca; fig = gcf; 
fig.Color = 'w';
ax.TickLength = [0 0];
ax.FontSize = 14;
% legend('Location','Best')
% legend('Box','Off')


end

