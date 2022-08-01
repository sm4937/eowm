% plot_bumps.m

% plot some number of bumps, each at a given location, size, height (z-axis
% scaled to 1). size of whole grid goes from -1 to 1 (X,Y)
%
% NOTE: all inputs should be equal size!!!! (TODO: check for that)
%
% Uses cos^7 function from Sprague reports - mostly because we know where
% it will reach 0 and it looks enough like a gaussian to make us happy
%
% example: plot_bumps([-0.5 0.3 0],[-0.5 0.5 0],[0.25 0.45 0.1],[.5 1.0 0.1])
%
% Tommy Sprague, 5/16/2018

function [this_fig,bump_img,bump_grid_x,bump_grid_y] = plot_bumps(bump_x,bump_y,bump_size,bump_height)

gridvals = linspace(-1,1,151);
[grid_x,grid_y] = meshgrid(gridvals,gridvals);

bump_img = zeros(size(grid_x));

COS_POW = 7;

FWHM_CONST = 0.3974; % scale factor from cos^7 to FWHM - multiply size (in FWHM) by 1/this

% NOTE: needs r already computed, but do that easily below; sig is in SIZE
% (convert below)
%bump_fcn = @(r,sig,amp) (r<=sig) .* amp * ( 0.5 * ( 1 + cos(pi*r/sig))).^COS_POW ;
bump_fcn = @(myr,mysig,myamp) (myr<=mysig)  .* myamp .* ((0.5*(1 + cos(myr*pi/mysig))).^COS_POW);

for bb = 1:length(bump_x)
    
    this_r = sqrt((grid_x-bump_x(bb)).^2+(grid_y-bump_y(bb)).^2);
    
    %this_bump =  bump_height(bb) * (0.5 * ( 1 + cos(pi*this_r/(bump_size(bb)*(1/FWHM_CONST))).^COS_POW ));
    this_bump = bump_fcn(this_r,bump_size(bb)*(1/FWHM_CONST),bump_height(bb));
    
    bump_img = bump_img+this_bump;
    clear this_bump;
    
    % get rid of points too far away
    %this_bump(this_r>=(bump_size*(1/FWHM_CONST))) = 0;
end
this_fig = figure;
thissurf = surf(gridvals,gridvals,bump_img,'MeshStyle','row','EdgeAlpha',0.2);

xlim([-1 1]);
ylim([-1 1]);
zlim([0 1]);
set(gca,'CLim',[0 1],'DataAspectRatioMode','manual','Projection','perspective');

view([0 30]);
axis off;

% to return if necessary
bump_grid_x = grid_x;
bump_grid_y = grid_y;

set(gcf,'Renderer','painters')

return