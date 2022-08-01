%% Make hypothetical probability distributions
% Sarah Master, first-year talk & paper content

condcolors = [190 0 110; 0 110 190]./255;

xs = -179:180;

hard_error = -1;
hard_precision = 25;

easy_error = 1;
easy_precision = 25;

figure; subplot(2,2,1)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trial')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trial')
xlabel('Representation error')
ax = gca; ax.FontSize = 12;
title('No changes based on trial type')

hard_error = 0;
hard_precision = 25;

easy_error = 70;
easy_precision = 25;

subplot(2,2,2)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trial')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trial')
xlabel('Representation error')
ax = gca; ax.FontSize = 12;
title('Hard trials more accurate')


hard_error = -1;
hard_precision = 25;

easy_error = 1;
easy_precision = 50;

subplot(2,2,3)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trial')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trial')
xlabel('Representation error')
ax = gca; ax.FontSize = 12;
legend('Location','Best')
title('Hard trials more precise')

hard_error = 0;
hard_precision = 25;

easy_error = 70;
easy_precision = 50;

subplot(2,2,4)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','Hard trial')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Easy trial')
xlabel('Representation error')
ax = gca; ax.FontSize = 12;
legend('Location','Best')
title('Hard trials more accurate & precise')
fig = gcf; fig.Color = 'w';


%% Make theory figure

plot_bumps([0.5],[0.5],[0.25],[0.5])
fig = gcf; fig.Color = 'w';
zlabel('Activity of neural population')
xlabel('Cortical surface')

plot_bumps([0.5],[0.5],[0.25],[1])
fig = gcf; fig.Color = 'w';
