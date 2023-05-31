%% Make hypothetical probability distributions
% Sarah Master, first-year talk & paper content

condcolors = [190 0 110; 0 110 190]./255;

xs = -179:180;

hard_error = 30;
hard_precision = 30;

easy_error = 30;
easy_precision = 30;

figure; subplot(2,2,1)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','More effort')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Less effort')
plot([0 0],ylim,'k--')
xlabel('Representation error')
clean_fig();


hard_error = 30;
hard_precision = 30;

easy_error = 45;
easy_precision = 30;

subplot(2,2,2)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','More effort')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Less effort')
plot([0 0],ylim,'k--')
xlabel('Representation error')
clean_fig();


hard_error = 30;
hard_precision = 30;

easy_error = 30;
easy_precision = 45;

subplot(2,2,3)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','More effort')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Less effort')
plot([0 0],ylim,'k--')
xlabel('Representation error')
clean_fig(); 


hard_error = 30;
hard_precision = 30;

easy_error = 45;
easy_precision = 45;

subplot(2,2,4)
hard_dist = normpdf(xs,hard_error,hard_precision);
easy_dist = normpdf(xs,easy_error,easy_precision); 
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','More effort')
hold on
plot(xs,easy_dist,'LineWidth',1.5,'Color',condcolors(2,:),'DisplayName','Less effort')
plot([0 0],ylim,'k--')
xlabel('Representation error')
clean_fig(); 

%% Make one Gaussian and change it around

location = 180;

hard_error = 0+location;
hard_precision = 30;

easy_error = 15+location;
easy_precision = 45;

figure();
subplot(3,1,1)
hard_dist = normpdf(xs,hard_error,hard_precision);
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','More effort')
hold on
plot([hard_error hard_error], [0,max(hard_dist)], '--k','DisplayName','\mu')
xcoord = hard_precision;
plot([hard_error xcoord],[hard_dist(xs==xcoord) hard_dist(xs==xcoord)],'--k','DisplayName','\sigma')
xlabel('Representation error'); ylabel('Posterior probability')
yticks([]); %xticks([-180 0 180])
clean_fig(); 

subplot(3,1,2)
hard_dist = normpdf(xs,easy_error,hard_precision);
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','More effort')
hold on
plot([easy_error easy_error], [0,max(hard_dist)], '--k','DisplayName','\mu')
xcoord = easy_error+hard_precision;
plot([easy_error xcoord],[hard_dist(xs==xcoord) hard_dist(xs==xcoord)],'--k','DisplayName','\sigma')
xlabel('Representation error'); ylabel('Posterior probability')
yticks([]); %xticks([-180 0 180])
clean_fig(); 

subplot(3,1,3)
hard_dist = normpdf(xs,easy_error,easy_precision);
plot(xs,hard_dist,'LineWidth',1.5,'Color',condcolors(1,:),'DisplayName','More effort')
hold on
plot([easy_error easy_error], [0,max(hard_dist)], '--k','DisplayName','\mu')
xcoord = easy_error+easy_precision;
plot([easy_error xcoord],[hard_dist(xs==xcoord) hard_dist(xs==xcoord)],'--k','DisplayName','\sigma')
xlabel('Representation error'); ylabel('Posterior probability')
yticks([]); %xticks([-180 0 180])
clean_fig(); 

%% Make theory figure

plot_bumps([0.5],[0.5],[0.25],[0.5])
fig = gcf; fig.Color = 'w';
zlabel('Activity of neural population')
xlabel('Cortical surface')

plot_bumps([0.5],[0.5],[0.25],[1])
fig = gcf; fig.Color = 'w';
