%% DRUMMING ARM MODEL VISUALIZATION
% Owen Douglas Pearl | CMU MBL | SPRING 2022

% Initialize drawing area
figh = figure;
figh.Color = [.95 .95 .95];
figh.Position(1:2) = [0 0];
figh.Position(3:4) = [1200 820];
figh.OuterPosition(3:4) = [1210, 850];
hold on

% Setup  Figure
pos2 = get(figh, 'Position'); % gives the position of current plot
new_pos2 = pos2 + [-.009 0 0.035 0];
set(figh, 'Position', new_pos2); % set new position of current plot

drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} )  ;     
RADtoDEG = 180/pi;

% Choose vector of k's to draw which frames
for k = frames
    % PLOT MODEL STATE
    plot([0, xe_viz_state(k)], [0, ye_viz_state(k)], 'r', 'LineWidth',18)
    plot([xe_viz_state(k), xw_viz_state(k)], [ye_viz_state(k), yw_viz_state(k)], 'g', 'LineWidth',18)
    plot([xw_viz_state(k), xf_viz_state(k)], [yw_viz_state(k), yf_viz_state(k)], 'b', 'LineWidth',18)
    plot([xf_viz_state(k), xb_viz_state(k)], [yf_viz_state(k), yb_viz_state(k)], 'Color', [.2, .2, .2, .5], 'LineWidth',7)
    plot([xf_viz_state(k), xbutt_viz_state(k)], [yf_viz_state(k), ybutt_viz_state(k)], 'Color', [.2, .2, .2, .5], 'LineWidth',7)
    plot(xb_viz_state(k), yb_viz_state(k),  '.', 'Color', [.2, .2, .2, .5], 'MarkerSize', 48)
    
    % Plot Circles on All Joints
    plot(0, 0, 'k.', 'MarkerSize', 60)
    plot(xe_viz_state(k), ye_viz_state(k), 'k.', 'MarkerSize', 60)
    plot(xw_viz_state(k), yw_viz_state(k), 'k.', 'MarkerSize', 60)
    plot(xf_viz_state(k), yf_viz_state(k), 'k.', 'MarkerSize', 60)

    % Edit Plot
    xlabel('X (m)')
    ylabel('Y (m)')
    set(gca,'Box','off')
    ylim(lim_y)
    xlim(lim_x)
    ax = figh.CurrentAxes; ax.FontSize = 14; ax.LineWidth = 1;
    set(gcf,'color','w');
end

hold on
% Plot Trajectories of all joint centers
plot(xe_viz_state, ye_viz_state, 'r-.', 'LineWidth',2)
plot(xw_viz_state, yw_viz_state, 'g-.', 'LineWidth',2)
plot(xf_viz_state, yf_viz_state, 'b-.', 'LineWidth',2)
plot(xb_viz_state, yb_viz_state, '-.', 'Color', [.2, .2, .2, .5], 'LineWidth',2)

yticks([-.6 -.4 -.2 0])
xticks([0 0.2 0.4 0.6])
