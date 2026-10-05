f = @(x) -x.^2+2*x.^3-6*x.^4;

fplot(f, [-0.5, 0.5])
ylabel("h(x)")
xlabel("x")

%%

dx = @(x,y) x.^2 + y.^2;

du = @(x) dx(x, f(x))

fplot(du, [-0.5, 0.5])

ylabel("u'")
xlabel("u")

%%

[x1, x2] = meshgrid(-0.5:0.02:0.5, -0.5:0.02:0.5);

u = x1.^2 + x2.^2;  
v = - x1.^2 - x2; 

[startX, startY] = meshgrid(-0.5:0.1:0.5, -0.5:0.1:0.5);

quiver(x1, x2, u, v, "Color", 'black', 'LineWidth', 0.7, 'AutoScale','on')

verts = stream2(x1, x2, u, v, startX, startY);

hold on
h_streams = streamline(verts, "Color", 'blue');

% Use the first streamline for legend reference
if ~isempty(h_streams)
    h2 = h_streams(1);
else
    h2 = plot(nan, nan);  % fallback dummy handle if no streamlines
end

xlabel('x')
ylabel('y')

h3 = fplot(f, 'r');

xlim([-0.5 0.5])

ylim([-0.5 0.5])

legend([h2, h3], {'streamline', 'h(x)'}, "Location","southeast");

grid on
hold off


%%

% Main plot
figure; hold on;

% Plot dx/dt curve
h1 = fplot(du, 'LineWidth', 2);

% Sample points and directions for arrows
x_arrows = linspace(-0.5, 0.5, 20);
f_arrows = du(x_arrows);
f_arrows_norm = 0.1 * sign(f_arrows);  % shorter arrows

% Plot arrows (no legend assigned here)
for i = 1:length(x_arrows)
    quiver(x_arrows(i), 0, f_arrows_norm(i), 0, ...
        'MaxHeadSize', 0.5, 'AutoScale', 'off', 'Color', 'k');
end

% Plot fixed points (no legend assigned here)
x_fixed = [0,0];
plot(x_fixed, zeros(size(x_fixed)), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

% Dummy plots for legend
quiver(nan, nan, 1, 0, 'Color', 'k', 'MaxHeadSize', 0.5);  % direction arrow
h3 = plot(nan, nan, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);   % fixed point

% Dummy arrow just outside axis limits
h2 = quiver(-999, 0, 0.1, 0, 'Color', 'k', 'MaxHeadSize', 0.5, 'AutoScale', 'off');


% Add legend using dummy handles
legend([h1, h2, h3], {'dx/dt', 'Direction', 'Fixed Point'}, "Location","northeast");

xlabel('u'); ylabel("u'");

xlim([-0.5 0.5])

ylim([-0.5 0.5])

title('1D Phase Portrait');
grid on;
ylim padded;

hold off

