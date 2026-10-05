syms t real;

A = sym([-1/10, -3,-3;24/5, -61/10,-51/10; -9/5, 6, 5])

[V, D] = eig(A)

%%

sol1 = (V(:,1))*sym(exp(D(1, 1)*t));

sol1(t) = symfun(sol1, t);

f = rewrite(V(:,2) * exp(imag(D(2,2))*1i*t), 'sincos');

f = simplify(f);

% Extract real and imaginary parts element-wise
sol2(t) = exp(real(D(2,2))*t) * simplify(real(f));
sol3(t) = exp(real(D(2,2))*t) * simplify(imag(f));

syms a b c;

sol(t,a,b,c) = a*sol1(t) + b*sol2(t) + c*sol3(t);


function dydt = ode1(t, y)
    A = [-1/10, -3,-3;24/5, -61/10,-51/10; -9/5, 6, 5];
   dydt = A * y;
end



[times,y] = ode45(@ode1,[0 1],[1; 0; -1]);


%%

eqns = sol(0,a,b,c) == [1; 0; -1];

% Solve for a, b, c
S = solve(eqns, [a, b, c]);

S_num = structfun(@double, S, 'UniformOutput', false);

ft_init = sol(t,S_num.a, S_num.b, S_num.c);

f1 = matlabFunction(ft_init(1), 'Vars', t); 
f2 = matlabFunction(ft_init(2), 'Vars', t); 
f3 = matlabFunction(ft_init(3), 'Vars', t); 

t_values = linspace(0, 1, 1000);

f1_values = f1(t_values);
f2_values = f2(t_values);
f3_values = f3(t_values);


%% 2-d f1,f2


subplot(1,3,1)

hold on

plot(f1_values, f2_values);

plot( y(:, 1), y(:,2), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("y")

xlabel("x")

hold off

subplot(1,3,2)

hold on

plot(f1_values, f3_values);

plot( y(:, 1), y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("z")

xlabel("x")

hold off

subplot(1,3,3)

hold on

plot(f2_values, f3_values);

plot( y(:, 2), y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("z")

xlabel("y")

hold off
%% 2-d t,f1


subplot(1,3,1)

hold on

plot(t_values, f1_values);

plot(times, y(:,1), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("x")

xlabel("t")

hold off

subplot(1,3,2)

hold on

plot(t_values, f2_values);

plot(times, y(:,2), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("y")

xlabel("t")

hold off

subplot(1,3,3)

hold on

plot(t_values, f3_values);

plot(times, y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("z")

xlabel("t")

hold off



%% 3d f1,f2,f3

hold on 

plot3(f1_values, f2_values, f3_values);

plot3(y(:,1), y(:,2), y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

zlabel("z")

ylabel("y")

xlabel("x")

hold off

%% new start point

[times,y] = ode45(@ode1,[0 1],[1; 1; 1]);


%%

eqns = sol(0,a,b,c) == [1; 1; 1];

% Solve for a, b, c
S = solve(eqns, [a, b, c])

S_num = structfun(@double, S, 'UniformOutput', false);

ft_init = sol(t,S_num.a, S_num.b, S_num.c);

f1 = matlabFunction(ft_init(1), 'Vars', t);
f2 = matlabFunction(ft_init(2), 'Vars', t); 
f3 = matlabFunction(ft_init(3), 'Vars', t);

t_values = linspace(0, 1, 1000); 

f1_values = f1(t_values);
f2_values = f2(t_values);
f3_values = f3(t_values);


%% 2-d f1,f2

figure;
subplot(1,3,1)

hold on

plot(f1_values, f2_values);

plot( y(:, 1), y(:,2), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("y")

xlabel("x")

hold off

subplot(1,3,2)

hold on

plot(f1_values, f3_values);

plot( y(:, 1), y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("z")

xlabel("x")

hold off

subplot(1,3,3)

hold on

plot(f2_values, f3_values);

plot( y(:, 2), y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("z")

xlabel("y")

hold off
%% 2-d t,f1

figure;
subplot(1,3,1)

hold on

plot(t_values, f1_values);

plot(times, y(:,1), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("x")

xlabel("t")

hold off

subplot(1,3,2)

hold on

plot(t_values, f2_values);

plot(times, y(:,2), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("y")

xlabel("t")

hold off

subplot(1,3,3)

hold on

plot(t_values, f3_values);

plot(times, y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

ylabel("z")

xlabel("t")

hold off



%% 3d f1,f2,f3

hold on 

plot3(f1_values, f2_values, f3_values);

plot3(y(:,1), y(:,2), y(:,3), "r", 'LineWidth', 2, 'Color', "r", "LineStyle","--")

zlabel("z")

ylabel("y")

xlabel("x")

hold off
%% part 2

syms a u v w u_prime v_prime w_prime real;

% Define the derivatives symbolically
u_prime = 2*(a - 6) + (a - 12)*u - 15*v - 3*w - 3*u^2 - 12*u*v - 3*u^2*v;
v_prime = 4*a + 5 + (2*a + 1)*u + (a - 6)*v + (a - 1)*w - u^2 + 4*u*w + u^2*w;
w_prime = -4*a - 2*(a + 1)*u + 9*v - w + 3*w^2 - 6*v*w + 3*v*w^2;

% Optionally display them
disp(u_prime)
disp(v_prime)
disp(w_prime)

dxdt = [u_prime, v_prime, w_prime];

ftot = matlabFunction(dxdt, 'Vars', [a, u, v, w]);

simplify(ftot(a, -2, v, -v))

ftot(a, -2, -1, 1)



%%

yline(-2 )
xlabel("a")
ylabel("u")

%%

function [times, y] = run_bifurcation_ode(dxdt, start_point, tspan)

    syms a u v w real
    syms da du dv dw real

    u_prime = dxdt(1);
    v_prime = dxdt(2);
    w_prime = dxdt(3);


    x = [da; du; dv; dw];

    Df = [diff(u_prime, a), diff(u_prime, u), diff(u_prime, v), diff(u_prime, w)];
    Dg = [diff(v_prime, a), diff(v_prime, u), diff(v_prime, v), diff(v_prime, w)];
    Dh = [diff(w_prime, a), diff(w_prime, u), diff(w_prime, v), diff(w_prime, w)];

    Dtot = [Df; Dg; Dh];

    sys_bif = [diff(dxdt, a); diff(dxdt, u); diff(dxdt, v); diff(dxdt, w)].';

    eqs = sys_bif * x == 0;

    sol_struct = vpasolve([eqs; da == 1], [da, du, dv, dw]);

    sol = [sol_struct.da; sol_struct.du; sol_struct.dv; sol_struct.dw];

    ode_rhs = matlabFunction(sol, 'Vars', [a; u; v; w]);
    ode_bif = @(t, y) ode_rhs(y(1), y(2), y(3), y(4));

    % Run ODE solver
    [times, y] = ode45(ode_bif, tspan, start_point);
end

%%

[times, y] = run_bifurcation_ode(dxdt, [1; -2; -1; 1], [0 1]);


hold on
yline(-2)

plot(y(:, 1), y(:, 2), 'r--', 'LineWidth', 1.3)
xlabel("a")
ylabel("u")

ylim([-4 0])

hold off




%%


sola0 = vpasolve(subs(dxdt, a, 0) == [0, 0, 0], [u,v,w]);

sola0u = double(sola0.u);
sola0v = double(sola0.v);
sola0w = double(sola0.w);


[times1,y1] = run_bifurcation_ode(dxdt, [0; sola0u(1); sola0v(1); sola0w(1)], [0 100]);
[times2,y2] = run_bifurcation_ode(dxdt, [0; sola0u(2); sola0v(2); sola0w(2)], [0 100]);
[times3,y3] = run_bifurcation_ode(dxdt, [0; sola0u(3); sola0v(3); sola0w(3)], [0 100]);
[times4,y4] = run_bifurcation_ode(dxdt, [0; sola0u(1); sola0v(1); sola0w(1)], [0 -100]);
[times5,y5] = run_bifurcation_ode(dxdt, [0; sola0u(2); sola0v(2); sola0w(2)], [0 -100]);
[times6,y6] = run_bifurcation_ode(dxdt, [0; sola0u(3); sola0v(3); sola0w(3)], [0 -100]);

figure;
hold on;

colors = lines(3);  % three distinct colors

% Plot each trajectory with grouping colors
plot(y1(:, 1), y1(:, 2), 'Color', colors(1, :))
plot(y4(:, 1), y4(:, 2), 'Color', colors(1, :))

plot(y2(:, 1), y2(:, 2), 'Color', colors(2, :))
plot(y5(:, 1), y5(:, 2), 'Color', colors(2, :))

plot(y3(:, 1), y3(:, 2), 'Color', colors(3, :))
plot(y6(:, 1), y6(:, 2), 'Color', colors(3, :))

xlabel('a')
ylabel('u')
title('Bifurcation Branches')


grid on;

%%

J = jacobian(dxdt, [u; v; w]);

J_func = matlabFunction(J, 'Vars', [a, u, v, w]);


figure;
hold on;

colors = lines(3);  % three distinct colors

% Plot each trajectory with grouping colors
plot(y1(:, 1), y1(:, 2), 'Color', colors(1, :))
plot(y4(:, 1), y4(:, 2), 'Color', colors(1, :))

plot(y2(:, 1), y2(:, 2), 'Color', colors(2, :))
plot(y5(:, 1), y5(:, 2), 'Color', colors(2, :))

plot(y3(:, 1), y3(:, 2), 'Color', colors(3, :))
plot(y6(:, 1), y6(:, 2), 'Color', colors(3, :))


y_x = {y1, y2, y3, y4, y5, y6};

for cont = 1:length(y_x)
    y_cont = y_x{cont};
    for ind = 1:size(y_cont, 1)
    
        eigen = eig(J_func(y_cont(ind, 1), y_cont(ind, 2), ...
            y_cont(ind, 3), y_cont(ind, 4)));
    
        if(all(eigen < 0 ))
            plot(y_cont(ind, 1), y_cont(ind, 2), 'r.');
        elseif(any(real(eigen) > 0))
            plot(y_cont(ind, 1), y_cont(ind, 2), 'b.');
        else
            plot(y_cont(ind, 1), y_cont(ind, 2), 'g.');
        end
    end
end
xlabel('a')
ylabel('u')
title('Bifurcation Branches')
grid on;

%%

disp(eig(J_func(-1e-4,-2,-1,1)))
disp(eig(J_func(0,-2,-1,1)))
disp(eig(J_func(1e-4,-2,-1,1)))

%% plot for negative a

item =  ftot(y4(27, 1), u, v, w).';

ftot_negative_a = matlabFunction(item, 'Vars', [u, v, w]);

ode_negative_a = @(t, y) ftot_negative_a(y(1), y(2), y(3));

[time_eq1, y_eq1] = ode45(ode_negative_a, [0 10], [y4(27, 2)+1e-4; y4(27, 3)+1e-4;y4(27, 4)+1e-4]);

[time_eq2, y_eq2] = ode45(ode_negative_a, [0 10], [y4(27, 2)-1e-4; y4(27, 3)-1e-4;y4(27, 4)-1e-4]);

item =  ftot(y5(27, 1), u, v, w).';

ftot_negative_a = matlabFunction(item, 'Vars', [u, v, w]);

ode_negative_a = @(t, y) ftot_negative_a(y(1), y(2), y(3));

[time_eq3, y_eq3] = ode45(ode_negative_a, [0 10], [y5(27, 2)+1e-4; y5(27, 3)+1e-4;y5(27, 4)+1e-4]);

[time_eq4, y_eq4] = ode45(ode_negative_a, [0 10], [y5(27, 2)-1e-4; y5(27, 3)-1e-4;y5(27, 4)-1e-4]);

item =  ftot(y6(27, 1), u, v, w).';

ftot_negative_a = matlabFunction(item, 'Vars', [u, v, w]);

ode_negative_a = @(t, y) ftot_negative_a(y(1), y(2), y(3));

[time_eq5, y_eq5] = ode45(ode_negative_a, [0 10], [y6(27, 2)+1e-4; y6(27, 3)+1e-4;y6(27, 4)+1e-4]);

[time_eq6, y_eq6] = ode45(ode_negative_a, [0 10], [y6(27, 2)-1e-4; y6(27, 3)-1e-4;y6(27, 4)-1e-4]);


hold on
plot(time_eq1, y_eq1(:, 1))
plot(time_eq2, y_eq2(:, 1))
plot(time_eq3, y_eq3(:, 1))
plot(time_eq4, y_eq4(:, 1))

plot(time_eq5, y_eq5(:, 1))
plot(time_eq6, y_eq6(:, 1))

% Plot one labeled equilibrium point
plot(0, y5(27, 2), 'r.', 'DisplayName', 'Equilibrium point');

% Plot the others without legend entries
plot(0, y4(27, 2), 'r.', 'HandleVisibility', 'off');
plot(0, y6(27, 2), 'r.', 'HandleVisibility', 'off');

hold off
legend(...
    'u1 + \delta', ...
    'u1 - \delta', ...
    'u2 + \delta', ...
    'u2 - \delta', ...
    'u3 + \delta', ...
    'u3 - \delta', ...
    'Equilibrium point', ...
    'Location', 'best');
xlim([0,2])

%% plot for positive a

item =  ftot(y1(27, 1), u, v, w).';

ftot_negative_a = matlabFunction(item, 'Vars', [u, v, w]);

ode_negative_a = @(t, y) ftot_negative_a(y(1), y(2), y(3));

[time_eq1, y_eq1] = ode45(ode_negative_a, [0 10], [y1(27, 2)+1e-4; y1(27, 3)+1e-4;y1(27, 4)+1e-4]);

[time_eq2, y_eq2] = ode45(ode_negative_a, [0 10], [y1(27, 2)-1e-4; y1(27, 3)-1e-4;y1(27, 4)-1e-4]);

item =  ftot(y2(27, 1), u, v, w).';

ftot_negative_a = matlabFunction(item, 'Vars', [u, v, w]);

ode_negative_a = @(t, y) ftot_negative_a(y(1), y(2), y(3));

[time_eq3, y_eq3] = ode45(ode_negative_a, [0 10], [y2(27, 2)+1e-4; y2(27, 3)+1e-4;y2(27, 4)+1e-4]);

[time_eq4, y_eq4] = ode45(ode_negative_a, [0 10], [y2(27, 2)-1e-4; y2(27, 3)-1e-4;y2(27, 4)-1e-4]);

item =  ftot(y3(27, 1), u, v, w).';

ftot_negative_a = matlabFunction(item, 'Vars', [u, v, w]);

ode_negative_a = @(t, y) ftot_negative_a(y(1), y(2), y(3));

[time_eq5, y_eq5] = ode45(ode_negative_a, [0 10], [y3(27, 2)+1e-4; y3(27, 3)+1e-4;y3(27, 4)+1e-4]);

[time_eq6, y_eq6] = ode45(ode_negative_a, [0 10], [y3(27, 2)-1e-4; y3(27, 3)-1e-4;y3(27, 4)-1e-4]);


hold on
plot(time_eq1, y_eq1(:, 1))
plot(time_eq2, y_eq2(:, 1))
plot(time_eq3, y_eq3(:, 1))
plot(time_eq4, y_eq4(:, 1))

plot(time_eq5, y_eq5(:, 1))
plot(time_eq6, y_eq6(:, 1))

plot(0, y5(27, 2), 'r.', 'DisplayName', 'Equilibrium point');

plot(0, y4(27, 2), 'r.', 'HandleVisibility', 'off');
plot(0, y6(27, 2), 'r.', 'HandleVisibility', 'off');
hold off
legend(...
    'u1 + \delta', ...
    'u1 - \delta', ...
    'u2 + \delta', ...
    'u2 - \delta', ...
    'u3 + \delta', ...
    'u3 - \delta', ...
    'Equilibrium point', ...
    'Location', 'best');

%%
hold on
plot3(y_eq5(:, 1), y_eq5(:, 2),y_eq5(:, 3))
plot3(y3(27, 2), y3(27, 3), y3(27, 4), 'r.')
hold off

%% 

