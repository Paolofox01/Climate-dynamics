function dydt = three_box_model(t, y, params)
    % y = [T1; T2; T3; S1; S2; S3]
    Tn = y(1); Ts = y(2);
    Sn = y(3); Ss = y(4);

    % Unpack parameters
    V = params.V;
    lambda = params.lambda;
    Tp = params.Tp;
    Te = params.Te;
    Phi_N = params.Phi_N;
    Phi_S = params.Phi_S;
    alpha = params.alpha;
    beta = params.beta;
    k = params.k;

    % Overturning strength (q)
    q = k* (alpha * (Tn - Ts) + beta * (Ss - Sn));

    % Preallocate derivative vector
    dydt = zeros(4,1);

    % Temperature equations
    if q >= 0
        dydt(1) = -(q/V)*Ts + lambda*(Te-Tp) - lambda*Tn - q*Tn;
        dydt(2) = -(q/V)*Ts + lambda*(Te - Tp) - q*(Ts - Tn) - lambda*Ts;
        
        dydt(3) = -(q/V)*Ss + (Phi_N + Phi_S)/V - q*Sn + Phi_N;
        dydt(4) = -(q/V)*Ss + (Phi_N + Phi_S)/V - q*(Ss - Sn) + Phi_S;
    else
        dydt(1) = -abs(q)/V*Tn + lambda*(Te - Tp) - abs(q)*(Tn - Ts) - lambda*Tn;
        dydt(2) = -abs(q)/V*Tn + lambda*(Te - Tp) - abs(q)*Ts - lambda*Ts;
        
        dydt(3) = -abs(q)/V*Sn + (Phi_N + Phi_S)/V - abs(q)*(Sn - Ss) + Phi_N;
        dydt(4) = -abs(q)/V*Sn + (Phi_N + Phi_S)/V - abs(q)*Ss + Phi_S;
    end
end


% Define parameters
params.V = 2;
params.lambda = 12.9 * 1e-10;
params.Tp = 0;
params.Te = 30;
params.Phi_N = 0.9 * 1e-10 * 1.5;
params.Phi_S = 0.9 * 1e-10;
params.alpha = 1.5 * 1e-4;
params.beta = 8 * 1e-4;
params.k = 1.5 * 1e-6;

y0 = [22; 26; 5; 20];

tspan = [0 100000000000];

[t_dim, y_dim] = ode45(@(t,y) three_box_model(t, y, params), tspan, y0);

figure;
subplot(2,1,1)
plot(t_dim, y_dim(:,1), t_dim, y_dim(:,2))
legend('T1','T2')
ylabel('Temperature')

subplot(2,1,2)
plot(t_dim, y_dim(:,3), t_dim, y_dim(:,4))
legend('S1','S2')
ylabel('Salinity')
xlabel('Time')

%%

function dydt = scaled_rooth_model(t, y, params)
    thetaN = y(1); thetaS = y(2);
    sN = y(3); sS = y(4);

    v = params.v;
    theta0 = params.theta0;
    phi1 = params.phi1;
    phi2 = params.phi2;

    xi = thetaN - thetaS + sS - sN;

    dydt = zeros(4, 1); 

    if xi >= 0
        dthetaN = -xi * (thetaS + v * thetaN) + theta0 - thetaN;
        dthetaS = xi * (v * thetaN - (1 + v) * thetaS) + theta0 - thetaS;
        
        dsN = -xi * (sS + v * sN) + phi1;
        dsS = xi * (v * sN - (1 + v) * sS) + phi2;
    else
        dthetaN = -xi * (v * thetaS - (1 + v) * thetaN) + theta0 - thetaN;
        dthetaS = xi * (thetaN + v * thetaS) + theta0 - thetaS;
        
        dsN = -xi * (v * sS - (1 + v) * sN) + phi1;
        dsS = xi * (sN + v * sS) + phi2;
    end

    dydt(1) = dthetaN;
    dydt(2) = dthetaS;
    dydt(3) = dsN;
    dydt(4) = dsS;
end

params.v = params.V;
params.theta0 = (params.k*params.alpha)/(params.lambda * params.V)*30;
params.phi1 = (params.beta * params.k * (params.v + 1) / (params.lambda^2 * params.v^2)) * params.Phi_N + ...
              (params.beta * params.k / (params.lambda^2 * params.v^2)) * params.Phi_S;

params.phi2 = (params.beta * params.k / (params.lambda^2 * params.v^2)) * params.Phi_N + ...
              (params.beta * params.k * (params.v + 1) / (params.lambda^2 * params.v^2)) * params.Phi_S;

y0 = [ (params.k*params.alpha)/(params.lambda * params.V)*22; 
    (params.k*params.alpha)/(params.lambda * params.V)*26; 
    (params.k*params.beta)/(params.lambda * params.V)*5; 
    (params.k*params.beta)/(params.lambda * params.V)*20];

tspan = [0 100000000000*params.lambda];

[t_nodim, y_nodim] = ode45(@(t,y) scaled_rooth_model(t, y, params), tspan, y0);

% Plot
figure;
subplot(2,1,1)
plot(t_nodim, y_nodim(:,1), t_nodim, y_nodim(:,2))
legend('\theta_N', '\theta_S')
ylabel('Temperature')

subplot(2,1,2)
plot(t_nodim, y_nodim(:,3), t_nodim, y_nodim(:,4))
legend('s_N', 's_S')
ylabel('Salinity')
xlabel('Time')


%%

y_transf = [ (params.k*params.alpha)/(params.lambda * params.V), ...
    (params.k*params.alpha)/(params.lambda * params.V),...
    (params.k*params.beta)/(params.lambda * params.V),...
    (params.k*params.beta)/(params.lambda * params.V)];

y_dim_to_nodim = y_dim .* y_transf;

t_dim_to_nodim = t_dim * params.lambda;



figure;


subplot(2,1,1)
hold on
plot(t_dim_to_nodim, y_dim_to_nodim(:,1), 'b-')
plot(t_dim_to_nodim, y_dim_to_nodim(:,2), 'black-')
plot(t_nodim, y_nodim(:,1), 'r--', 'LineWidth', 1.2);
plot(t_nodim, y_nodim(:,2), 'g--', 'LineWidth', 1.2);
legend('\theta_N', '\theta_S', 'scaled T_N', 'scaled T_S')
ylabel('Temperature')
hold off

subplot(2,1,2)
hold on
plot(t_dim_to_nodim, y_dim_to_nodim(:,3), 'b-')
plot(t_dim_to_nodim, y_dim_to_nodim(:,4), 'black-')
plot(t_nodim, y_nodim(:,3), 'r--', 'LineWidth', 1.2);
plot(t_nodim, y_nodim(:,4), 'g--', 'LineWidth', 1.2);
legend('s_N', 's_S', 'scaled S_N', 'scaled S_S')
ylabel('Salinity')
xlabel('Time')
hold off



%%

syms x y z w mu real

ode1 = mu + x * w;
ode2 = y - x;
ode3 = z - y;
ode4 = w - z;

dxdt = [ode1, ode2, ode3, ode4];

%%

function [times, y] = integrate_bifurcation(dxdt, vars, start_point, tspan)
    
    dvars = sym('d', size(vars));  
    assume(dvars, 'real');

    sys_bif = jacobian(dxdt, vars);

    eqs = sys_bif * dvars == 0;
    
    eqs = [eqs; dvars(1) == 1];

    sol = solve(eqs, dvars);

    sol_vec = [sol.d1 ; sol.d2; sol.d3; sol.d4; sol.d5];
    
    ode_rhs = matlabFunction(sol_vec, 'Vars', {vars});

    ode_bif = @(t, y) ode_rhs(y);

    [times, y] = ode45(ode_bif, tspan, start_point);
end



%%

sola0 = vpasolve(subs(dxdt, mu, -1) == [0, 0, 0, 0], [x,y,z,w]);

sola0x = double(sola0.x)
sola0y = double(sola0.y)
sola0z = double(sola0.z)
sola0w = double(sola0.w)

%%

[times1, y1] = integrate_bifurcation(dxdt, [mu; x; y; z; w], ...
    [-1; sola0x(1); sola0y(1); sola0z(1); sola0w(1)], [0 1]);

[times2, y2] = integrate_bifurcation(dxdt, [mu; x; y; z; w], ...
    [-1; sola0x(2); sola0y(2); sola0z(2); sola0w(2)], [0 1]);

[times3, y3] = integrate_bifurcation(dxdt, [mu; x; y; z; w], ...
    [-1; sola0x(1); sola0y(1); sola0z(1); sola0w(1)], [0 -3]);

[times4, y4] = integrate_bifurcation(dxdt, [mu; x; y; z; w], ...
    [-1; sola0x(2); sola0y(2); sola0z(2); sola0w(2)], [0 -3]);

%% exact solution

exact_sol1 = @(u) [sqrt(-u); sqrt(-u); sqrt(-u); sqrt(-u)]
exact_sol2 = @(u) [-sqrt(-u); -sqrt(-u); -sqrt(-u); -sqrt(-u)]



%%
figure;
hold on;

colors = lines(2); 

p1 = plot(y1(:, 1), y1(:, 2), 'Color', colors(1, :));
plot(y3(:, 1), y3(:, 2), 'Color', colors(1, :));

p2 = plot(y2(:, 1), y2(:, 2), 'Color', colors(2, :));
plot(y4(:, 1), y4(:, 2), 'Color', colors(2, :));

exact_sol1_plot = exact_sol1(-4:0.01:0);
exact_sol2_plot = exact_sol2(-4:0.01:0);

p3 = plot(-4:0.01:0, exact_sol1_plot(1,:), 'g--', 'LineWidth', 2);
p4 = plot(-4:0.01:0, exact_sol2_plot(1,:), 'y--', 'LineWidth', 2);


xlabel('mu')
ylabel('x')
title('Bifurcation Branches')

legend([p1, p2, p3, p4], ...
       {'Branch 1 (numerical)', ...
        'Branch 2 (numerical)', ...
        'Branch 1 (exact)', ...
        'Branch 2 (exact)'},"Location" ,'west');

grid on;


%% our model

syms thetaN thetaS sN sS real
syms phi1 real

phi_min = (params.beta * params.k / (params.lambda^2 * params.v^2)) * params.Phi_S;

xi = thetaN - thetaS + sS - sN;

phi2 = phi1/(params.v + 1) + params.Phi_S * params.beta * params.k * (params.v + 2)/(params.lambda.^2 * params.v * (params.v + 1));

dthetaN_pos = -xi*(thetaS + params.v*thetaN) + params.theta0 - thetaN;

% θS'
dthetaS_pos = xi*(params.v*thetaN - (1 + params.v)*thetaS) + params.theta0 - thetaS;

% sN'
dsN_pos = -xi*(sS + params.v*sN) + phi1;

% sS'
dsS_pos = xi*(params.v*sN - (1 + params.v)*sS) + phi2;

%%

dxdt_pos = [dthetaN_pos, dthetaS_pos, dsN_pos, dsS_pos]

sol_mu001_pos = vpasolve(subs(dxdt_pos, phi1, 0.1) == [0, 0, 0, 0], [thetaN, thetaS, sN, sS]);

solthN_pos = double(sol_mu001_pos.thetaN)
solthS_pos = double(sol_mu001_pos.thetaS)
solsN_pos = double(sol_mu001_pos.sN)
solsS_pos = double(sol_mu001_pos.sS)


%% now with xi negative 

dthetaN_neg = -xi*(params.v*thetaS - (1 + params.v)*thetaN) + params.theta0 - thetaN;

% θS'
dthetaS_neg = xi*(thetaN + params.v*thetaS) + params.theta0 - thetaS;

% sN'
dsN_neg = -xi*(params.v*sS - (1 + params.v)*sN) + phi1;

% sS'
dsS_neg = xi*(sN + params.v*sS) + phi2;

%%

dxdt_neg = [dthetaN_neg, dthetaS_neg, dsN_neg, dsS_neg]


sol_mu001_neg = vpasolve(subs(dxdt_neg, phi1, 0.1) == [0, 0, 0, 0], [thetaN, thetaS, sN, sS]);

solthN_neg = double(sol_mu001_neg.thetaN)
solthS_neg = double(sol_mu001_neg.thetaS)
solsN_neg = double(sol_mu001_neg.sN)
solsS_neg = double(sol_mu001_neg.sS)


%%

[times1,y1] = integrate_bifurcation(dxdt_pos, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_pos(1); solthS_pos(1); solsN_pos(1); solsS_pos(1)], 0:0.002:0.5);

[times2,y2] = integrate_bifurcation(dxdt_pos, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_pos(2); solthS_pos(2); solsN_pos(2); solsS_pos(2)], 0:0.002:0.5);

[times3,y3] = integrate_bifurcation(dxdt_pos, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_pos(1); solthS_pos(1); solsN_pos(1); solsS_pos(1)], 0:-0.002:-0.5);

[times4,y4] = integrate_bifurcation(dxdt_pos, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_pos(2); solthS_pos(2); solsN_pos(2); solsS_pos(2)], 0:-0.002:-0.5);

[times5,y5] = integrate_bifurcation(dxdt_neg, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_neg(1); solthS_neg(1); solsN_neg(1); solsS_neg(1)], 0:0.002:0.5);

[times6,y6] = integrate_bifurcation(dxdt_neg, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_neg(2); solthS_neg(2); solsN_neg(2); solsS_neg(2)], 0:0.002:0.5);

[times7,y7] = integrate_bifurcation(dxdt_neg, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_neg(1); solthS_neg(1); solsN_neg(1); solsS_neg(1)], 0:-0.002:-0.5);

[times8,y8] = integrate_bifurcation(dxdt_neg, [phi1; thetaN; thetaS; sN; sS], ...
    [0.1; solthN_neg(2); solthS_neg(2); solsN_neg(2); solsS_neg(2)], 0:-0.002:-0.5);

y1_q = [y1(:,1) , y1(:,2) - y1(:,3) - y1(:,4) + y1(:,5)];
y2_q = [y2(:,1) , y2(:,2) - y2(:,3) - y2(:,4) + y2(:,5)];
y3_q = [y3(:,1) , y3(:,2) - y3(:,3) - y3(:,4) + y3(:,5)];
y4_q = [y4(:,1) , y4(:,2) - y4(:,3) - y4(:,4) + y4(:,5)];
y5_q = [y5(:,1) , y5(:,2) - y5(:,3) - y5(:,4) + y5(:,5)];
y6_q = [y6(:,1) , y6(:,2) - y6(:,3) - y6(:,4) + y6(:,5)];
y7_q = [y7(:,1) , y7(:,2) - y7(:,3) - y7(:,4) + y7(:,5)];
y8_q = [y8(:,1) , y8(:,2) - y8(:,3) - y8(:,4) + y8(:,5)];

index_1 = y1_q(:,2) > 0 & y1_q(:,1) > phi_min;
index_2 = y2_q(:,2) > 0 & y2_q(:,1) > phi_min;
index_3 = y3_q(:,2) > 0 & y3_q(:,1) > phi_min;
index_4 = y4_q(:,2) > 0 & y4_q(:,1) > phi_min;
index_5 = y5_q(:,2) < 0 & y5_q(:,1) > phi_min;
index_6 = y6_q(:,2) < 0 & y6_q(:,1) > phi_min;
index_7 = y7_q(:,2) < 0 & y7_q(:,1) > phi_min;
index_8 = y8_q(:,2) < 0 & y8_q(:,1) > phi_min;



% For q > 0 and phi > phi_min
y1_q = y1_q(index_1, :);
y2_q = y2_q(index_2, :);
y3_q = y3_q(index_3, :);
y4_q = y4_q(index_4, :);

% For q < 0 and phi > phi_min
y5_q = y5_q(index_5, :);
y6_q = y6_q(index_6, :);
y7_q = y7_q(index_7, :);
y8_q = y8_q(index_8, :);



%%
figure;
hold on;

colors = lines(4); 

plot(y1_q(:, 1), y1_q(:, 2), 'Color', colors(1, :))
plot(y3_q(:, 1), y3_q(:, 2), 'Color', colors(1, :))

plot(y2_q(:, 1), y2_q(:, 2), 'Color', colors(2, :))
plot(y4_q(:, 1), y4_q(:, 2), 'Color', colors(2, :))

plot(y5_q(:, 1), y5_q(:, 2), 'Color', colors(3, :))
plot(y7_q(:, 1), y7_q(:, 2), 'Color', colors(3, :))

plot(y6_q(:, 1), y6_q(:, 2), 'Color', colors(4, :))
plot(y8_q(:, 1), y8_q(:, 2), 'Color', colors(4, :))

xlabel('\phi_1')
ylabel('\xi')
title('Bifurcation Branches')

xlim([phi_min 0.6])
grid on;

figure;
hold on;

colors = lines(4);  
plot(y1_q(:, 1), y1(index_1, 2), 'Color', colors(1, :))
plot(y3_q(:, 1), y3(index_3, 2), 'Color', colors(1, :))

plot(y2_q(:, 1), y2(index_2, 2), 'Color', colors(2, :))
plot(y4_q(:, 1), y4(index_4, 2), 'Color', colors(2, :))

plot(y5_q(:, 1), y5(index_5, 2), 'Color', colors(3, :))
plot(y7_q(:, 1), y7(index_7, 2), 'Color', colors(3, :))

plot(y6_q(:, 1), y6(index_6, 2), 'Color', colors(4, :))
plot(y8_q(:, 1), y8(index_8, 2), 'Color', colors(4, :))

xlabel('\phi_1')
ylabel('\Theta_N')
title('Bifurcation Branches')

xlim([phi_min 0.6])
grid on;

%%

figure;
hold on;

colors = lines(4); 
plot(y1_q(:, 1), y1_q(:, 2), 'Color', colors(1, :))
plot(y3_q(:, 1), y3_q(:, 2), 'Color', colors(1, :))

plot(y2_q(:, 1), y2_q(:, 2), 'Color', colors(2, :))
plot(y4_q(:, 1), y4_q(:, 2), 'Color', colors(2, :))

plot(y5_q(:, 1), y5_q(:, 2), 'Color', colors(3, :))
plot(y7_q(:, 1), y7_q(:, 2), 'Color', colors(3, :))

plot(y6_q(:, 1), y6_q(:, 2), 'Color', colors(4, :))
plot(y8_q(:, 1), y8_q(:, 2), 'Color', colors(4, :))

xlabel('\phi_1')
ylabel('\xi')
title('Bifurcation Branches')

xlim([phi_min 0.6])
grid on;


J_pos = jacobian(dxdt_pos, [thetaN, thetaS, sN, sS]);

J_func_pos = matlabFunction(J_pos, 'Vars', [phi1, thetaN, thetaS, sN, sS]);

J_neg = jacobian(dxdt_neg, [thetaN, thetaS, sN, sS]);

J_func_neg = matlabFunction(J_neg, 'Vars', [phi1, thetaN, thetaS, sN, sS]);

y_pos = {y1, y2, y3, y4};
y_neg = {y5, y6, y7, y8};


for cont = 1:length(y_pos)
    y_cont = y_pos{cont};
    for ind = 1:size(y_cont, 1)

        xi_cont = y_cont(ind, 2) - y_cont(ind, 3) - y_cont(ind, 4) + y_cont(ind, 5);
        if (xi_cont >= 0 && y_cont(ind,1) > phi_min)    
            eigen = eig(J_func_pos(y_cont(ind, 1), ...
                y_cont(ind, 2), y_cont(ind, 3), y_cont(ind, 4), y_cont(ind, 5)));
        
            if(all(eigen < 0 ))
                plot(y_cont(ind, 1), xi_cont, 'r.');
            elseif(any(real(eigen) > 0))
                plot(y_cont(ind, 1), xi_cont, 'b.');
            else
                plot(y_cont(ind, 1), xi_cont, 'g.');
            end
        end
    end
end

for cont = 1:length(y_neg)
    y_cont = y_neg{cont};
    for ind = 1:size(y_cont, 1)
        xi_cont = y_cont(ind, 2) - y_cont(ind, 3) - y_cont(ind, 4) + y_cont(ind, 5);
        if (xi_cont < 0 && y_cont(ind,1) > phi_min)    
            eigen = eig(J_func_neg(y_cont(ind, 1), ...
                y_cont(ind, 2), y_cont(ind, 3), y_cont(ind, 4), y_cont(ind, 5)));
        
            if(all(eigen < 0 ))
                plot(y_cont(ind, 1), xi_cont, 'r.');
            elseif(any(real(eigen) > 0))
                plot(y_cont(ind, 1), xi_cont, 'b.');
            else
                plot(y_cont(ind, 1), xi_cont, 'g.');
            end
        end
    end
end

%%

n1=35;

n2=33;

eig(J_func_pos(y1(n1, 1), ...
                y1(n1, 2), y1(n1, 3), y1(n1, 4), y1(n1, 5)))

eig(J_func_pos(y1(n2, 1), ...
                y1(n2, 2), y1(n2, 3), y1(n2, 4), y1(n2, 5)))

%%

n1=35;

n2=36;

eig(J_func_neg(y7(n1, 1), ...
                y7(n1, 2), y7(n1, 3), y7(n1, 4), y7(n1, 5)))

eig(J_func_neg(y7(n2, 1), ...
                y7(n2, 2), y7(n2, 3), y7(n2, 4), y7(n2, 5)))

%% 

plot(y1_q(100,1), y1_q(100,2) + 0.01, 'x', 'Color', 'red');
text(y1_q(100,1) + 0.01, y1_q(100,2) + 0.01, 'y_1', 'Color', 'red', 'FontSize', 10);

plot(y1_q(10,1), y1_q(10,2) + 0.01, 'x', 'Color', 'red');
text(y1_q(10,1) + 0.01, y1_q(10,2) + 0.01, 'y_2', 'Color', 'red', 'FontSize', 10);


plot(y3_q(42,1), y3_q(42,2) + 0.01, 'x', 'Color', 'red');
text(y3_q(42,1) + 0.01, y3_q(42,2) + 0.01, 'y_3', 'Color', 'red', 'FontSize', 10);


plot(y5_q(100,1), y5_q(100,2) + 0.01, 'x', 'Color', 'red');
text(y5_q(100,1) + 0.01, y5_q(100,2) + 0.01, 'y_4', 'Color', 'red', 'FontSize', 10);


plot(y5_q(10,1), y5_q(10,2) + 0.01, 'x', 'Color', 'red');
text(y5_q(10,1) + 0.01, y5_q(10,2) + 0.01, 'y_5', 'Color', 'red', 'FontSize', 10);

plot(y7_q(42,1), y7_q(42,2) + 0.01, 'x', 'Color', 'red');
text(y7_q(42,1) + 0.01, y7_q(42,2) + 0.01, 'y_6', 'Color', 'red', 'FontSize', 10);

%%

y0 = [y1(100,2); y1(100,3); y1(100,4); y1(100,5)] + 0.01;

params.phi1 = y1(100,1);

params.phi2 = params.phi1/(params.v + 1) + params.Phi_S * params.beta * params.k * (params.v + 2)/(params.lambda.^2 * params.v * (params.v + 1));

% Solve
[t, y] = ode45(@(t,y) scaled_rooth_model(t, y, params), tspan , y0);

% Plot
figure;
subplot(2,3,1)
plot(t, (y(:,1) - y(:,2)- y(:,3) + y(:,4)))
ylabel('\xi')
yline(y1_q(100,2), '--r', 'eq. with \xi>0')
yline(y5_q(100,2), '--r', 'eq. with \xi<0')
title('starting point y_1')


xlabel('Time')

%%

y0 = [y1(10,2); y1(10,3); y1(10,4); y1(10,5)] + 0.01;

params.phi1 = y1(10,1);
params.phi2 = params.phi1/(params.v + 1) + params.Phi_S * params.beta * params.k * (params.v + 2)/(params.lambda.^2 * params.v * (params.v + 1));

% Solve
[t, y] = ode45(@(t,y) scaled_rooth_model(t, y, params), tspan, y0);

% Plot
subplot(2,3,2)
plot(t, (y(:,1) - y(:,2)- y(:,3) + y(:,4)))
ylabel('\xi')
yline(y1_q(10,2), '--r', 'eq. with \xi>0')
yline(y5_q(10,2), '--r', 'eq. with \xi<0')
title('starting point y_2')

xlabel('Time')


%%

y0 = [y3(42,2); y3(42,3); y3(42,4); y3(42,5)] + 0.01;

params.phi1 = y3(42,1);
params.phi2 = params.phi1/(params.v + 1) + params.Phi_S * params.beta * params.k * (params.v + 2)/(params.lambda.^2 * params.v * (params.v + 1));

% Solve
[t, y] = ode45(@(t,y) scaled_rooth_model(t, y, params), tspan, y0);

% Plot
subplot(2,3,3)
plot(t, (y(:,1) - y(:,2)- y(:,3) + y(:,4)))
ylabel('\xi')
yline(y3_q(42,2), '--r', 'eq. with \xi>0')
yline(y7_q(42,2), '--r', 'eq. with \xi<0')
title('starting point y_3')

xlabel('Time')

%%

y0 = [y5(100,2); y5(100,3); y5(100,4); y5(100,5)] + 0.01;

params.phi1 = y5(100,1);

params.phi2 = params.phi1/(params.v + 1) + params.Phi_S * params.beta * params.k * (params.v + 2)/(params.lambda.^2 * params.v * (params.v + 1));

% Solve
[t, y] = ode45(@(t,y) scaled_rooth_model(t, y, params), tspan , y0);

% Plot
subplot(2,3,4)
plot(t, (y(:,1) - y(:,2)- y(:,3) + y(:,4)))
ylabel('\xi')
yline(y1_q(100,2), '--r', 'eq. with \xi>0')
yline(y5_q(100,2), '--r', 'eq. with \xi<0')
title('starting point y_4')

xlabel('Time')

%%

y0 = [y5(10,2); y5(10,3); y5(10,4); y5(10,5)] + 0.01;

params.phi1 = y5(10,1);
params.phi2 = params.phi1/(params.v + 1) + params.Phi_S * params.beta * params.k * (params.v + 2)/(params.lambda.^2 * params.v * (params.v + 1));

% Solve
[t, y] = ode45(@(t,y) scaled_rooth_model(t, y, params), tspan, y0);

% Plot
% subplot(2,3,5)
plot(t, (y(:,1) - y(:,2)- y(:,3) + y(:,4)))
ylabel('\xi')
yline(y1_q(10,2), '--r', 'eq. with \xi>0')
yline(y5_q(10,2), '--r', 'eq. with \xi<0')
title('starting point y_5')

xlabel('Time')


%%

y0 = [y7(27,2); y7(27,3); y7(27,4); y7(27,5)] + 0.01;

params.phi1 = y7(27,1);
params.phi2 = params.phi1/(params.v + 1) + params.Phi_S * params.beta * params.k * (params.v + 2)/(params.lambda.^2 * params.v * (params.v + 1));

% Solve
[t, y] = ode45(@(t,y) scaled_rooth_model(t, y, params), tspan, y0);

% Plot
hold on
% subplot(2,3,6)
plot(t, (y(:,1) - y(:,2)- y(:,3) + y(:,4)))
ylabel('\xi')
yline(y3_q(27,2), '--r', 'eq. with \xi>0')
yline(y7_q(27,2), '--r', 'eq. with \xi<0')
title('starting point y_6')

xlabel('Time')

%%


syms omega real
lambda = 1i * omega;

find_eigen = det(J_pos - lambda * eye(4));

real_part = simplify(real(find_eigen));
imag_part = simplify(imag(find_eigen));
 
dxdt_hopf = [dthetaN_pos; dthetaS_pos; dsN_pos; dsS_pos; real_part; imag_part];
vars = [phi1, thetaN, thetaS, sN, sS, omega];

hopf_jac = jacobian(dxdt_hopf, vars); 

f_fun = matlabFunction(dxdt_hopf, 'Vars', {vars});
J_fun = matlabFunction(hopf_jac, 'Vars', {vars});

x0 = [y1(60, :), 1.0]; 

for i = 1:100
    F = f_fun(x0);     
    J = J_fun(x0);      
    x0 = x0 - (J \ F)';   
end

hopf_1 = x0 

mu_1 = x0(1)

%%


syms omega real
lambda = 1i * omega;

find_eigen = det(J_neg - lambda * eye(4));

real_part = simplify(real(find_eigen));
imag_part = simplify(imag(find_eigen));
 
dxdt_hopf = [dthetaN_neg; dthetaS_neg; dsN_neg; dsS_neg; real_part; imag_part];
vars = [phi1, thetaN, thetaS, sN, sS, omega]; 

hopf_jac = jacobian(dxdt_hopf, vars); 

f_fun = matlabFunction(dxdt_hopf, 'Vars', {vars});
J_fun = matlabFunction(hopf_jac, 'Vars', {vars});

x0 = [y7(35,:), 1.0]; 

for i = 1:100
    F = f_fun(x0);  
    J = J_fun(x0);      
    x0 = x0 - (J \ F)';  
end

hopf_2 = x0 

mu_2 = x0(1)

%%

phi_N_1 = (mu_1 - params.beta * params.k * params.Phi_S/(params.lambda*params.v)^2) ... 
    * (params.lambda*params.v)^2/(params.beta * params.k * (params.v + 1))

phi_N_2 = (mu_2 - params.beta * params.k * params.Phi_S/(params.lambda*params.v)^2) ...
    * (params.lambda*params.v)^2/(params.beta * params.k * (params.v + 1))


%%

y0 = [ (params.lambda * params.V)/(params.k*params.alpha)*hopf_1(2); 
    (params.lambda * params.V)/(params.k*params.alpha)*hopf_1(3); 
    (params.lambda * params.V)/(params.k*params.beta)*hopf_1(4); 
    (params.lambda * params.V)/(params.k*params.beta)*hopf_1(5)];

params.Phi_N = phi_N_1 + 0.015e-10;

[t, y] = ode45(@(t,y) three_box_model(t, y, params), [0 2000000000000] , y0);

figure;
plot(t/(365*24*60*60), params.k*(params.alpha*(y(:,1) - y(:,2)) - params.beta*(y(:,3) - y(:,4))))
ylabel('q')


xlabel('Years')


%%

y0 = [ (params.lambda * params.V)/(params.k*params.alpha)*hopf_2(2); 
    (params.lambda * params.V)/(params.k*params.alpha)*hopf_2(3); 
    (params.lambda * params.V)/(params.k*params.beta)*hopf_2(4); 
    (params.lambda * params.V)/(params.k*params.beta)*hopf_2(5)];

params.Phi_N = phi_N_2 - 0.01e-11;

[t, y] = ode45(@(t,y) three_box_model(t, y, params), [0 10000000000000] , y0);

% Plot
figure;
plot(t/(365*24*60*60), params.k*(params.alpha*(y(:,1) - y(:,2)) - params.beta*(y(:,3) - y(:,4))))
ylabel('q')


xlabel('Years')