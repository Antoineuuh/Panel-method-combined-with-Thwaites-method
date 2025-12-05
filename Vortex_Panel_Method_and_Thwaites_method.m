clc; clear all; close all; warning('off','all');

% =========================================================================
% ------------------------ Custom Functions -------------------------------
% =========================================================================
function [x, y] = projection_of_vector(X, Y, Phi)
    x = X*cos(Phi) + Y*sin(Phi);
    y = -X*sin(Phi) + Y*cos(Phi);
end

function [index_k_sep, theta_pow_2, K, S, F, H, delta] = compute_thwaite_param(X, Ue, nu)

    K_sep = -0.090;
    c = log(0.45/2)/log(-K_sep);

    dUe_dx = fnval(fnder(spapi(5, X, Ue)), X); % Derivate using cubic spline interpolation.

    % Calculation of theta power 2.
    theta_pow_2 = zeros(1, length(X));    
    S = 0.0;
    for i = 2:length(X)
        S = S + (Ue(i-1)^5 + Ue(i)^5) * (X(i) - X(i-1));  % Integral via rectangle method.
        theta_pow_2(i) = nu *  (0.45 * S + Ue(1)^6*theta_pow_2(1)/nu) / Ue(i)^6;
    end
    
    % Calculation of K.
    K = theta_pow_2.*dUe_dx./nu;
    [~, index_k_sep] = min(abs(K - K_sep)); % Calculation of the index where the B.L separates.
    
    % Calculation of S.
    S = (K - K_sep).^c;    
    
    % Calculation of F.
    F = 0.45 - 6*K;
    
    % Calculation of H.
    H = (-(0.5*F - S)./K) - 2;
    
    % Calculation of δ*.
    delta = sqrt(theta_pow_2).*H;
    
end

% =========================================================================
% ------------------------ Flow Parameters --------------------------------
% =========================================================================
rho = 1.225; % Air density [kg/m^3]
mu = 1.802e-5; % Viscosity dynamic [Pa.s]
nu = mu/rho; % Viscosity cinetic [m^2/s]
U_inf = 13; % Upstream velocity [m/s]
AoA_deg = [0 5 10]; % Angle of attack [deg]
Re = U_inf/nu;

% =========================================================================
% ---------------- Airfoil Geometry & Discretization ----------------------
% =========================================================================
c = 1; % Chord [m]
function y_star = y_star(x_star)
    y_star = -0.12/0.20*(0.2969*(sqrt(x_star)) - 0.1260*x_star - 0.3516*(x_star.^2) + 0.2843*(x_star.^3) - 0.1036*(x_star.^4));
end

N = 200; % Number of panels

theta = linspace(0, 0.95*pi, ceil((N + 1)/2)); % [rad] Starting from the lower contour trailing edge.
x_airfoil_intrados = (1 - cos(theta))/(1 - cos(0.95*pi)); % Discretization of the intrados part.
y_airfoil_intrados = -y_star(x_airfoil_intrados);
x_airfoil_extrados = flip(x_airfoil_intrados); % Discretization of the extrados part.
y_airfoil_extrados = y_star(x_airfoil_extrados);
x_airfoil = [x_airfoil_extrados x_airfoil_intrados(2:end)]; % Concatenate both intrados/extrados parts.
y_airfoil = [y_airfoil_extrados y_airfoil_intrados(2:end)]; % The extrados part starts from the 2nd element because it shares the same first element as the intrados part.

figure; set(gca, 'FontName', 'Pt Serif'); hold on;
plot(x_airfoil, y_airfoil, "ko", "Markersize", 3);
plot(x_airfoil, y_airfoil, "k-", "linewidth", 0.5);
grid off; axis equal; xlim([-0.1 c+0.1]); ylim([-0.1 0.1]);
xlabel({"$\frac{x}{c}$"}, "interpreter", "latex"); ylabel({"$\frac{y}{c}$"}, "interpreter", "latex");
legend("Number of panels : " + N, "Airfoil equation curve"); hold off;

% =========================================================================
% ------------------------ Vortex Panel Method ----------------------------
% =========================================================================
[slope, b, x_mid, y_mid] = deal(zeros(1, N)); % Placing control points at the middle and calculating slope and semi length "b" of each panel.

for i = 1:N
    slope(i) = atan2(y_airfoil(i) - y_airfoil(i+1), x_airfoil(i) - x_airfoil(i+1)); % Slope of the i-th panel regarding the principal frame base.
    b(i) = sqrt((x_airfoil(i+1) - x_airfoil(i))^2 + (y_airfoil(i+1) - y_airfoil(i))^2) / 2; % Calculation of b = L/2 of each panel.
    x_mid(i) = (x_airfoil(i) + x_airfoil(i+1)) / 2; % Control points coordinates.
    y_mid(i) = (y_airfoil(i) + y_airfoil(i+1)) / 2;
end

gamma = zeros(length(AoA_deg), N + 1);

for alpha = 1:length(AoA_deg) % For each angle alpha for which we want to compute gamma.
    A = zeros(N + 1); % We initialize the matrix A & vector B of our system of N + 1 equations.
    B = zeros(N + 1, 1);
    s = zeros(1, N + 1); % Parametrization of the curve with s ∈ [0,1].

    for i = 1:N % For each panel in the airfoil.
        for j = 1:N % For each other panel including himself.
            % Calculation of the vector X-X' in the general frame where X is on the i-th panel and X' on the j-th panel.
            dx = x_mid(i) - x_mid(j);
            dy = y_mid(i) - y_mid(j);
            Phi = slope(j) - slope(i);

            % Obtaining the X vector coordinates in the frame of the j-th panel.
            [x, y] = projection_of_vector(dx, dy, slope(j));

            log_plus_B = log((x - b(j))^2 + y^2);
            log_minus_B = log((x + b(j))^2 + y^2);

            arctan_plus_B = atan((x - b(j)) / y);
            arctan_minus_B = atan((x + b(j)) / y);

            v_special_plus_B = -(x - b(j)) + y * arctan_plus_B;
            v_special_minus_B = -(x + b(j)) + y * arctan_minus_B;
            
            u_Right = -1/(2*pi)*(y/(4*b(j))*(log_plus_B - log_minus_B)-(x/(2*b(j))+1/2)*(arctan_plus_B - arctan_minus_B));
            u_Left = -1/(2*pi)*(-y/(4*b(j))*(log_plus_B - log_minus_B)-(x/(-2*b(j))+1/2)*(arctan_plus_B - arctan_minus_B));

            v_Right = -1/(2*pi)*(1/(2*b(j))*(v_special_plus_B - v_special_minus_B)+1/2*(x/(2*b(j))+1/2)*(log_plus_B - log_minus_B));
            v_Left = -1/(2*pi)*(-1/(2*b(j))*(v_special_plus_B - v_special_minus_B)+1/2*(x/(-2*b(j))+1/2)*(log_plus_B - log_minus_B));

            V_Right = u_Right*sin(Phi)+ v_Right*cos(Phi);
            V_Left = u_Left*sin(Phi)+ v_Left*cos(Phi);

            % Filling the A matrix with system coefficients.
            A(i,j) = A(i,j) + V_Right;
            A(i,j+1) = A(i,j+1) + V_Left;
        end

        B(i) = U_inf*(sin(slope(i) - deg2rad(AoA_deg(alpha))));
        s(i+1) = s(i) + 2*b(i);
    end

    A(N + 1, 1) = 1;
    A(N + 1, N + 1) = 1;

    gamma(alpha, :) = A\B;
end

% =========================================================================
% ------------------------ Vortex sheet strength --------------------------
% =========================================================================
figure; set(gca, 'FontName', 'Pt Serif'); hold on; grid on;
color = ["#0D7D87" "#8CC5E3" "#F45F74"];

for alpha = 1:length(AoA_deg)
    plot(s, gamma(alpha,:)/U_inf, "-","Linewidth", 1, "Color", color(alpha));
end

xlabel({"$\frac{s}{c}$"}, "interpreter", "latex"); ylabel({"$\frac{\gamma}{U_\infty}$"}, "interpreter", "latex");
legend(arrayfun(@(x) sprintf('α = %d°', x), AoA_deg, 'UniformOutput', false));
hold off;

% =========================================================================
% ------------------------ Pressure coefficient ---------------------------
% =========================================================================
figure; set(gca, 'FontName', 'Pt Serif'); hold on; grid on;
for alpha = 1:length(AoA_deg)
    Cp = zeros(1, N);
    for i = 1:N
        Cp(i) = 1 - ((((gamma(alpha, i+1) + gamma(alpha, i)) / 2) / U_inf)^2);
    end
    plot(s(1:end-1), Cp, "-","Linewidth", 1, "Color", color(alpha));
end
xlabel({"$\frac{s}{c}$"}, "interpreter", "latex"); ylabel({"$C_p$"}, "interpreter", "latex");
legend(arrayfun(@(x) sprintf('α = %d°', x), AoA_deg, 'UniformOutput', false));
hold off;

% =========================================================================
% ------------------ Streamlines & Stagnation point -----------------------
% =========================================================================
for aoa = 1:length(AoA_deg)
    
    X_OBS = linspace(-1.5*c, 1.5*c, 600);
    Y_OBS = linspace(-1.5*c, 1.5*c, 600);
    [X_OBS, Y_OBS] = meshgrid(X_OBS, Y_OBS);

    U = U_inf*cosd(AoA_deg(aoa));
    V = U_inf*sind(AoA_deg(aoa));

    alpha = zeros(1, N);
    beta = alpha;

    for i = 1:N

        alpha(i) = (gamma(aoa, i) - gamma(aoa, i+1)) / (2*b(i));
        beta(i) = (gamma(aoa, i) + gamma(aoa, i+1)) / 2;

        dx = X_OBS - x_mid(i);
        dy = Y_OBS - y_mid(i);
        [x, y] = projection_of_vector(dx, dy, slope(i));

        log_plus_B = log((x - b(i)).^2 + y.^2);
        log_minus_B = log((x + b(i)).^2 + y.^2);
        arctan_plus_B = atan((x - b(i)) ./ y);
        arctan_minus_B = atan((x + b(i)) ./ y);

        s_term = (x - b(i)) - (x + b(i));

        u = -1/(2*pi).*(alpha(i)/2.*y.*(log_plus_B - log_minus_B) - (alpha(i).*x + beta(i)).*(arctan_plus_B - arctan_minus_B));
        v = -1/(2*pi).*(alpha(i).*(-s_term+y.*(arctan_plus_B - arctan_minus_B)) + 0.5.*(alpha(i).*x+beta(i)).*(log_plus_B - log_minus_B));

        U_global = u.*cos(slope(i)) - v.*sin(slope(i));
        V_global = u.*sin(slope(i)) + v.*cos(slope(i));

        U = U + U_global;
        V = V + V_global;
    end

    Streamline_figures(aoa) = figure; set(gca, 'FontName', 'Pt Serif'); hold on;
    starty = linspace(-c, c, 100);
    startx = -0.4*c*ones(1, length(starty));
    streamline(X_OBS, Y_OBS, U, V, startx, starty, 'Color', color(2));
    
    fill(x_airfoil, y_airfoil, "k"); axis equal;
    
    [~, stag_index] = min(abs(gamma(aoa, 2:N)));
    stag_point = plot(x_airfoil(stag_index + 1), y_airfoil(stag_index + 1),'o','Markersize',3, "Color", color(3));
    div_stream = streamline(X_OBS, Y_OBS, -U, -V, x_airfoil(stag_index + 1), y_airfoil(stag_index + 1), 'Color', color(3), 'LineWidth', 1.2);
    
    grid off; axis equal; xlim([-0.2 c+0.2]); ylim([-0.2 0.2]);
    xlabel({"$\frac{x}{c}$"}, "interpreter", "latex"); ylabel({"$\frac{y}{c}$"}, "interpreter", "latex");
    legend([stag_point, div_stream], "Stagnation point", "Dividing streamline");
    hold off;
end

% =========================================================================
% ------------------------ Lift coefficient -------------------------------
% =========================================================================
Cl = zeros(1, length(AoA_deg));
for alpha = 1:length(AoA_deg)
    Gamma = 0;
    for i = 1:N
        Gamma = Gamma + 0.5*(gamma(alpha, i) + gamma(alpha, i+1))*2*b(i);
    end
    l = -rho*U_inf*Gamma;
    Cl(alpha) = l/(0.5*rho*c*(U_inf^2));

end
figure; set(gca, "FontName", "Pt Serif"); xlim([0, 10]); ylim([0, 1.3]); hold on;
x2 = [0 2 4 6 8 10]; y2 = [0 0.21 0.44 0.66 0.88 1.12];
p1 = polyfit(AoA_deg, Cl, 1); p2 = polyfit(x2, y2, 1);
eq1 = sprintf("C_l = %.3fα + %.3f", p1(1), p1(2));
eq2 = sprintf("C_l = %.3fα  %.3f", p2(1), p2(2));
plot(AoA_deg, Cl, "Color", color(2), "DisplayName", "Panel method");
plot(AoA_deg, Cl, "x", "Markersize", 10, "Color", color(2), "DisplayName", eq1);
plot([0 2 4 6 8 10], [0 0.21 0.44 0.66 0.88 1.12], "Color", color(3), "DisplayName", "Thin airfoil theory");
plot([0 2 4 6 8 10], [0 0.21 0.44 0.66 0.88 1.12], "x", "Markersize", 10, "Color", color(3), "DisplayName", eq2);
grid on; legend(); xlabel("α (deg)"); ylabel("$C_l$", "interpreter", "latex"); hold off;

% =========================================================================
% ------------------------ Thwaites' Method -------------------------------
% =========================================================================
fig_Ue = figure; set(gca, 'FontName', 'Pt Serif'); grid on; legend("α = 0°", "α = 5°", "α = 10°"); xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\frac{U_e}{U_\infty}$", "interpreter", "latex");
fig_Theta = figure; set(gca, 'FontName', 'Pt Serif'); grid on; xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\sqrt{Re}\frac{\theta}{c}$", "interpreter", "latex");
fig_K = figure; set(gca, 'FontName', 'Pt Serif'); grid on; xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\frac{K}{c}$", "interpreter", "latex");
fig_S = figure; set(gca, 'FontName', 'Pt Serif'); grid on; xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\frac{S}{c}$", "interpreter", "latex");
fig_F = figure; set(gca, 'FontName', 'Pt Serif'); grid on; xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\frac{F}{c}$", "interpreter", "latex");
fig_H = figure; set(gca, 'FontName', 'Pt Serif'); grid on; xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\frac{H}{c}$", "interpreter", "latex");
fig_Delta = figure; set(gca, 'FontName', 'Pt Serif'); grid on; xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\sqrt{Re}\frac{\delta *}{c}$", "interpreter", "latex");
fig_Delta_rc = figure; set(gca, 'FontName', 'Pt Serif'); grid on; xlabel("$\frac{x}{c}$", "interpreter", "latex"); ylabel("$\frac{\delta *}{r_c}$", "interpreter", "latex");

% Rayon de courvature.
dx = gradient(x_airfoil);
dy = gradient(y_airfoil);
d2x = gradient(dx);
d2y = gradient(dy);
rc = 1 ./ (abs(dx .* d2y - dy .* d2x) ./ (dx.^2 + dy.^2).^(3/2));

% Normal vector to each panel.
normale_x = -dy ./ sqrt(dx.^2 + dy.^2); 
normale_y = dx ./ sqrt(dx.^2 + dy.^2);

for alpha = 1:length(AoA_deg)

    [~, stag_index] = min(abs(gamma(alpha, 2:N)));
    
    % =====================================================================
    % ---------------------------- Upper part -----------------------------
    % =====================================================================
    X_up = s(stag_index+1:end) - s(stag_index+1);
    Ue_up = -gamma(alpha, stag_index+1:end);

    [i_k_sep_up, theta_pow_2_up, K_up, S_up, F_up, H_up, delta_up] = compute_thwaite_param(X_up, Ue_up, nu);

    % =====================================================================
    % ---------------------------- Lower part -----------------------------
    % =====================================================================
    X_low = s(stag_index+1) - s(1:stag_index+1); X_low = X_low(end:-1:1); 
    Ue_low = gamma(alpha, 1:stag_index+1); Ue_low = Ue_low(end:-1:1);

    [i_k_sep_low, theta_pow_2_low, K_low, S_low, F_low, H_low, delta_low] = compute_thwaite_param(X_low, Ue_low, nu);
    
    % =====================================================================
    % -- Profile augmented by the boundary layer displacement ticknesses --
    % =====================================================================
    
    x_airfoil_up = x_airfoil(stag_index+1:end) + delta_up .* normale_x(stag_index+1:end);  % Augmenter le profil supérieur
    y_airfoil_up = y_airfoil(stag_index+1:end) + delta_up .* normale_y(stag_index+1:end);  % Augmenter le profil supérieur
    
    x_airfoil_low = flip(x_airfoil(1:stag_index+1) + delta_low .* normale_x(1:stag_index+1));  % Augmenter le profil supérieur
    y_airfoil_low = flip(y_airfoil(1:stag_index+1) + delta_low .* normale_y(1:stag_index+1));  % Augmenter le profil supérieur
    
    % Plot different figures.
    figg = figure;
    copyobj(get(Streamline_figures(alpha), "Children"), figg); hold on;
    plot(x_airfoil_up(1:i_k_sep_up), y_airfoil_up(1:i_k_sep_up), "g--", "LineWidth", 1.8, "DisplayName", "Upper boundary layer");  
    plot(x_airfoil_low(1:i_k_sep_low), y_airfoil_low(1:i_k_sep_low), "c--", "LineWidth", 1.8, "DisplayName", "Lower boundary layer");  
    plot([x_airfoil_up(i_k_sep_up) x_airfoil_low(i_k_sep_low)], [y_airfoil(stag_index+1+i_k_sep_up) y_airfoil(stag_index+1-i_k_sep_low)], "mx", "LineWidth", 1.5, "DisplayName", "Separation points")
    figure(figg);
    
    figure(fig_Ue); hold on; 
    plot(X_up, abs(gamma(alpha,stag_index+1:end))./U_inf, "Color", color(alpha));
    plot(X_low, abs(gamma(alpha,1:stag_index+1))./U_inf, "--", "Color", color(alpha));

    figure(fig_Theta); hold on;
    plot(X_up(1:i_k_sep_up), sqrt(Re)*theta_pow_2_up(1:i_k_sep_up)/c, "Color", color(alpha));
    plot(X_up(i_k_sep_up), sqrt(Re)*theta_pow_2_up(i_k_sep_up)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(X_low(1:i_k_sep_low), sqrt(Re)*theta_pow_2_low(1:i_k_sep_low)/c, "--", "Color", color(alpha));
    plot(X_low(i_k_sep_low), sqrt(Re)*theta_pow_2_low(i_k_sep_low)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    figure(fig_K); hold on;
    plot(X_up(1:i_k_sep_up), K_up(1:i_k_sep_up)/c, "Color", color(alpha));
    plot(X_up(i_k_sep_up), K_up(i_k_sep_up)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(X_low(1:i_k_sep_low), K_low(1:i_k_sep_low)/c, "--", "Color", color(alpha));
    plot(X_low(i_k_sep_low), K_low(i_k_sep_low)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    figure(fig_S); hold on;
    plot(X_up(1:i_k_sep_up), S_up(1:i_k_sep_up)/c, "Color", color(alpha));
    plot(X_up(i_k_sep_up), S_up(i_k_sep_up)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(X_low(1:i_k_sep_low), S_low(1:i_k_sep_low)/c, "--", "Color", color(alpha));
    plot(X_low(i_k_sep_low), S_low(i_k_sep_low)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    figure(fig_F); hold on;
    plot(X_up(1:i_k_sep_up), F_up(1:i_k_sep_up)/c, "Color", color(alpha));
    plot(X_up(i_k_sep_up), F_up(i_k_sep_up)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(X_low(1:i_k_sep_low), F_low(1:i_k_sep_low)/c, "--", "Color", color(alpha));
    plot(X_low(i_k_sep_low), F_low(i_k_sep_low)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    figure(fig_H); hold on;
    plot(X_up(1:i_k_sep_up), H_up(1:i_k_sep_up), "Color", color(alpha));
    plot(X_up(i_k_sep_up), H_up(i_k_sep_up)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(X_low(1:i_k_sep_low), H_low(1:i_k_sep_low), "--", "Color", color(alpha));
    plot(X_low(i_k_sep_low), H_low(i_k_sep_low)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    figure(fig_Delta); hold on;
    plot(X_up(1:i_k_sep_up), sqrt(Re)*delta_up(1:i_k_sep_up)/c, "Color", color(alpha));
    plot(X_up(i_k_sep_up), sqrt(Re)*delta_up(i_k_sep_up)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(X_low(1:i_k_sep_low), sqrt(Re)*delta_low(1:i_k_sep_low)/c, "--", "Color", color(alpha));
    plot(X_low(i_k_sep_low), sqrt(Re)*delta_low(i_k_sep_low)/c, 'kx', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    if AoA_deg(alpha) == 0
        delta_rc_up = delta_up./rc(stag_index+1:end);
        delta_rc_low = delta_low./rc(1:stag_index+1);
        leg = sprintf("Re = %.3f", Re);
        figure(fig_Delta_rc); legend(); hold on;  
        plot(X_up(1:i_k_sep_up), delta_rc_up(1:i_k_sep_up), "Color", color(alpha), "DisplayName", leg);
    end

end

plots = findobj(fig_Ue, "Type", "Line"); 
fig_Ue.Position = [100, 100, 1000, 300];
legend([plots(6), plots(4), plots(2)], {"α = 0°", "α = 5°", "α = 10°"});

plots = findobj(fig_Theta, "Type", "Line"); 
fig_Theta.Position = [100, 100, 1000, 300];
legend([plots(12), plots(8), plots(4), plots(1)], {"α = 0°", "α = 5°", "α = 10°", "Separation point"}, "Location", "southeast");

plots = findobj(fig_K, "Type", "Line"); 
fig_K.Position = [100, 100, 1000, 300];
legend([plots(12), plots(8), plots(4), plots(1)], {"α = 0°", "α = 5°", "α = 10°", "Separation point"}, "Location", "northeast");
    
plots = findobj(fig_S, "Type", "Line");
fig_S.Position = [100, 100, 1000, 300];
legend([plots(12), plots(8), plots(4), plots(1)], {"α = 0°", "α = 5°", "α = 10°", "Separation point"}, "Location", "northeast");

plots = findobj(fig_F, "Type", "Line"); 
fig_F.Position = [100, 100, 1000, 300];
legend([plots(12), plots(8), plots(4), plots(1)], {"α = 0°", "α = 5°", "α = 10°", "Separation point"}, "Location", "southeast");

plots = findobj(fig_H, "Type", "Line");
fig_H.Position = [100, 100, 1000, 300];
legend([plots(12), plots(8), plots(4), plots(1)], {"α = 0°", "α = 5°", "α = 10°", "Separation point"}, "Location", "southeast");

plots = findobj(fig_Delta, "Type", "Line"); 
fig_Delta.Position = [100, 100, 1000, 300];
legend([plots(12), plots(8), plots(4), plots(1)], {"α = 0°", "α = 5°", "α = 10°", "Separation point"}, "Location", "southeast");