%% ME-474 Numerical Flow Simulation — Assessment #1 (Fall 2025)
% Solution Q1→Q22 (focalisé Q17–Q20 avec QUICK borné en deferred correction)
% -----------------------------------------------------------------------
% PDE: div(rho*u*T) = div((k/cp)*grad T),  u=(u_x(y),0), T passif
% FVM structured grid. Convection: UD (base), QUICK (Q17) via deferred correction.
% Diffusion: CD. BC: Inlet Dirichlet, Outlet Neumann dT/dx=0,
%  Walls: Dirichlet (Q9..16), puis Neumann flux au DOF de mur (Q20..Q22).
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========= PARAMETRES =========a
rho = 1.0;                 % [kg/m^3]
cp  = 10.0;                % [J/(kg.K)]
k   = 0.12;                % [W/(m.K)]
nu  = 1e-2;                % [m^2/s]
Gamma = k / cp;            % [kg/(m.s)]   (consistant avec div(rho u T) = div((k/cp) grad T))
L = 10.0;                  % [m]
H = 1.0;                   % [m]
Pe_target = 16.5;

% Pe = (2*rho*H*u_mean)/Gamma  => u_mean = Pe*Gamma/(2*rho*H)
u_mean = Pe_target * Gamma / (2*rho*H);
Re     = 2*H*u_mean/nu;
fprintf('Q10: Gamma=%.6g [kg/(m.s)], u_mean=%.6g [m/s], Re=%.2f (laminar=%s)\n', ...
    Gamma,u_mean,Re, ternary(Re<2300,'Oui','Non'));

Tinlet = 50;             % [°C]
Twall  = 100;            % [°C]

nx = 50; ny = 5;
[xc,yc,dx,dy] = buildGrid(L,H,nx,ny);

u_prof = laminarUx_profile(yc,H,u_mean);

BC.walls       = 'Dirichlet';   % 'Dirichlet' or 'NeumannFlux'
BC.inlet       = 'Dirichlet';
BC.outlet      = 'Neumann0';
BC.q_wall_val  = 10.0;          % (Q20+)
BC.bottomFluxWindow = [];       % [x1 x2 q] for Q22

%% ========= Assemblage & solve (Q11) =========
[A,b,maskDir,Tpresc] = assembleSystem(nx,ny,dx,dy,u_prof,rho,Gamma,k, ...
                                      Tinlet,Twall,BC,'UD',[]);
T = A\b;
Tf = reshape(T,[ny,nx]);

%% Q12 quick checks
checkBCs_post(Tf,xc,yc,dx,dy,k,u_prof,Tinlet,Twall);

%% Q13
[To, Tc, Tmean, xe] = compute_profiles_and_xe(Tf,xc,yc,u_prof,Tinlet,Twall);
fprintf('Q13: x_e ≈ %.3f m (90%%)\n', xe);

%% Q14 (identique)
omegas=[1.0 1.5]; Ntot=numel(T); T0=100*ones(Ntot,1);
for w=omegas
    [~,hist]=SOR_solve(A,b,T0,w,1e-5,1e4);
    fprintf('Q14: SOR w=%.1f -> %d iters\n',w,hist.nit);
end

%% Q15 Nu_T
NuT = compute_Nu_T(Tf,xc,yc,dy,k,u_prof,Twall);
fprintf('Q15: Nu_T(outlet) ≈ %.3f (→7.54)\n', NuT(end));

%% Q16 (comme avant, non affiché ici)

%% ========= Q17: QUICK (deferred correction bornée) vs UD =========
meshes = [50 5; 100 11; 200 21; 400 41];
xe_UD = zeros(size(meshes,1),1);
xe_QK = zeros(size(meshes,1),1);

for m=1:size(meshes,1)
    nxm=meshes(m,1); nym=meshes(m,2);
    [xc,yc,dx,dy] = buildGrid(L,H,nxm,nym);
    u_prof = laminarUx_profile(yc,H,u_mean);

    % --- UD
    BC.walls='Dirichlet';
    [A,b,~,~] = assembleSystem(nxm,nym,dx,dy,u_prof,rho,Gamma,k,Tinlet,Twall,BC,'UD',[]);
    T_UD = A\b;
    Tf_UD = reshape(T_UD,[nym,nxm]);
    [~,~,~,xe_UD(m)] = compute_profiles_and_xe(Tf_UD,xc,yc,u_prof,Tinlet,Twall);

    % --- QUICK via deferred correction sur base UD
    Tprev = T_UD;                     % initialisation
    max_dc_iter = 5; tol_dc = 1e-6;
    for it=1:max_dc_iter
        [Aq,bq,~,~] = assembleSystem(nxm,nym,dx,dy,u_prof,rho,Gamma,k, ...
                                     Tinlet,Twall,BC,'QUICK',Tprev);
        Tnew = Aq\bq;
        if norm(Tnew-Tprev,2)/max(1e-16,norm(Tprev,2)) < tol_dc
            Tprev = Tnew; break;
        end
        Tprev = Tnew;
    end
    Tf_QK = reshape(Tprev,[nym,nxm]);
    [~,~,~,xe_QK(m)] = compute_profiles_and_xe(Tf_QK,xc,yc,u_prof,Tinlet,Twall);
end

xe_ref_UD = xe_UD(end); xe_ref_QK = xe_QK(end);
rel_UD = abs(xe_UD - xe_ref_UD)/xe_ref_UD*100;
rel_QK = abs(xe_QK - xe_ref_QK)/xe_ref_QK*100;

figure('Color','w');
subplot(1,2,1); plot(meshes(:,1),xe_UD,'-ob',meshes(:,1),xe_QK,'-sr'); grid on;
xlabel('n_x'); ylabel('x_e [m]'); title('x_e — UD vs QUICK (DC)');
legend('UD','QUICK-DC','Location','best');
subplot(1,2,2); plot(meshes(:,1),rel_UD,'-ob',meshes(:,1),rel_QK,'-sr'); grid on;
xlabel('n_x'); ylabel('erreur relative x_e [%]'); title('Convergence');
legend('UD','QUICK-DC','Location','best');

%exportgraphics(gcf, 'NFSQ17.png', 'Resolution', 300);

%% ========= Q18: maillages anisotropes (asse = numero di CV) =========
% Mesh families:
%  - varY: n_x fisso (400), cambiamo n_y
%  - varX: n_y fisso (41),  cambiamo n_x
meshes_varY = [400 5; 400 11; 400 21];
meshes_varX = [100 41; 200 41; 400 41];

ny_list = meshes_varY(:,2);   % asse per varY
nx_list = meshes_varX(:,1);   % asse per varX

xe_varY   = zeros(size(ny_list));
xe_varX   = zeros(size(nx_list));
NuMin_varY = zeros(size(ny_list));
NuMin_varX = zeros(size(nx_list));

% ---- famiglia con n_x=400, vario n_y ----
BC.walls = 'Dirichlet';
for m = 1:numel(ny_list)
    nxm = 400; nym = ny_list(m);
    [xc,yc,dx,dy] = buildGrid(L,H,nxm,nym);
    u_prof = laminarUx_profile(yc,H,u_mean);

    [A,b,~,~] = assembleSystem(nxm,nym,dx,dy,u_prof,rho,Gamma,k,Tinlet,Twall,BC,'UD',[]);
    Tf = reshape(A\b,[nym,nxm]);

    % x_e
    [~,~,~,xe_varY(m)] = compute_profiles_and_xe(Tf,xc,yc,u_prof,Tinlet,Twall);
    % min Nu_T(x)
    NuT = compute_Nu_T(Tf,xc,yc,dy,k,u_prof,Twall);
    NuMin_varY(m) = NuT(end-50);
end

% ---- famiglia con n_y=41, vario n_x ----
for m = 1:numel(nx_list)
    nxm = nx_list(m); nym = 41;
    [xc,yc,dx,dy] = buildGrid(L,H,nxm,nym);
    u_prof = laminarUx_profile(yc,H,u_mean);

    [A,b,~,~] = assembleSystem(nxm,nym,dx,dy,u_prof,rho,Gamma,k,Tinlet,Twall,BC,'UD',[]);
    Tf = reshape(A\b,[nym,nxm]);

    % x_e
    [~,~,~,xe_varX(m)] = compute_profiles_and_xe(Tf,xc,yc,u_prof,Tinlet,Twall);
    % min Nu_T(x)
    NuT = compute_Nu_T(Tf,xc,yc,dy,k,u_prof,Twall);
    NuMin_varX(m) = NuT(end-50);
end

% ---- Figure A: x_e vs n_y (n_x=400) ----
figure('Color','w');
subplot(2,2,1);
plot(ny_list, xe_varY,'-o'); grid on;
xlabel('n_y');
ylabel('x_e [m]');
title('x_e vs n_y (n_x = 400)');

% ---- Figure B: x_e vs n_x (n_y=41) ----
subplot(2,2,2);
plot(nx_list, xe_varX,'-s'); grid on;
xlabel('n_x');
ylabel('x_e [m]');
title('x_e vs n_x (n_y = 41)');

% ---- Figure C: min Nu_T(x) vs n_y (n_x=400) ----
subplot(2,2,3);
plot(ny_list, NuMin_varY,'-o'); grid on;
xlabel('n_y');
ylabel('min_x Nu_T(x) [-]');
title('min Nu_T(x) vs n_y (n_x = 400)');

% ---- Figure D: min Nu_T(x) vs n_x (n_y=41) ----
subplot(2,2,4);
plot(nx_list, NuMin_varX,'-s'); grid on;
xlabel('n_x');
ylabel('min_x Nu_T(x) [-]');
title('min Nu_T(x) vs n_x (n_y = 41)');

%exportgraphics(gcf, 'q18_nxny.png', 'Resolution', 300);
%% ========= Q19: Pe = 50 et 100 (UD) =========
for PeV = [50 100]
    u_mean_hi = PeV * Gamma / (2*rho*H);
    [xc,yc,dx,dy] = buildGrid(L,H,200,21);
    u_prof = laminarUx_profile(yc,H,u_mean_hi);
    BC.walls='Dirichlet';
    [A,b,~,~] = assembleSystem(200,21,dx,dy,u_prof,rho,Gamma,k,Tinlet,Twall,BC,'UD',[]);
    Tf = reshape(A\b,[21,200]);
    figure('Color','w'); imagesc(xc,yc,Tf); set(gca,'YDir','normal'); colorbar;
    title(sprintf('Q19: T(x,y) — Pe=%g (UD, 200x21)',PeV)); xlabel('x [m]'); ylabel('y [m]');
end

%% ========= Q20: murs flux imposé (Neumann au DOF de mur), Nu_q =========
BC_flux = BC; BC_flux.walls='NeumannFlux'; BC_flux.q_wall_val=10.0; BC_flux.bottomFluxWindow=[];
[xc,yc,dx,dy] = buildGrid(L,H,400,41);
u_prof = laminarUx_profile(yc,H,u_mean);
[A,b,~,~] = assembleSystem(400,41,dx,dy,u_prof,rho,Gamma,k,Tinlet,Twall,BC_flux,'UD',[]);
Tf20 = reshape(A\b,[41,400]);

% --- Plot 2D della temperatura T(x,y) ---
figure('Color','w');
imagesc(xc, yc, Tf20);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title('Temperature field T(x,y) with Neumann boundary conditions (q_{wall} = 10 W/m²)');
caxis([min(Tf20(:)) max(Tf20(:))]);  % range colori

%exportgraphics(gcf, 'q19TempField.png', 'Resolution', 300);

% --- Check flux ---
[qb,qt,eb,et] = checkBCs_flux_post(Tf20, xc, yc, dy, k, BC_flux.q_wall_val);
fprintf('Q20: check flux -> max|err| bottom=%.3g, top=%.3g\n', eb, et);

% --- Plot Nu_q(x) ---
Nuq = compute_Nu_q(Tf20, xc, yc, dy, k, u_mean, H, BC_flux.q_wall_val);
figure('Color','w'); 
plot(xc(3:end), Nuq(3:end), '-o', 'LineWidth', 1.5); 
grid on;
xlabel('x [m]'); 
ylabel('Nu_q(x) [-]'); 
title('Nu_q(x) along the lower wall (→ 8.24 downstream)');
ylim([0 12]);  % range ragionevole per Nu
hold on;
yline(8.24, '--r', 'LineWidth', 2, 'Label', 'Theoretical limit = 8.24');
hold off;

fprintf('\nOK: Q17–Q20 avec QUICK-DC borné & BC flux correctes.\n');


%% ========= Q21: Bilan d'énergie global =========
[Qin,Qout,Qwalls,imb] = global_energy_balance(Tf20, yc, u_mean, H, rho, cp, BC_flux.q_wall_val, L);
fprintf('Q_in = %.2f W, Q_out = %.3g W, Q_walls = %.3g W, Imbalance = %.3e W\n', Qin,Qout,Qwalls,imb);

%% ========= Q22: Riscaldamento localizzato (2 ≤ x ≤ 5) =========
% Obiettivo: trovare q_wall per ottenere T_outlet_mean = 60°C, 70°C, 80°C
% (cioè ΔT = +10°C, +20°C, +30°C rispetto a T_inlet = 50°C)

target_temps = [60, 70, 80];  % [°C]
q_values = zeros(size(target_temps));

fprintf('Localised heating (2 ≤ x ≤ 5 m) ===\n');

for idx = 1:length(target_temps)
    T_target = target_temps(idx);
    
    % Ricerca iterativa del valore di q_wall
    % Metodo: bisezione o Newton semplificato
    q_min = 0;
    q_max = 100;  % valore massimo ragionevole
    tol = 0.1;    % tolleranza sulla temperatura [°C]
    max_iter = 30;
    
    for iter = 1:max_iter
        q_test = (q_min + q_max) / 2;
        
        % Setup BC con riscaldamento localizzato
        BC_q22 = BC;
        BC_q22.walls = 'NeumannFlux';
        BC_q22.q_wall_val = 0.0;  % default (isolato)
        BC_q22.bottomFluxWindow = [2.0, 5.0, q_test];  % [x1, x2, q]
        
        % Risolvi
        [xc, yc, dx, dy] = buildGrid(L, H, 200, 21);
        u_prof = laminarUx_profile(yc, H, u_mean);
        [A, b, ~, ~] = assembleSystem(200, 21, dx, dy, u_prof, rho, Gamma, k, ...
                                      Tinlet, Twall, BC_q22, 'UD', []);
        Tf_q22 = reshape(A\b, [21, 200]);
        
        % Calcola temperatura media all'outlet
        T_outlet_mean = trapz(yc, u_prof .* Tf_q22(:, end)') / trapz(yc, u_prof);
        
        % Aggiorna limiti per bisezione
        if T_outlet_mean < T_target
            q_min = q_test;
        else
            q_max = q_test;
        end
        
        % Check convergenza
        if abs(T_outlet_mean - T_target) < tol
            q_values(idx) = q_test;
            fprintf('T_outlet = %.1f°C: q_wall = %.2f W/m² (iter=%d)\n', ...
                    T_target, q_test, iter);
            break;
        end
        
        if iter == max_iter
            q_values(idx) = q_test;
            fprintf('T_outlet = %.1f°C: q_wall ≈ %.2f W/m² (max iter, T=%.2f°C)\n', ...
                    T_target, q_test, T_outlet_mean);
        end
    end
end

%% Plot delle soluzioni per Q22
figure('Color', 'w');
for idx = 1:length(target_temps)
    T_target = target_temps(idx);
    q_opt = q_values(idx);
    
    % Risolvi con q ottimale
    BC_q22 = BC;
    BC_q22.walls = 'NeumannFlux';
    BC_q22.q_wall_val = 0.0;
    BC_q22.bottomFluxWindow = [2.0, 5.0, q_opt];
    
    [xc, yc, dx, dy] = buildGrid(L, H, 200, 21);
    u_prof = laminarUx_profile(yc, H, u_mean);
    [A, b, ~, ~] = assembleSystem(200, 21, dx, dy, u_prof, rho, Gamma, k, ...
                                  Tinlet, Twall, BC_q22, 'UD', []);
    Tf_q22 = reshape(A\b, [21, 200]);
    
    % Subplot: campo di temperatura
    subplot(2, 3, idx);
    imagesc(xc, yc, Tf_q22);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('x [m]');
    ylabel('y [m]');
    title(sprintf('T(x,y): T_{out}=%.0f°C, q=%.1f W/m²', T_target, q_opt));
    % Evidenzia zona riscaldata
    hold on;
    rectangle('Position', [2, 0, 3, H], 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
    hold off;
    
    % Subplot: temperatura media T_mean(x)
    subplot(2, 3, 3 + idx);
    den = trapz(yc, u_prof);
    Tmean_x = zeros(1, 200);
    for j = 1:200
        Tmean_x(j) = trapz(yc, u_prof .* Tf_q22(:, j)') / den;
    end
    plot(xc, Tmean_x, '-b', 'LineWidth', 1.5);
    grid on;
    xlabel('x [m]');
    ylabel('T_{mean}(x) [°C]');
    title(sprintf('T_{mean}(x): outlet = %.1f°C', Tmean_x(end)));
    % Evidenzia zona riscaldata
    xline(2, '--r', 'LineWidth', 1.5);
    xline(5, '--r', 'LineWidth', 1.5);
    yline(T_target, '--g', 'Target', 'LineWidth', 1.5);
end
exportgraphics(gcf, 'q21.png', 'Resolution', 300);

% Tabella riassuntiva
fprintf('\n=== Table Q22: Required q_wall values ===\n');
fprintf('ΔT [°C] | T_outlet [°C] | q_wall [W/m²]\n');
fprintf('--------|---------------|---------------\n');
for idx = 1:length(target_temps)
    fprintf('  +%.0f   |     %.0f       |    %.2f\n', ...
            target_temps(idx) - Tinlet, target_temps(idx), q_values(idx));
end
fprintf('==============================================\n\n');

%% =================== FONCTIONS ===================
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end

function [xc,yc,dx,dy] = buildGrid(L,H,nx,ny)
    xc = linspace(0,L,nx);
    yc = linspace(0,H,ny);
    dx = L/(nx-1);
    dy = H/(ny-1);
end

function u = laminarUx_profile(yc,H,u_mean)
    u = (3/2)*u_mean*(1 - ( (2*yc/H - 1).^2 ));
end

function [A,b,maskDir,Tpresc] = assembleSystem(nx,ny,dx,dy,u_prof,rho,Gamma,k, ...
                                               Tinlet,Twall,BC,scheme,Tprev)
    % Matrix = UD. Si scheme='QUICK', on ajoute une deferred correction bornée
    % au RHS, évaluée avec Tprev (vecteur N) en évitant tout index hors borne.
    N = nx*ny;
    A = spalloc(N,N,7*N);
    b = zeros(N,1);
    maskDir = false(N,1);
    Tpresc  = zeros(N,1);

    useDC = strcmpi(scheme,'QUICK') && ~isempty(Tprev);
    if useDC
        Tpf = reshape(Tprev,[ny,nx]);     % champ précédent
        tmin = min(Tinlet,Twall);         % bornes physiques sûres
        tmax = max(Tinlet,Twall);
    end

    paroisDir = strcmpi(BC.walls,'Dirichlet');

    for j=1:nx
        for i=1:ny
            n = (j-1)*ny + i;
            isInlet  = (j==1);
            isOutlet = (j==nx);
            isBottom = (i==1);
            isTop    = (i==ny);

            % --- Inlet (prioritaire)
            if isInlet
                maskDir(n)=true; Tpresc(n)=Tinlet;
                A(n,n)=1; b(n)=Tinlet; continue;
            end

            % --- Outlet Neumann dT/dx=0
            if isOutlet
                nL = (j-2)*ny + i;
                A(n,n)=1; A(n,nL)=-1; b(n)=0; continue;
            end

            % --- Parois bas/haut
            if (isBottom || isTop)
                if paroisDir
                    maskDir(n)=true; Tpresc(n)=Twall;
                    A(n,n)=1; b(n)=Twall; continue;
                else
                    % Neumann flux AU DOF DE MUR (Q20+)
                    xP = (j-1)*dx;
                    q_local = BC.q_wall_val; % par défaut
                    if isBottom && ~isempty(BC.bottomFluxWindow)
                        x1=BC.bottomFluxWindow(1); x2=BC.bottomFluxWindow(2); qwin=BC.bottomFluxWindow(3);
                        if xP>=x1 && xP<=x2, q_local = qwin; else, q_local = 0.0; end
                    end
                    if isBottom
                        nN = n + 1;  % -k/dy*T(1,j) + k/dy*T(2,j) = q
                        A(n,n)  = k/dy; 
                        A(n,nN) = -k/dy;  
                        b(n) = q_local;
                    else
                        nS = n - 1;  %  k/dy*T(ny,j) - k/dy*T(ny-1,j) = q
                        A(n,n)  =  k/dy;  A(n,nS) = -k/dy;  b(n) = q_local;
                    end
                    continue;
                end
            end

            % --- Cellule intérieure
            De = Gamma*dy/dx; Dw=De; Dn = Gamma*dx/dy; Ds=Dn;
            uP = u_prof(i);
            Fe = rho*uP*dy; Fw = Fe;

            nW = n - ny; nE = n + ny; nS = n - 1; nN = n + 1;
            aP = 0;

            % Diffusion
            aP = aP + Dw; A(n,nW) = A(n,nW) - Dw;
            aP = aP + De; A(n,nE) = A(n,nE) - De;
            aP = aP + Ds; A(n,nS) = A(n,nS) - Ds;
            aP = aP + Dn; A(n,nN) = A(n,nN) - Dn;

            % Convection UD (base M-matrix)
            if Fe > 0
                aP = aP + Fw; A(n,nW) = A(n,nW) - Fw;   % face Ouest ← T_W
            else
                aP = aP - Fe; A(n,nE) = A(n,nE) + Fe;   % face Est  ← T_E  (Fe<0)
            end

            % QUICK deferred correction (bornée) → RHS
            if useDC
                % Centres connus du champ précédent (sécurité bords assurée par tests)
                TP = Tpf(i, j);
                % j∈[2..nx-1] ici (ni inlet ni outlet) ⇒ TW et TE existent
                TW = Tpf(i, j-1);
                TE = Tpf(i, j+1);

                % Points WW/EE avec garde explicite (pas d’index 0 ou nx+1)
                if j >= 3
                    TWW = Tpf(i, j-2);
                else
                    TWW = TW;           % fallback amont
                end
                if j <= nx-2
                    TEE = Tpf(i, j+2);
                else
                    TEE = TE;           % fallback aval
                end

                % UD face values (à partir de Tprev)
                if Fe > 0
                    phiW_UD = TW;   phiE_UD = TP;
                    % QUICK non-borné (formule) puis bornage
                    phiW_Q = (3*TP + 6*TW - TWW)/8;
                    phiE_Q = (3*TE + 6*TP - TW)/8;
                else
                    phiW_UD = TP;   phiE_UD = TE;
                    phiW_Q = (3*TW + 6*TP - TE)/8;
                    phiE_Q = (3*TP + 6*TE - TEE)/8;
                end
                % Bornage physique
                phiW_Q = min(max(phiW_Q, tmin), tmax);
                phiE_Q = min(max(phiE_Q, tmin), tmax);

                % Correction “déférée” ajoutée à b (divergence des flux convectifs)
                delta_b = Fe*(phiE_Q - phiE_UD) - Fw*(phiW_Q - phiW_UD);
                b(n) = b(n) + delta_b;
            end

            A(n,n) = aP;
        end
    end
end

function checkBCs_post(Tf,xc,yc,dx,dy,k,u_prof,Tinlet,Twall)
    fprintf('Q12: checks BCs → inlet[%.2f,%.2f], bottom[%.2f,%.2f], top[%.2f,%.2f]\n', ...
        min(Tf(:,1)),max(Tf(:,1)), min(Tf(1,:)),max(Tf(1,:)), min(Tf(end,:)),max(Tf(end,:)));
    outDiff = Tf(:,end)-Tf(:,end-1);
    fprintf('     outlet max |ΔT| = %.3e\n', max(abs(outDiff)));
end

function [To,Tc,Tmean,xe] = compute_profiles_and_xe(Tf,xc,yc,u_prof,Tinlet,Twall)
    ny = numel(yc); nx = numel(xc);
    To = Tf(:,end);
    ic = round((ny+1)/2);
    Tc = Tf(ic,:);
    den = trapz(yc,u_prof);
    Tmean = zeros(1,nx);
    for j=1:nx
        Tmean(j) = trapz(yc, u_prof .* Tf(:,j)')/den;
    end
    ratio = (Tc - Tinlet)/(Twall - Tinlet);
    j90 = find(ratio>=0.9,1);
    if isempty(j90), xe = xc(end);
    else
        if j90==1, xe = xc(1);
        else
            x0=xc(j90-1); x1=xc(j90);
            r0=ratio(j90-1); r1=ratio(j90);
            xe = x0 + (0.9-r0)*(x1-x0)/(r1-r0);
        end
    end
end

function [T_sor, hist] = SOR_solve(A,b,T0,omega,tol,maxit)
    T = T0; N = numel(T);
    hist.err=[]; hist.res=[]; hist.nit=0;
    for k=1:maxit
        Told=T;
        for n=1:N
            sumL = A(n,1:n-1)*T(1:n-1);
            sumU = A(n,n+1:end)*Told(n+1:end);
            T_GS = ( b(n) - sumL - sumU ) / A(n,n);
            T(n) = (1-omega)*Told(n) + omega*T_GS;
        end
        rel = norm(T-Told,2)/max(1e-16,norm(Told,2));
        rr  = A*T - b; nres = norm(rr,2)/max(1e-16,norm(diag(A).*T,2));
        hist.err(end+1)=rel; hist.res(end+1)=nres;
        if rel<tol && nres<tol, hist.nit=k; break; end
    end
    if hist.nit==0, hist.nit=maxit; end
    T_sor=T;
end

function NuT = compute_Nu_T(Tf,xc,yc,dy,k,u_prof,Twall)
    nx=numel(xc); den = trapz(yc,u_prof);
    Tmean = zeros(1,nx);
    for j=1:nx
        Tmean(j) = trapz(yc, u_prof .* Tf(:,j)')/den;
    end
    NuT = zeros(1,nx);
    for j=1:nx
        dTdy = (Tf(2,j) - Tf(1,j))/dy;
        q = -k*dTdy;
        NuT(j) = (2*(yc(end)-yc(1)) * q)/(k*(Twall - Tmean(j)));
    end
end

function Nuq = compute_Nu_q(Tf, xc, yc, dy, k, u_mean, H, qwall)
    u_prof = laminarUx_profile(yc, H, u_mean);
    nx = numel(xc); 
    den = trapz(yc, u_prof);
    
    % Calcolo Tmean(x)
    Tmean = zeros(1, nx);
    for j = 1:nx
        Tmean(j) = trapz(yc, u_prof .* Tf(:,j)') / den;
    end
    
    Nuq = zeros(1, nx);
    for j = 1:nx
        % Il flusso qwall è definito come entrante nel fluido (+y direction)
        % Per la formula Nu serve il flusso dalla parete, quindi stesso segno
        q_bottom = qwall;
        
        % Temperatura alla parete (lower wall)
        Tw = Tf(1, j);
        
        % Formula Nu_q = (2H * q) / (k * (Tw - Tmean))
        % ATTENZIONE: se Tw < Tmean, abbiamo denominatore negativo!
        denom = k * (Tw - Tmean(j));
        
        if abs(denom) < 1e-12
            Nuq(j) = NaN;  % evita divisione per zero
        else
            Nuq(j) = (2*H * q_bottom) / denom;
        end
    end
end

function [q_bottom, q_top, err_bottom, err_top] = checkBCs_flux_post(Tf, xc, yc, dy, k, qwall)
    [ny, nx] = size(Tf);
    q_bottom = zeros(1, nx);
    for j = 1:nx
        dTdy = (Tf(2,j) - Tf(1,j)) / dy;
        q_bottom(j) = -k * dTdy;
    end
    q_top = zeros(1, nx);
    for j = 1:nx
        dTdy = (Tf(end,j) - Tf(end-1,j)) / dy;
        q_top(j) =  k * dTdy;
    end
    err_bottom = max(abs(q_bottom - qwall), [], 'omitnan');
    err_top    = max(abs(q_top    - qwall), [], 'omitnan');
end

function [Qin,Qout,Qwalls,imb] = global_energy_balance(Tf, yc, u_mean, H, rho, cp, qwall, L)
   u_prof = laminarUx_profile(yc, H, u_mean);
    
    % Portata massica per unità di profondità [kg/(s·m)]
    m_dot = rho * trapz(yc, u_prof);
    
    % Temperature medie pesate con velocità all'inlet e outlet
    Tin_bulk  = trapz(yc, u_prof .* Tf(:,1)')  / trapz(yc, u_prof);
    Tout_bulk = trapz(yc, u_prof .* Tf(:,end)')/ trapz(yc, u_prof);
    
    % Flussi energetici convettivi [W/m]
    Qin  = m_dot * cp * Tin_bulk;
    Qout = m_dot * cp * Tout_bulk;
    
    % Calore dalle pareti [W/m]
    % Q_walls = integrale di q lungo entrambe le pareti
    % Per parete uniforme: Q_bottom + Q_top = qwall * L + qwall * L
    Qwalls = 2 * qwall * L;   % [W/m] per unità di profondità
    
    % Bilancio: Qin + Qwalls - Qout dovrebbe essere ≈ 0
    imb = Qin + Qwalls - Qout;
end