%% MDM2 GP1 - Phugoid Model

function [t,thta_SOLNs, thtadot_SOLNs, thtaddot_SOLNs, X_SOLNs, Xdot_SOLNs, Xddot_SOLNs] = phugoidModel(thta0, N)

%% SOLUTION PARAMATERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial Values
thta0vals = thta0;   % pitch, rad
thtadot0 = 0;   % angular velocity, rad/s

X0 = [0;0];     % displacement, m
Xdot0 = [50;0]; % velocity, m/s

% Interval
t0 = 0;         % interval start, s
tN = 30;        % interval end, s
N = N;      % number of steps
h = (tN-t0)/N;  % stepsize, s
t = linspace(t0,tN,N+1);

%% PLOT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 1; % linewidth
TitleSize = 18;
LabelSize = 16;
FigSize = 800;
LegLabels = append('$\theta_0=',string(thta0),'$ $rad$');
LegSize = 14;

%% CONSTANTS (all in SI units) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physical Conditions
g = 9.81;       
                % acceleration due to gravitym m^2/s
rho = 1.112;    
                % air density at 1000m, kg/m^3

% Aircraft 
l = 8.2;        
                % length of fuselage, m
d = 2.7432;     
                % diameter of fuselage, m
r = l*(3/4);    
                % distance from CoM to tail stabilisers, m
m = 900;        
                % mass, kg
Ic = 2130;      % based on C172 Skyhawk, and scaled to Cessna by mass (similar planes)
                % MoI for cylinder (with a hole) about CoM, kgm^2
A1 = 16;        
                % wing area, m^2
A2 = 2;         
                % tail stabiliser area, m^2
lam1 = 2*pi;
                % Cf:AoA ratio for wings
lam2 = 2*pi;    
                % Cf:AoA ratio for stabilisers
Cd = 1.28;      
                % coefficient of drag for elevators
alph = (2*g*m)/(rho*dot(Xdot0,Xdot0)*lam1*A1); 
                % wing AoI (calculated for equilibrium stability), rad

% Combined Constants

LAM1 = lam1*rho*A1/2;
LAM2 = lam2*rho*A2/2;
GAM = Cd*rho*r*A2/2; %NOT R^2 because of linear case when |theta|<1, only R^1



%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Air Direction
phi = @(Xdot) atan(Xdot(2)/Xdot(1));

% Forces
L1 = @(thta,Xdot) LAM1.*(alph + thta - phi(Xdot)) .* dot(Xdot,Xdot);
L2 = @(thta,Xdot) LAM2.*(thta - phi(Xdot)) .* dot(Xdot,Xdot);
F = @(thtadot) GAM .* thtadot.*(abs(r*thtadot).^(abs(thtadot)>1));

% ODEs
Xddot = { @(thta,thtadot,Xdot) -(1/m).*((L1(thta,Xdot) + L2(thta,Xdot)) .* sin(phi(Xdot)) + F(thtadot).*sin(thta)),...
          @(thta,thtadot,Xdot) (1/m).*((L1(thta,Xdot) + L2(thta,Xdot)) .* cos(phi(Xdot)) + F(thtadot).*cos(thta) - m*g) };

thtaddot = @(thta,thtadot,Xdot) (-r/Ic).*(L2(thta,Xdot) .* cos(thta - phi(Xdot)) + F(thtadot));



%% SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thta_SOLNs = [];
thtadot_SOLNs = [];
thtaddot_SOLNs = [];

X_SOLNs = [];
Xdot_SOLNs = [];
Xddot_SOLNs = [];

for thta0=thta0vals

    % Initialise

    thta_SOLN = thta0;
    thtadot_SOLN = thtadot0;
    thtaddot_SOLN = thtaddot(thta0,thtadot0,Xdot0);

    X_SOLN = X0;
    Xdot_SOLN = Xdot0;
    Xddot_SOLN = [Xddot{1}(thta0,thtadot0,Xdot0);Xddot{2}(thta0,thtadot0,Xdot0)];

    % Forward Euler Iteration 

    fEuler3 = @(f,prv,xprv,yprv,zprv) prv + h*(f(xprv,yprv,zprv));

    for n=2:N+1
        Xdot_SOLN(1,n) = fEuler3(Xddot{1},Xdot_SOLN(1,n-1),thta_SOLN(n-1),thtadot_SOLN(n-1),Xdot_SOLN(:,n-1));
        Xdot_SOLN(2,n) = fEuler3(Xddot{2},Xdot_SOLN(2,n-1),thta_SOLN(n-1),thtadot_SOLN(n-1),Xdot_SOLN(:,n-1));
        thtadot_SOLN(n) = fEuler3(thtaddot,thtadot_SOLN(n-1),thta_SOLN(n-1),thtadot_SOLN(n-1),Xdot_SOLN(:,n-1));
        thta_SOLN = cumtrapz(t(1:n),thtadot_SOLN(1:n),2)+thta0;
        X_SOLN = [ cumtrapz(t(1:n),Xdot_SOLN(1,1:n),2); cumtrapz(t(1:n),Xdot_SOLN(2,1:n),2) ];
        thtaddot_SOLN(n) = thtaddot(thta_SOLN(n),thtadot_SOLN(n),Xdot_SOLN(:,n));
        Xddot_SOLN(1,n) = Xddot{1}(thta_SOLN(n),thtadot_SOLN(n),Xdot_SOLN(:,n));
        Xddot_SOLN(2,n) = Xddot{2}(thta_SOLN(n),thtadot_SOLN(n),Xdot_SOLN(:,n));
    end
    
    thta_SOLNs = [thta_SOLNs; thta_SOLN];
    thtadot_SOLNs = [thtadot_SOLNs; thtadot_SOLN];
    thtaddot_SOLNs = [thtaddot_SOLNs; thtaddot_SOLN];

    X_SOLNs = [X_SOLNs; X_SOLN];
    Xdot_SOLNs = [Xdot_SOLNs; Xdot_SOLN];
    Xddot_SOLNs = [Xddot_SOLNs; Xddot_SOLN];
    
end

    %% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Picth Dynamics Plots

    ftheta = figure;
    ftheta.Position = [ 0 0 FigSize FigSize];
    tlTheta = tiledlayout(2,1);
    
    plotEnd = 3; %time to end pitch plots
    axTh1 = nexttile;
    plot(t,thta_SOLNs,'linewidth',w)     % pitch over time
    title('\textbf{Pitch - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\theta$, radians','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','location','southeast','fontsize',LegSize)
    
    axTh2 = nexttile;
    plot(t,thtadot_SOLNs,'linewidth',w)  % angular velocity over time
    title('\textbf{Angular Velocity - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\dot\theta$, radians s$^{-1}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','location','southeast','FontSize',LegSize)
    
    
    
    % Angular Accel
    fthetaddot = figure;
    fthetaddot.Position = [ 0 0 FigSize FigSize];
    
    axTh3S = nexttile;
    plot(t,thtaddot_SOLNs,'linewidth',w)  % angular accel over time
    title('\textbf{Angular Acceleration - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\ddot\theta$, radians s$^{-2}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)
    
    axTh3L = nexttile;
    plot(t,thtaddot_SOLNs,'linewidth',w)  % angular accel over time
    title('\textbf{Angular Acceleration - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\ddot\theta$, radians s$^{-2}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)
    
    % Velocity Plots

    fXdot = figure;
    fXdot.Position = [ 0 0 FigSize FigSize];
    tlXdot = tiledlayout(2,1);
    
    axVel1 = nexttile;
    plot(t,transpose(Xdot_SOLNs(1:2:end,:)),'linewidth',w) % x velocity over time
    title('\textbf{Horizontal Velocity - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\dot x$, metres s$^{-1}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','location','southeast','fontsize',LegSize)

    axVel2 = nexttile;
    plot(t,transpose(Xdot_SOLNs(2:2:end,:)),'linewidth',w) % y velocity over time
    title('\textbf{Vertical Velocity - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\dot y$, metres s$^{-1}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)
    
    % Accel Plots

    fXddot = figure;
    fXddot.Position = [ 0 0 FigSize FigSize];
    tlXddot = tiledlayout(2,1);
    
    axAcc1 = nexttile;
    plot(t,transpose(Xddot_SOLNs(1:2:end,:)),'linewidth',w) % x accel over time
    title('\textbf{Horizontal Acceleration - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\ddot x$, metres s$^{-2}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','location','southeast','fontsize',LegSize)

    axAcc2 = nexttile;
    plot(t,transpose(Xddot_SOLNs(2:2:end,:)),'linewidth',w) % y accel over time
    title('\textbf{Vertical Acceleration - Time}','Interpreter','latex','fontsize',TitleSize)
    xlabel('t, seconds','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\ddot y$, metres s$^{-2}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)

    
    % Trajectory Plot

    ftrajectory = figure;
    ftrajectory.Position = [ 0 0 FigSize FigSize];
    
    axTrajS = nexttile;
    plot(transpose(X_SOLNs(1:2:end,:)),transpose(X_SOLNs(2:2:end,:)),'linewidth',w) % trajectory of aircraft
    title('\textbf{Trajectory} $y$ - $x$','Interpreter','latex','fontsize',TitleSize)
    xlabel('$x$, metres','Interpreter','latex','fontsize',LabelSize)
    ylabel('$y$, metres','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)
    
    axTrajL = nexttile;
    plot(transpose(X_SOLNs(1:2:end,:)),transpose(X_SOLNs(2:2:end,:)),'linewidth',w) % trajectory of aircraft
    title('\textbf{Trajectory} $y$ - $x$','Interpreter','latex','fontsize',TitleSize)
    xlabel('$x$, metres','Interpreter','latex','fontsize',LabelSize)
    ylabel('$y$, metres','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)
    
    % Vector Space Plots

    vecX = figure;
    vecX.Position = [ 0 0 FigSize FigSize];
    axVecX = axes;
    plot(transpose(X_SOLNs(2:2:end,:)),transpose(Xdot_SOLNs(2:2:end,:)),'linewidth',w)
    title('\textbf{Vertical Velocity - Displacement }','Interpreter','latex','fontsize',TitleSize)
    xlabel('$y$, metres','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\dot y$, metres s$^{-1}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)
    
    vecThta = figure;
    vecThta.Position = [ 0 0 FigSize FigSize];
    axVecTh = axes;
    plot(transpose(thta_SOLNs),transpose(thtadot_SOLNs),'linewidth',w)
    title('\textbf{Angular Acceleration - Displacement}','Interpreter','latex','fontsize',TitleSize)
    xlabel('$\theta$, radians','Interpreter','latex','fontsize',LabelSize)
    ylabel('$\dot\theta$, radians s$^{-1}$','Interpreter','latex','fontsize',LabelSize)
    a = gca;
    a.TickLabelInterpreter = 'latex';
    grid on
    legend(LegLabels,'Interpreter','latex','fontsize',LegSize)
    
    
    %SHORT PERIOD scale - ALL SET
    axis(axTh1,[0 plotEnd -0.09 0.18])
    legend(axTh1,'location','northeast')
    
    axis(axTh2, [0 plotEnd -0.9 0.5])
    
    axis(axVel1, [0 4 49.82 50.02])
    legend(axVel1,'location','southeast')
    
    axis(axVel2,[0 4 -1.5 3])
    
    axis(axAcc1,[0 plotEnd -1.15 0.35])
    
    axis(axAcc2,[0 plotEnd -15 30])
    
    exportgraphics(tlTheta,'pitchplotshort.pdf','ContentType','vector')
    exportgraphics(tlXdot,'velplotshort.pdf','ContentType','vector')
    exportgraphics(tlXddot,'accplotshort.pdf','ContentType','vector')


    % PHUGOID scale - ALL SET
    axis(axTh1,[0 tN -4*10^-3 4*10^-3])
    legend(axTh1,'location','northeast')
    axis(axTh2, [0 tN -1*10^-3 1*10^-3])
    legend(axTh2,'location','northeast')
    
    axis(axVel1, [0 tN 49.85 50.15 ])
    legend(axVel1,'location','northeast')
    
    axis(axVel2,[0 tN -0.2 0.2])
    legend(axVel2,'location','northeast')
    
    axis(axAcc1,[0 tN -0.035 0.035])
    legend(axAcc1,'location','southeast')
    
    axis(axAcc2,[0 tN -0.05 0.05])
    legend(axAcc2,'location','northeast')
    
    exportgraphics(tlTheta,'pitchplotlong.pdf','ContentType','vector')
    exportgraphics(tlXdot,'velplotlong.pdf','ContentType','vector')
    exportgraphics(tlXddot,'accplotlong.pdf','ContentType','vector')
    
    % Short and Long in One
    
    axis(axTh3S,[0 plotEnd -8 5])
    legend(axTh3S,'location','southeast')
    axis(axTh3L, [0 tN -4*10^-4 4*10^-4])
    legend(axTh3L,'location','southeast')
    
    axis(axTrajS,[0 300 -0.05 0.85])
    legend(axTrajS,'location','northeast')
    axis(axTrajL,[0 1500 -0.6 0.85])
    legend(axTrajL,'location','north')
    
    exportgraphics(ftrajectory,'trajectoryplot.pdf','ContentType','vector')
    exportgraphics(fthetaddot,'angaccelplot.pdf','ContentType','vector')

    % VECTOR SPACE
    axis(axVecX,[-0.55 0.85 -1.35 2.9])
    legend(axVecX,'location','northwest')
    axis(axVecTh,[-0.08 0.17 -0.9 0.5])
    
    exportgraphics(vecX,'veldispplot.pdf','ContentType','vector')
    exportgraphics(vecThta,'pitchangvelplot.pdf','ContentType','vector')
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%END%