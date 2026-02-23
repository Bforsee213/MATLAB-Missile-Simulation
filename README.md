% missile_intercept_3D Sim
% 3D missile vs evasive target intercept simulation
clearvars; close all; clc;

%% === Config (change variables as needed) ===
rng('shuffle');

% simulation
dt = 0.02;           % time step (s)
t_max = 200;         % simulation duration (s)

% spawn constraints
spawn_min_dist = 8000;    % minimum distance from missile (origin) (m)
spawn_max_dist = 10000;   % maximum distance from missile (origin) (m)
min_alt = 500;            % minimum altitude (m)
max_alt = 10000;          % maximum altitude (m)

% target (evader) kinematic params (3D-capable)
target_init_speed = 300;        % initial speed (m/s)
target_speed_min = 50;          % minimum possible speed (m/s)
max_target_speed = 500;         % **absolute max speed** (m/s)
target_long_accel_max = 20;     % max longitudinal accel (m/s^2) (how fast target can speed up / slow down)
target_max_turn_deg = 120;       % deg/s (heading turn-rate limit)
target_max_turn = deg2rad(target_max_turn_deg);

% evasion behavior
evasion_detect_radius = 3000;   % m - when missile inside this, target performs stronger 3D evasion
evasion_gain = 1.4;
vertical_evasion_bias = 0.8;

% missile parameters
missile_speed = 750;            % m/s (assumed constant for missile)
missile_max_turn_deg = 30;      % deg/s
missile_max_turn = deg2rad(missile_max_turn_deg);

% intercept
hit_radius = 5;                % meters considered a hit

% visualization
show_quivers = true;
pause_scale = 0.5;             % animation speed (smaller -> faster)
quiver_scale = 100;             % length multiplier for heading quivers

%% === Initial spawn (3D) with radial distance constraint ===
max_attempts = 20000;
attempt = 0;
target.pos = [];
while attempt < max_attempts
    attempt = attempt + 1;
    dir = randn(1,3); dir = dir / norm(dir);
    r = spawn_min_dist + (spawn_max_dist - spawn_min_dist)*rand;
    pos = r * dir;
    if pos(3) >= min_alt && pos(3) <= max_alt
        target.pos = pos;
        break;
    end
end
if isempty(target.pos)
    warning('Fallback spawn used.');
    target.pos = [spawn_min_dist, 0, (min_alt+max_alt)/2];
end

% initial heading roughly away from origin plus small noise
init_dir = unit(target.pos - [0 0 0]) + 0.2*randn(1,3);
target.heading = unit(init_dir);

% target speed state (variable)
target.speed = target_init_speed;  % current scalar speed (m/s)

% initial velocity vector
target.vel = target.speed * target.heading;

% missile initial state
missile.pos = [0 0 0];
missile.heading = unit(target.pos - missile.pos);
missile.vel = missile_speed * missile.heading;

% preallocate trajectories
T = 0:dt:t_max;
Nt = numel(T);
tar_traj = nan(Nt,3);
mis_traj = nan(Nt,3);

hit = false; hit_t = NaN; hit_idx = NaN;

%% === Simulation loop ===
for k = 1:Nt
    t = T(k);
    tar_traj(k,:) = target.pos;
    mis_traj(k,:) = missile.pos;
    
    R = target.pos - missile.pos;
    dist = norm(R);
    if dist <= hit_radius
        hit = true; hit_t = t; hit_idx = k;
        fprintf('Intercept at t = %.2f s, separation = %.2f m\n', t, dist);
        break;
    end
    
    % ---------- TARGET EVASION (3D) ----------
    % baseline small random lateral/3D acceleration (affects heading only)
    a_rand = 6 * unit(randn(1,3));
    if norm(target.vel) > 1e-6
        a_rand = a_rand - (dot(a_rand,target.vel)/norm(target.vel)^2)*target.vel;
    end
    
    if dist <= evasion_detect_radius
        away = unit(target.pos - missile.pos);
        perp = cross(away, unit(randn(1,3)));
        if norm(perp) < 1e-6
            perp = null(away)'; perp = perp(:,1)';
        end
        perp = unit(perp);
        vertical_component = [0 0 away(3)] * vertical_evasion_bias;
        evade_dir = unit( evasion_gain*away + 0.9*perp + 0.6*vertical_component + 0.4*randn(1,3) );
        desired_heading = evade_dir;
        % also bias longitudinal accel when evading: try to speed up to escape
        desired_speed_change = +30; % m/s (target attempts to increase speed when evading)
    else
        desired_heading = unit( target.heading + 0.02*randn(1,3) );
        % occasionally vary speed slightly when not evading
        desired_speed_change = (-5) + 10*rand; % in range [-5, +5] m/s when not evading
    end
    
    % --------- Enforce target turn-rate (3D) ----------
    cur_dir = target.heading;
    cosang_t = max(min(dot(cur_dir, desired_heading),1),-1);
    ang_t = acos(cosang_t);
    max_rot_t = target_max_turn * dt;
    if ang_t <= 1e-12
        new_heading = desired_heading;
    elseif ang_t <= max_rot_t
        new_heading = desired_heading;
    else
        axis_t = cross(cur_dir, desired_heading);
        if norm(axis_t) < 1e-12
            axis_t = null(cur_dir)'; axis_t = axis_t(:,1)';
        end
        axis_t = unit(axis_t);
        new_heading = rodrigue_rotation(cur_dir, axis_t, max_rot_t);
        new_heading = unit(new_heading);
    end
    target.heading = unit(new_heading);
    
    % --------- Longitudinal (speed) update ----------
    % simple longitudinal accel model: choose a_long within allowed bounds
    % make a_long try to achieve desired_speed_change over 1 second-ish (simple heuristic)
    target_desired_speed = target.speed + desired_speed_change;
    % clamp desired target desired speed to realistic range before computing accel
    target_desired_speed = min(max(target_desired_speed, target_speed_min), max_target_speed);
    % compute needed accel to reach desired in ~1.0s (time-constant)
    tau = 1.0; % seconds to reach desired approximately
    a_long = (target_desired_speed - target.speed) / tau;
    % limit acceleration magnitude
    a_long = max(min(a_long, target_long_accel_max), -target_long_accel_max);
    % apply acceleration
    target.speed = target.speed + a_long * dt;
    % finally clamp to absolute speed bounds
    if target.speed > max_target_speed
        target.speed = max_target_speed;
    elseif target.speed < target_speed_min
        target.speed = target_speed_min;
    end
    
    % --------- update velocity & position ----------
    target.vel = target.speed * target.heading;
    target.pos = target.pos + target.vel * dt;
    
    % softly keep altitude inside band
    if target.pos(3) < min_alt
        target.pos(3) = min_alt;
        target.heading(3) = abs(target.heading(3)) + 0.05;
        target.heading = unit(target.heading);
        target.vel = target.speed * target.heading;
    elseif target.pos(3) > max_alt
        target.pos(3) = max_alt;
        target.heading(3) = -abs(target.heading(3)) - 0.05;
        target.heading = unit(target.heading);
        target.vel = target.speed * target.heading;
    end
    
    % ---------- MISSILE (3D pure pursuit with turn limit) ----------
    desired_dir = unit( target.pos - missile.pos );
    current_dir = missile.heading;
    cosang = max(min(dot(current_dir, desired_dir),1),-1);
    ang = acos(cosang);
    max_rot = missile_max_turn * dt;
    if ang <= 1e-12
        new_dir = desired_dir;
    elseif ang <= max_rot
        new_dir = desired_dir;
    else
        axis = cross(current_dir, desired_dir);
        if norm(axis) < 1e-12
            axis = null(current_dir)'; axis = axis(:,1)';
        end
        axis = unit(axis);
        new_dir = rodrigue_rotation(current_dir, axis, max_rot);
        new_dir = unit(new_dir);
    end
    missile.heading = unit(new_dir);
    missile.vel = missile_speed * missile.heading;
    missile.pos = missile.pos + missile.vel * dt;
    
    % safety bailout
    if any(abs(missile.pos) > 5e6) || any(abs(target.pos) > 5e6)
        warning('Objects left simulation bounds - aborting.');
        break;
    end
end

% trim arrays
if hit
    tar_traj = tar_traj(1:hit_idx,:);
    mis_traj = mis_traj(1:hit_idx,:);
    T = T(1:hit_idx);
else
    last = find(~isnan(tar_traj(:,1)),1,'last');
    if isempty(last), last = 1; end
    tar_traj = tar_traj(1:last,:);
    mis_traj = mis_traj(1:last,:);
    T = T(1:last);
    fprintf('No intercept within %.1f s\n', t_max);
end

%% === 3D Plot & animation ===
figure('Color','w','Position',[200 150 1100 750]);
ax = axes('NextPlot','add'); hold(ax,'on'); grid(ax,'on'); view(3); axis equal;

allpos = [tar_traj; mis_traj; 0 0 0];
maxabs = max(abs(allpos(:)));
lim = max(500, maxabs*1.2);
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); zlim(ax, [0 max(max_alt, lim)]);

hTarPath = plot3(ax, nan, nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.85 0.1 0.1]);
hMisPath = plot3(ax, nan, nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.06 0.2 0.8]);
hTar = plot3(ax, tar_traj(1,1), tar_traj(1,2), tar_traj(1,3), 'o', 'MarkerSize',11,'MarkerFaceColor',[0.85 0.1 0.1]);
hMis = plot3(ax, mis_traj(1,1), mis_traj(1,2), mis_traj(1,3), 's', 'MarkerSize',11,'MarkerFaceColor',[0.06 0.2 0.8]);
hHit = plot3(ax, NaN, NaN, NaN, 'p', 'MarkerSize', 18, 'MarkerFaceColor',[0.9 0.6 0]);

if show_quivers
    hTarQ = quiver3(ax, 0,0,0, 0,0,0, 'MaxHeadSize',1.2, 'LineWidth',1.2, 'Color',[0.6 0.1 0.1]);
    hMisQ = quiver3(ax, 0,0,0, 0,0,0, 'MaxHeadSize',1.2, 'LineWidth',1.2, 'Color',[0.06 0.2 0.8]);
end

xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Evasive target (3D) — spawn dist in [%d, %d] m, z in [%d, %d] m — max speed = %d m/s', ...
    spawn_min_dist, spawn_max_dist, min_alt, max_alt, max_target_speed));

% small axis markers for reference
plot3(ax, [0 500], [0 0], [0 0], 'k-', 'LineWidth', 1); text(500,0,0,'X','FontWeight','bold');
plot3(ax, [0 0], [0 500], [0 0], 'k-', 'LineWidth', 1); text(0,500,0,'Y','FontWeight','bold');
plot3(ax, [0 0], [0 0], [0 500], 'k-', 'LineWidth', 1); text(0,0,500,'Z','FontWeight','bold');

rotate3d(ax, 'on');

nSteps = size(tar_traj,1);
for k = 1:nSteps
    set(hTarPath, 'XData', tar_traj(1:k,1), 'YData', tar_traj(1:k,2), 'ZData', tar_traj(1:k,3));
    set(hMisPath, 'XData', mis_traj(1:k,1), 'YData', mis_traj(1:k,2), 'ZData', mis_traj(1:k,3));
    set(hTar, 'XData', tar_traj(k,1), 'YData', tar_traj(k,2), 'ZData', tar_traj(k,3));
    set(hMis, 'XData', mis_traj(k,1), 'YData', mis_traj(k,2), 'ZData', mis_traj(k,3));
    
    if show_quivers
        if k>1
            tdir = unit(tar_traj(min(k,nSteps),:) - tar_traj(max(k-1,1),:));
        else
            tdir = unit(target.heading);
        end
        set(hTarQ, 'XData', tar_traj(k,1), 'YData', tar_traj(k,2), 'ZData', tar_traj(k,3), ...
            'UData', tdir(1)*quiver_scale, 'VData', tdir(2)*quiver_scale, 'WData', tdir(3)*quiver_scale);
        
        if k>1
            mdir = unit(mis_traj(min(k,nSteps),:) - mis_traj(max(k-1,1),:));
        else
            mdir = unit(missile.heading);
        end
        set(hMisQ, 'XData', mis_traj(k,1), 'YData', mis_traj(k,2), 'ZData', mis_traj(k,3), ...
            'UData', mdir(1)*quiver_scale, 'VData', mdir(2)*quiver_scale, 'WData', mdir(3)*quiver_scale);
    end
    
    if hit && k==hit_idx
        set(hHit, 'XData', tar_traj(k,1), 'YData', tar_traj(k,2), 'ZData', tar_traj(k,3));
        title(ax, sprintf('Intercept at t = %.2f s', hit_t), 'Color', [0.8 0.2 0]);
        drawnow; break;
    end
    
    drawnow limitrate;
    pause(dt * pause_scale);
end

if ~hit
    title(ax, 'No intercept within simulation time', 'Color', [0.2 0.2 0.2]);
end

%% === Helper functions ===
function v = unit(v)
    if isempty(v) || all(abs(v) < 1e-12)
        return
    end
    v = v / norm(v);
end

function vr = rodrigue_rotation(v, k, theta)
    % Rotate vector v about axis k (unit) by angle theta (rad) using Rodrigues' formula
    vr = v*cos(theta) + cross(k,v)*sin(theta) + k*(dot(k,v))*(1-cos(theta));
    vr = vr(:)';
end
