# 3D Missile Intercept Simulation (MATLAB)

This project simulates a 3D missile interception scenario where a missile pursues an evasive airborne target using turn rate limited pursuit. The target performs adaptive evasive maneuvers when the missile enters a detection radius and operates under realistic speed and acceleration constraints. The simulation outputs a real time 3D animation of both trajectories and reports whether an intercept occurs within the simulation time.

---

```matlab
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
target_long_accel_max = 20;     % max longitudinal accel (m/s^2)
target_max_turn_deg = 120;       % deg/s (heading turn-rate limit)
target_max_turn = deg2rad(target_max_turn_deg);

% evasion behavior
evasion_detect_radius = 3000;   % m
evasion_gain = 1.4;
vertical_evasion_bias = 0.8;

% missile parameters
missile_speed = 750;            % m/s (assumed constant)
missile_max_turn_deg = 30;      
missile_max_turn = deg2rad(missile_max_turn_deg);

% intercept
hit_radius = 5;                % meters considered a hit

% visualization
show_quivers = true;
pause_scale = 0.5;
quiver_scale = 100;

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

init_dir = unit(target.pos - [0 0 0]) + 0.2*randn(1,3);
target.heading = unit(init_dir);
target.speed = target_init_speed;
target.vel = target.speed * target.heading;

missile.pos = [0 0 0];
missile.heading = unit(target.pos - missile.pos);
missile.vel = missile_speed * missile.heading;

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
    
    % ---------- TARGET EVASION ----------
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
        desired_speed_change = +30;
    else
        desired_heading = unit( target.heading + 0.02*randn(1,3) );
        desired_speed_change = (-5) + 10*rand;
    end
    
    % turn-rate enforcement
    cur_dir = target.heading;
    cosang_t = max(min(dot(cur_dir, desired_heading),1),-1);
    ang_t = acos(cosang_t);
    max_rot_t = target_max_turn * dt;
    
    if ang_t <= max_rot_t
        new_heading = desired_heading;
    else
        axis_t = cross(cur_dir, desired_heading);
        if norm(axis_t) < 1e-12
            axis_t = null(cur_dir)'; axis_t = axis_t(:,1)';
        end
        axis_t = unit(axis_t);
        new_heading = rodrigue_rotation(cur_dir, axis_t, max_rot_t);
    end
    
    target.heading = unit(new_heading);
    
    % speed update
    target_desired_speed = target.speed + desired_speed_change;
    target_desired_speed = min(max(target_desired_speed, target_speed_min), max_target_speed);
    
    tau = 1.0;
    a_long = (target_desired_speed - target.speed) / tau;
    a_long = max(min(a_long, target_long_accel_max), -target_long_accel_max);
    
    target.speed = target.speed + a_long * dt;
    target.speed = min(max(target.speed, target_speed_min), max_target_speed);
    
    target.vel = target.speed * target.heading;
    target.pos = target.pos + target.vel * dt;
    
    % altitude bounds
    if target.pos(3) < min_alt
        target.pos(3) = min_alt;
        target.heading(3) = abs(target.heading(3)) + 0.05;
    elseif target.pos(3) > max_alt
        target.pos(3) = max_alt;
        target.heading(3) = -abs(target.heading(3)) - 0.05;
    end
    
    target.heading = unit(target.heading);
    target.vel = target.speed * target.heading;
    
    % ---------- MISSILE ----------
    desired_dir = unit( target.pos - missile.pos );
    current_dir = missile.heading;
    
    cosang = max(min(dot(current_dir, desired_dir),1),-1);
    ang = acos(cosang);
    max_rot = missile_max_turn * dt;
    
    if ang <= max_rot
        new_dir = desired_dir;
    else
        axis = cross(current_dir, desired_dir);
        if norm(axis) < 1e-12
            axis = null(current_dir)'; axis = axis(:,1)';
        end
        axis = unit(axis);
        new_dir = rodrigue_rotation(current_dir, axis, max_rot);
    end
    
    missile.heading = unit(new_dir);
    missile.vel = missile_speed * missile.heading;
    missile.pos = missile.pos + missile.vel * dt;
end

%% === Helper functions ===
function v = unit(v)
    if isempty(v) || all(abs(v) < 1e-12)
        return
    end
    v = v / norm(v);
end

function vr = rodrigue_rotation(v, k, theta)
    vr = v*cos(theta) + cross(k,v)*sin(theta) + k*(dot(k,v))*(1-cos(theta));
    vr = vr(:)';
end
```

---
