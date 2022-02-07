function ball_collision(balls,ball_pos,ball_v,varargin)
%% PARAMETROS
% Constantes
G_newton  = 0;          %Constante de Newton
g_grav    = 0;          %Aceleracion de la gravedad
k_coulomb = 0;          %Constante de Coulomb
Cr        = 1;          %Coeficiente de restitucion
n         =15;          %Numero de bolas (particulas)
dflt_axis_lim = [-10,10;-10,10;-10,10];
%% Default
GravityVec  = [0 0 -1];
GravitySrc  = [];
CoulombSrc  = [];
GroundPoint = [0 0 0];
grav_ball = 0;
friction_factor = 0;
axis_lim = dflt_axis_lim;
dt = 0.1;
tmax = 100;
t_fig = 50;
debug_mode = true;
dflt_layout = '3d';
show_legend = false;
show_grid = 'on';
w_pln_v = {};
w_sph_v = {};
pln_pnt = cell(1,6); %{[pnt1;pnt2;pnt3;pnt4],...}
%% caras de la caja:
pln_pnt{1} = [-10 -10  10;-10  10  10; 10  10  10; 10 -10  10];
pln_pnt{2} = [-10 -10 -10;-10  10 -10; 10  10 -10; 10 -10 -10];
pln_pnt{3} = [-10 -10 -10;-10 -10  10;-10  10  10;-10  10 -10];
pln_pnt{4} = [ 10 -10 -10; 10 -10  10; 10  10  10; 10  10 -10];
pln_pnt{5} = [-10  10 -10;-10  10  10; 10  10  10; 10  10 -10];
pln_pnt{6} = [-10 -10 -10;-10 -10  10; 10 -10  10; 10 -10 -10];
wall_pln = []; %[0 0 -1 10;0 0 -1 -10;-1 0 0 -10;-1 0 0 10;0 1 0 -10;0 1 0 10];
wall_sph = [];
ball_color = [];
pln_color = [];
sph_color = [];
% allocate input arguments
n_arg_in = nargin;
switch n_arg_in
    case 0
        var_add = {};
    case 1
        var_add = {balls};
    case 2
        var_add = {balls,ball_pos};
    otherwise
        var_add = {balls,ball_pos,ball_v};
end
k_var = 1; %##
while k_var <= length(var_add)
    if ischar(var_add{k_var})
        break
    end
    k_var = k_var+1;
end
switch k_var
    case 1
        balls = [];
        ball_pos = [];
        ball_v = [];
    case 2
        ball_pos = [];
        ball_v = [];
    case 3
        ball_v = [];
end
varargin = [var_add(k_var:end) varargin];
if rem(n_arg_in-k_var+1,2) == 1
    fprintf('Not enough input arguments.\nValue of "%s" is not specified.\n',...
        varargin{end})
    return
end
for kk = 1:2:length(varargin)
    switch upper(varargin{kk})
        case {'BALLS','BALL','BALL_PROPERTY','BALL_PROPERTIES',...
                'PROPERTY','PROPERTIES','PROP'}
            balls = varargin{kk+1};
        case {'BALL_POS','POS','POSITION'}
            ball_pos = varargin{kk+1};
        case {'BALL_V','BALL_VELOCITY','BALL_VELOCITIES',...
                'VELOCITY','VELOCITIES','V'}
            ball_v = varargin{kk+1};
        case {'BALL_COLOR','COLOR'}
            ball_color = varargin{kk+1};
        case {'PLN_COLOR','PLANE_COLOR','WALL_COLOR'}
            pln_color = varargin{kk+1};
        case {'SPH_COLOR','SPHERE_COLOR'}
            sph_color = varargin{kk+1};
        case {'WALL_PLANE','PLANE'}
            wall_pln = varargin{kk+1};
        case {'PLANE_POINTS','PLANE_POINT','PLN_PNT','PLANE_PNT','PLN_POINT'}
            pln_pnt = varargin{kk+1};
        case {'WALL_SPHERE','SPHERE'}
            wall_sph = varargin{kk+1};
        case {'PLANE_VELOCITY','PLANE_V'}
            w_pln_v = varargin{kk+1};
        case {'SPHERE_VELOCITY','SPHERE_V'}
            w_sph_v = varargin{kk+1};
        case 'GRAVITYVEC'
            GravityVec = varargin{kk+1};
        case {'GRAVITYSRC','GRAVITY SRC','GRAVITY_SRC',...
                'GRAVITYSOURCE','GRAVITY SOURCE','GRAVITY_SOURCE',...
                'GRAVSRC','GRAV SRC','GRAV_SRC'}
            GravitySrc = varargin{kk+1};
        case {'COULOMBSRC','COULOMB SRC','COULOMB_SRC',...
                'COULOMBSOURCE','COULOMB SOURCE','COULOMB_SOURCE'}
            CoulombSrc = varargin{kk+1};
        case {'GROUNDPOINT','GROUND POINT','GROUND_POINT'}
            GroundPoint = varargin{kk+1};
        case {'GRAV_BALL','BALL_GRAV','GRAV BALL','BALL GRAV',...
                'GRAVBALL','BALLGRAV'}
            if ischar(varargin{kk+1})
                switch upper(varargin{kk+1})
                    case 'ON'
                        grav_ball = 1;
                    case 'OFF'
                        grav_ball = 0;
                    otherwise
                        fprintf('''%s'' is not a property value of "debug_mode".\n',...
                            varargin{kk+1})
                        return
                end
            else
                grav_ball = 1;
            end
        case 'DT'
            dt = varargin{kk+1};
        case {'TMAX','T'}
            tmax = varargin{kk+1};
        case {'DEBUG_MODE','DEBUG MODE','DEBUGMODE','DEBUG'}
            if ischar(varargin{kk+1})
                switch upper(varargin{kk+1})
                    case 'ON'
                        debug_mode = true;
                    case 'OFF'
                        debug_mode = false;
                    otherwise
                        fprintf('''%s'' is not a property value of "debug_mode".\n',...
                            varargin{kk+1})
                        return
                end
            else
                debug_mode = varargin{kk+1};
            end
        case 'G_NEWTON'
            G_newton = varargin{kk+1};
        case 'G_GRAV'
            g_grav = varargin{kk+1};
        case 'G'
            switch varargin{kk}
                case 'G'
                    G_newton = varargin{kk+1};
                case 'g'
                    g_grav = varargin{kk+1};
            end
        case 'K_COULOMB'
            k_coulomb = varargin{kk+1};
        case 'K'
            switch varargin{kk}
                case 'k'
                    k_coulomb = varargin{kk+1};
                case 'K'
                    fprintf('"K" is not a property name.\nDid you mean "k"?\n')
                    return
            end
        case {'CR','E'}
            Cr = varargin{kk+1};
        case {'FRICTION','FRICTION_FACTOR','FRICTION FACTOR','FRICT'}
            friction_factor = varargin{kk+1};
        case {'DFLT_LAYOUT','DFLT LAYOUT','DEFAULT_LAYOUT',...
                'DEFAULT LAYOUT','LAYOUT'}
            if ischar(varargin{kk+1})
                switch upper(varargin{kk+1})
                    case {'2D','3D'}
                        dflt_layout = varargin{kk+1};
                    otherwise
                        fprintf('''%s'' is not a property value of "dflt_layout".\n',...
                            varargin{kk+1})
                        return
                end
            else
                fprintf('''%f'' is not a property value of "show_legend".\n',...
                    varargin{kk+1})
                return
            end
        case {'SHOW_LEGEND','LEGEND'}
            if ischar(varargin{kk+1})
                switch upper(varargin{kk+1})
                    case 'ON'
                        show_legend = true;
                    case 'OFF'
                        show_legend = false;
                    otherwise
                        fprintf('''%s'' is not a property value of "show_legend".\n',...
                            varargin{kk+1})
                        return
                end
            else
                show_legend = varargin{kk+1};
            end
        case {'GRID','SHOW_GRID'}
            if ~ischar(varargin{kk+1})
                if varargin{kk+1}
                    show_grid = 'on';
                else
                    show_grid = 'off';
                end
            end
        case {'AXIS_LIM','AXIS_LIMIT','LIM','LIMIT','AXIS'}
            if ischar(varargin{kk+1})
                switch upper(varargin{kk+1})
                    case {'AUTO','CALC'}
                        axis_lim = 'AUTO';
                    case {'DEFAULT','DFLT'}
                        axis_lim = 'DEFAULT';
                    otherwise
                        fprintf('''%s'' is not a property value of "axis_lim".\n',...
                            varargin{kk+1})
                        return
                end
            else
                axis_lim = varargin{kk+1};
            end
        otherwise
            fprintf('''%s'' is not a property name.\n',varargin{kk})
            return
    end
end
%% propiedades de bolas
if isempty(balls)
    %balls = [0.5 0.5 0; 0.5 0.5 0; 0.5 0.5 0;0.5 0.5 0;0.5 0.5 0;0.5 0.5 0;0.5 0.5 0];
    radio=1;
    for i=1:n
        balls(i,1) = radio;
        balls(i,2) = radio;
        balls(1,3) = 0;
    end
end
N = size(balls,1);
if isempty(ball_pos)
    ball_pos = zeros(N,3);
    rmax = max(balls(:,1))+0.5;
    lx = dflt_axis_lim(1,1)+1+rmax:2*rmax:dflt_axis_lim(1,2)-1-rmax;
    ly = dflt_axis_lim(2,1)+1+rmax:2*rmax:dflt_axis_lim(2,2)-1-rmax;
    lz = dflt_axis_lim(3,1)+1+rmax:2*rmax:dflt_axis_lim(3,2)-1-rmax;
    [mx,my,mz] = meshgrid(lx,ly,lz);
    for kk = 1:N
        ball_pos(kk,:) = [mx(kk),my(kk),mz(kk)];
    end
end
%% CREACION DE DATOS PARA PROYECTO 
if isempty(ball_v)
    ball_v = 2.5*rand(N,3);
    arr_modulos = modulo_vel(ball_v);
    value = calculo_vel_media(arr_modulos);
    [vect_uniq,num_rep] = contar_repetidos(arr_modulos);
    array = horzcat(vect_uniq, num_rep);
    arr_velocidades = vector_velocidades(vect_uniq, value);
    arr_num_moles = vector_moles(num_rep);
%% IMPRESION DE DATOS RELEVANTES
    disp("Modulos de velocidades");
    disp(arr_modulos);
    disp("Valor de velocidad media");
    disp(value);
    disp("Modulos repetidos");
    disp(array);
    disp("Vector de velocidades");
    disp(arr_velocidades);
    disp("Vector de numero de moles");
    disp(arr_num_moles);
%% GRAFICO DE FUNCION DE DISTRIBUCION
    f = fit(arr_velocidades,arr_num_moles,'gauss1');
    plot(f,arr_velocidades,arr_num_moles),grid;
    xlabel('Velocidad');
    ylabel('Numero de Moleculas');   
end
if isempty(ball_color)
    ball_color = repmat([0 0.5 1 0.75],N,1);
end
if strcmpi(dflt_layout,'2d')
    ball_v(:,3) = 0;
end
if ~debug_mode && show_legend
    fprintf('The property "show_legend" is only used in debug mode.\n')
end
tlen = tmax/dt;
cldmax = 25000;
if isempty(wall_pln)
    n_pln = length(pln_pnt);
    wall_pln = generate_walls;
else
    n_pln = size(wall_pln,1);
end
n_sph = size(wall_sph,1);
if isempty(pln_color)
    pln_color = repmat([0.85 0.4 0.24 0.2],n_pln,1);
end
if isempty(sph_color)
    sph_color = repmat([0.85 0.4 0.24 0.2],n_sph,1);
end
n_grav = size(GravitySrc,1);
n_clmb = size(CoulombSrc,1);
if ischar(axis_lim)
    switch axis_lim
        case 'AUTO'
            axis_lim = find_axislim;
        case 'DEFAULT'
            axis_lim = dflt_axis_lim;
    end
elseif isvector(axis_lim)
    axis_lim = reshape(axis_lim,2,3)';
end
%% MODIFY
if n_pln ~= 0
    plane_magn = sqrt(sum(wall_pln(:,1:3).^2,2));
    wall_pln = wall_pln./repmat(plane_magn,1,4);
end
ball_pos(:,4) = 1;
ball_pos_new = ball_pos;
ball_v_new = ball_v;
GravityVec = GravityVec/sqrt(GravityVec*GravityVec');
GravityField = g_grav*GravityVec;
GravityPlane = [GravityVec -GravityVec*GroundPoint'];

[~,kc,~] = check_overlap;
if kc > 0
    fprintf(['ERROR: There is an overlap in the current sphere layout.'...
        '\nPlease check ''ball_pos'' and try again.\n\n'])
    return
end
if isempty(w_pln_v)
    w_pln_v = cell(1,n_pln);
    for ikv = 1:n_pln
        w_pln_v{ikv} = {@(x) 0, @(x) 0, @(x) 0};
    end
end
if isempty(w_sph_v)
    w_sph_v = cell(1,n_sph);
    for ikv = 1:n_sph
        w_sph_v{ikv} = {@(x) 0, @(x) 0, @(x) 0};
    end
end
if n_sph ~= 0
    sph_cntr = wall_sph(:,2:4);
end

%% PLOTs & SIMULACION

t_dura = max([min([10^(round(log10(dt/N^2))),0.01]),1e-4]);
h_f = figure('name','Sphere Collider','menubar','none','numbertitle','off',...
    'position',[30 270 600 450],'color',[1 1 1],...
    'closerequestfcn',@my_closereq);
close_requested = 'none';
h_a = axes('parent',h_f);
set(h_a,'box','on','projection','perspective','dataaspectratio',[1 1 1])
set(h_a,'xlim',axis_lim(1,:),'ylim',axis_lim(2,:),'zlim',axis_lim(3,:))
set(h_a,'xlimmode','manual','ylimmode','manual','zlimmode','manual')
set(h_a,'xgrid',show_grid,'ygrid',show_grid,'zgrid',show_grid)
% xlabel(h_a,'x')
% ylabel(h_a,'y')
% zlabel(h_a,'z')
rotate3d(h_a)
%Generar Balls
sph_n = 16;
[bx,by,bz] = sphere(sph_n);
xyzdata = cell(1,N);
sf_h = zeros(1,N);
%Inicial xyz_data of balls
for ikb = 1:N
    xyzdata{ikb} = zeros(sph_n+1,sph_n+1,3);
    xyzdata{ikb}(:,:,1) = balls(ikb)*bx;
    xyzdata{ikb}(:,:,2) = balls(ikb)*by;
    xyzdata{ikb}(:,:,3) = balls(ikb)*bz;
    sf_h(ikb) = surface(xyzdata{ikb}(:,:,1),xyzdata{ikb}(:,:,2),xyzdata{ikb}(:,:,3),...
        'FaceColor',ball_color(ikb,1:3),'FaceAlpha',ball_color(ikb,4),...
        'linestyle','none','visible','off','parent',h_a);
end %Inicial xyz_data: balls
for ikb = 1:N
    set(sf_h(ikb),'xdata',xyzdata{ikb}(:,:,1)+ball_pos(ikb,1),...
        'ydata',xyzdata{ikb}(:,:,2)+ball_pos(ikb,2),...
        'zdata',xyzdata{ikb}(:,:,3)+ball_pos(ikb,3),'visible','on')
end
h_pln = zeros(1,n_pln);
for ikw = 1:n_pln
    h_pln(ikw) = patch('xdata',pln_pnt{ikw}(:,1),'ydata',pln_pnt{ikw}(:,2),'zdata',pln_pnt{ikw}(:,3),...
        'parent',h_a,'facecolor',pln_color(ikw,1:3),'facealpha',pln_color(ikw,4),'linestyle','none');
end
h_sph = zeros(1,n_sph);
for ikw = 1:n_sph
    h_sph(ikw) = surface(wall_sph(ikw,1)*bx+wall_sph(ikw,2),...
        wall_sph(ikw,1)*by+wall_sph(ikw,3),...
        wall_sph(ikw,1)*bz+wall_sph(ikw,4),...
        'parent',h_a,'facecolor',sph_color(ikw,1:3),'facealpha',sph_color(ikw,4),'linestyle','none');
end

%% SIMULACION
dt_pre = dt;
dt_c = 0;
t = dt;
cnt_t = 1;
cnt_c = 0;
cnt_error = 0;
while cnt_c < cldmax && cnt_t < tlen
    cnt_c = cnt_c+1;
    ball_pos_new(:,1:3) = ball_pos(:,1:3)+dt_pre*ball_v;
    [dt_int,kc,~] = check_overlap;
    if dt_int < inf
        dt_c = dt_c+dt_int;
        dt_pre = dt-dt_c;
        ball_pos(:,1:3) = ball_pos(:,1:3)+dt_int*ball_v;
        ball_v(kc,:) = ball_v_new(kc,:);
    else
        cnt_t = cnt_t+1;
        t = t+dt;
        for kb1 = 1:N
            a = [0 0 0];
            for kb2 = [1:kb1-1,kb1+1:N]
                norm_balls = ball_pos(kb2,1:3)-ball_pos(kb1,1:3);
                dist_balls = sqrt(norm_balls*norm_balls');
                a_grav_magn = G_newton*balls(kb1,2)*balls(kb2,2)/dist_balls^2;
                a_coulomb_magn = -k_coulomb*balls(kb1,3)*balls(kb2,3)/dist_balls^2;
                norm_balls = norm_balls/dist_balls;
                a = a+(a_grav_magn*grav_ball+a_coulomb_magn)*norm_balls;
            end
            for ksg = 1:n_grav
                norm_src = ball_pos(kb1,1:3)-GravitySrc(ksg,2:4);
                dist_src = sqrt(norm_src*norm_src');
                a_grav_magn = G_newton*balls(kb1,2)*GravitySrc(ksg,1)/dist_src^2;
                norm_src = norm_src/dist_src;
                a = a+a_grav_magn*norm_src;
            end
            for ksc = 1:n_clmb
                norm_src = ball_pos(kb1,1:3)-CoulombSrc(ksc,2:4);
                dist_src = sqrt(norm_src*norm_src');
                a_coulomb_magn = k_coulomb*balls(kb1,3)*CoulombSrc(ksc,1)/dist_src^2;
                norm_src = norm_src/dist_src;
                a = a+a_coulomb_magn*norm_src;
            end
            a = a+GravityField;
            a_frict = -friction_factor/balls(kb1,2)*ball_v(kb1,:);
            ball_v(kb1,:) = ball_v(kb1,:)+dt*(a+a_frict);
            set(sf_h(kb1),'xdata',xyzdata{kb1}(:,:,1)+ball_pos_new(kb1,1),...
                'ydata',xyzdata{kb1}(:,:,2)+ball_pos_new(kb1,2),...
                'zdata',xyzdata{kb1}(:,:,3)+ball_pos_new(kb1,3))
        end
        ball_pos = ball_pos_new;
        for kwv = 1:n_pln
            vw = [w_pln_v{kwv}{1}(t),w_pln_v{kwv}{2}(t),w_pln_v{kwv}{3}(t)];
            pln_pnt{kwv} = pln_pnt{kwv}+([1;1;1;1]*vw)*dt;
            set(h_pln(kwv),...
                'xdata',pln_pnt{kwv}(:,1),...
                'ydata',pln_pnt{kwv}(:,2),...
                'zdata',pln_pnt{kwv}(:,3));
        end
        for kwv = 1:n_sph
            vw = [w_sph_v{kwv}{1}(t),w_sph_v{kwv}{2}(t),w_sph_v{kwv}{3}(t)];
            sph_cntr(kwv,:) = sph_cntr(kwv,:)+vw*dt;
            set(h_sph(kwv),...
                'xdata',wall_sph(kwv,1)*bx+sph_cntr(kwv,1),...
                'ydata',wall_sph(kwv,1)*by+sph_cntr(kwv,2),...
                'zdata',wall_sph(kwv,1)*bz+sph_cntr(kwv,3));
        end
        pause(t_dura);
        dt_pre = dt;
        dt_c = 0;
    end
end

%% REPORTE
if strcmp(close_requested,'abort')
    delete(h_f)
end
fprintf('Numero de ERRORES: %d\n',cnt_error)
close_requested = 'done';

%% FUNCTIONS
    function [dt_int,kc,kc_w] = check_overlap
        dt_int = inf;
        kc = -1;
        kc_w = -1;
        for fk1 = 1:N
            for fkw = 1:n_pln
                dist_plane = abs(wall_pln(fkw,:)*ball_pos_new(fk1,:)');
                iscross = (wall_pln(fkw,:)*ball_pos_new(fk1,:)')*...
                    (wall_pln(fkw,:)*ball_pos(fk1,:)') < 0;
                if dist_plane < balls(fk1,1) || iscross %!!
                    vbn = ball_v(fk1,:)*wall_pln(fkw,1:3)'*wall_pln(fkw,1:3);
                    vwn = ([w_pln_v{fkw}{1}(t),w_pln_v{fkw}{2}(t),w_pln_v{fkw}{3}(t)]*...
                        wall_pln(fkw,1:3)')*wall_pln(fkw,1:3);
                    dist_plane_old = abs(wall_pln(fkw,:)*ball_pos(fk1,:)');
                    dt_int_c = (dist_plane_old-balls(fk1,1))/...
                        sqrt((vbn-vwn)*(vbn-vwn)');
                    if dt_int_c < dt_int
                        dt_int = dt_int_c;
                        ball_v_new(fk1,:) = ball_v(fk1,:)-(1+Cr)*vbn+2*vwn;
                        kc = fk1;
                        kc_w = fkw;
                    end
                end
            end
            for fkw = 1:n_sph
                sph_norm = ball_pos_new(fk1,1:3)-wall_sph(fkw,2:4);
                sph_norm_old = ball_pos(fk1,1:3)-wall_sph(fkw,2:4);
                dist_sphere = sqrt(sph_norm*sph_norm');
                dist_sphere_old = sqrt(sph_norm_old*sph_norm_old');
                iscross = ((dist_sphere-wall_sph(fkw,1))*...
                    (dist_sphere_old-wall_sph(fkw,1))) < 0;
                if abs(dist_sphere-wall_sph(fkw,1)) < balls(fk1,1) || iscross
                    sph_norm = sph_norm/dist_sphere;
                    vbn = (ball_v(fk1,:)*sph_norm')*sph_norm;
                    vwn = ([w_sph_v{fkw}{1}(t),w_sph_v{fkw}{2}(t),w_sph_v{fkw}{3}(t)]*...
                        sph_norm')*sph_norm;
                    dt_int_c = (abs(dist_sphere_old-wall_sph(fkw,1))...
                        -balls(fk1,1))/sqrt((vbn-vwn)*(vbn-vwn)');
                    if dt_int_c < dt_int
                        dt_int = dt_int_c;
                        ball_v_new(fk1,:) = ball_v(fk1,:)-(1+Cr)*vbn+2*vwn;
                        kc = fk1;
                        kc_w = fkw+n_pln;
                    end
                end
            end
            for fk2 = fk1+1:N
                
                norm_balls = ball_pos_new(fk2,1:3)-ball_pos_new(fk1,1:3);
                if norm_balls*norm_balls' < (balls(fk1,1)+balls(fk2,1))^2
                    v2_1 = ball_v(fk2,:)-ball_v(fk1,:);
                    r2_1 = ball_pos(fk2,1:3)-ball_pos(fk1,1:3);
                    d_sq = (det([r2_1([1,2]);v2_1([1,2])])^2+det([r2_1([2,3]);v2_1([2,3])])^2+...
                        det([r2_1([3,1]);v2_1([3,1])])^2)/(v2_1*v2_1');
                    d_c = sqrt(r2_1*r2_1'-d_sq)-sqrt((balls(fk1,1)+balls(fk2,1))^2-d_sq);
                    dt_int_c = abs(d_c)/sqrt(v2_1*v2_1');
                    if d_sq >= (balls(fk1,1)+balls(fk2,1))^2
                        set(sf_h([fk1,fk2]),'facecolor',[0.75 0.25 0.75])
                        disp('ERROR: "d" >= R1+R2')
                        cnt_error = cnt_error+1;
                    end
                    if dt_int_c < 0
                        set(sf_h([fk1,fk2]),'facecolor',[0.75 0.25 0.75])
                        disp('ERROR: "dt_int_c" < 0')
                        cnt_error = cnt_error+1;
                    end
                    ball_pos_new(fk1,1:3) = ball_pos(fk1,1:3)+dt_int_c*ball_v(fk1,1:3);
                    ball_pos_new(fk2,1:3) = ball_pos(fk2,1:3)+dt_int_c*ball_v(fk2,1:3);
                    norm_balls_c = (ball_pos_new(fk2,1:3)-ball_pos_new(fk1,1:3))/...
                        (balls(fk1,1)+balls(fk2,1));
                    if dt_int_c < dt_int
                        dt_int = dt_int_c;
                        vn1 = (ball_v(fk1,:)*norm_balls_c')*norm_balls_c;
                        vn2 = (ball_v(fk2,:)*norm_balls_c')*norm_balls_c;
                        va = (balls(fk1,2)*vn1+balls(fk2,2)*vn2+balls(fk2,2)*...
                            Cr*(vn2-vn1))/(balls(fk1,2)+balls(fk2,2));
                        vb = (balls(fk1,2)*vn1+balls(fk2,2)*vn2+balls(fk1,2)*...
                            Cr*(vn1-vn2))/(balls(fk1,2)+balls(fk2,2));
                        ball_v_new(fk1,:) = ball_v(fk1,:)-vn1+va;
                        ball_v_new(fk2,:) = ball_v(fk2,:)-vn2+vb;
                        kc = [fk1 fk2];
                    end
                end
            end
        end
    end

    function axislim = find_axislim
        pln_idx = combination_recs(1,n_pln,3);
        pt = zeros(size(pln_idx,1),3);
        for fkw = 1:size(pln_idx,1)
            pt(fkw,:) = -wall_pln(pln_idx(fkw,:),1:3)\wall_pln(pln_idx(fkw,:),4);
        end
        pt(isinf(pt) | isnan(pt)) = 0;
        axislim = [min(pt(:,1)),max(pt(:,1));min(pt(:,2)),max(pt(:,2));...
            min(pt(:,3)),max(pt(:,3))];
    end

    function c = combination_recs(m,n,r)
        if (n-m+1) < r
            c = [];
            return;
        end
        if r == 1
            c = (m:n)';
            return;
        end
        nCr = factorial(n-m+1)/(factorial(r)*factorial(n-m+1-r));
        c = zeros(nCr,r);
        cstart = 1;
        for fk = m:(n-r+1)
            csub = combination_recs((fk+1),n,(r-1));
            cend = cstart+size(csub,1)-1;
            c(cstart:cend,1) = fk;
            c(cstart:cend,2:r) = csub;
            cstart = cend+1;
        end
    end

    function my_closereq(varargin)
        if strcmp(close_requested,'none')
            selection = questdlg('Abort Simulation?',...
                'Close Figure',...
                'Yes','No','Yes');
            switch selection
                case 'Yes'
                    close_requested = 'abort';
                    cnt_c = cldmax;
                case 'No'
                    return
            end
        else
            delete(h_f)
        end
    end

    function wall_pln = generate_walls
        wall_pln = zeros(n_pln,4);
        for fkw = 1:n_pln
            vec1 = pln_pnt{fkw}(2,:)-pln_pnt{fkw}(1,:);
            vec2 = pln_pnt{fkw}(3,:)-pln_pnt{fkw}(1,:);
            vecn = cross(vec1,vec2);
            vecn = vecn/sqrt(vecn*vecn');
            wall_pln(fkw,:) = [vecn -vecn*pln_pnt{fkw}(1,:)'];
        end
    end
%% FUNCIONES PROPIAS
    function[arr_modulos] = modulo_vel(ball_v)
        arr_modulos = zeros(N,1);
        for fila=1:N
            arr_modulos(fila,1) = round(sqrt(ball_v(fila,1)^2) + round(ball_v(fila,2)^2) + round(ball_v(fila,3)^2));
        end
    end
    
    function [vect_uniq, num_rep] = contar_repetidos(arr_modulos)
        vect_uniq = unique(arr_modulos);
        num_rep = zeros(length(vect_uniq),1);
        
        for uniq=1:length(vect_uniq)
            for mod=1:length(arr_modulos)
               if vect_uniq(uniq) == arr_modulos(mod)
                   num_rep(uniq) = num_rep(uniq) +1;
               end
            end
        end    
    end

    function[value] = calculo_vel_media(arr_modulos)
        sumatoria = 0;
        for index=1:length(arr_modulos)
            sumatoria = sumatoria + arr_modulos(index,1)^2;
        end
        value = round(sqrt(sumatoria/n));
    end

    function[arr_velocidades] = vector_velocidades(vect_uniq, value)
        arr_velocidades = zeros(length(vect_uniq)+1,1);
        arr_velocidades(1) = value;
        
        for index=2:length(arr_velocidades)
            arr_velocidades(index) = vect_uniq(index-1);
        end    
    end    

    function[arr_num_moles] = vector_moles(num_rep)
        arr_num_moles = zeros(length(num_rep)+1,1);
        arr_num_moles(1) = n;
        
        for index=2:length(arr_num_moles)
            arr_num_moles(index) = num_rep(index-1);
        end    
    end
end


