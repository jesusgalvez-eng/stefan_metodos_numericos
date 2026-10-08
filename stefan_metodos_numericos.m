%% METODOS NUMERICOS PARA LA ECUACION DEL CALOR Y EL PROBLEMA DE STEFAN
% Autor: Jesus

clear; close all; clc;

guardarFiguras = true;
carpetaSalida = 'resultados_stefan';
if guardarFiguras && ~exist(carpetaSalida, 'dir')
    mkdir(carpetaSalida);
end

%% PARTE A. ECUACION DEL CALOR EN UN DOMINIO FIJO
% Modelo: u_t = alpha*u_xx, 0 < x < L
% Bordes: u(0,t) = u(L,t) = 0
% Condicion inicial con dos modos de Fourier.

L = 1.0;
alpha = 0.10;
Nx = 50;                         % Numero de intervalos espaciales
x = linspace(0, L, Nx + 1)';
dx = L/Nx;
tFinal = 0.50;

% FTCS necesita r = alpha*dt/dx^2 <= 1/2.
rObjetivo = 0.45;
dtTentativo = rObjetivo*dx^2/alpha;
Nt = ceil(tFinal/dtTentativo);
dt = tFinal/Nt;
r = alpha*dt/dx^2;

u0 = @(x) sin(pi*x/L) + 0.40*sin(3*pi*x/L);
q0 = @(x) (pi/L)*cos(pi*x/L) + 0.40*(3*pi/L)*cos(3*pi*x/L);
uExacta = @(x,t) sin(pi*x/L).*exp(-alpha*(pi/L)^2*t) ...
    + 0.40*sin(3*pi*x/L).*exp(-alpha*(3*pi/L)^2*t);

fprintf('\nPARTE A: ECUACION DEL CALOR\n');
fprintf('Nx = %d, dx = %.5f, Nt = %d, dt = %.6f, r = %.4f\n', ...
    Nx, dx, Nt, dt, r);

if r > 0.5
    warning('FTCS puede ser inestable porque r > 0.5.');
end

[Uftcs, errorFTCS, tiempo] = resolverCalor( ...
    'FTCS', x, dt, Nt, alpha, u0, uExacta);
[Ubtcs, errorBTCS, ~] = resolverCalor( ...
    'BTCS', x, dt, Nt, alpha, u0, uExacta);
[Ucn, errorCN, ~] = resolverCalor( ...
    'CN', x, dt, Nt, alpha, u0, uExacta);
[Ukeller, Qkeller, errorKeller, ~] = resolverKellerBox( ...
    x, dt, Nt, alpha, u0, q0, uExacta); 

uFinExacta = uExacta(x, tFinal);
erroresFinales = [errorFTCS(end); errorBTCS(end); ...
    errorCN(end); errorKeller(end)];
metodos = {'FTCS'; 'BTCS'; 'Crank-Nicolson'; 'Keller-Box'};
tablaErrores = table(metodos, erroresFinales, ...
    'VariableNames', {'Metodo', 'Error_L2_final'});
disp(tablaErrores);

% Grafico 1: perfiles al tiempo final.
fig1 = figure('Color', 'w');
plot(x, uFinExacta, 'k-', 'LineWidth', 2.2); hold on;
plot(x, Uftcs(:,end), '--', 'LineWidth', 1.5);
plot(x, Ubtcs(:,end), '-.', 'LineWidth', 1.5);
plot(x, Ucn(:,end), ':', 'LineWidth', 2.0);
plot(x, Ukeller(:,end), '--', 'LineWidth', 1.3);
grid on; box on;
xlabel('Posicion x'); ylabel('Temperatura adimensional u(x,t)');
title(sprintf('Comparacion al tiempo t = %.2f', tFinal));
legend('Solucion exacta', 'FTCS', 'BTCS', 'Crank-Nicolson', ...
    'Keller-Box', ...
    'Location', 'best');
if guardarFiguras
    print(fig1, fullfile(carpetaSalida, '01_perfiles_calor.png'), ...
        '-dpng', '-r200');
end

% Grafico 2: evolucion del error L2.
fig2 = figure('Color', 'w');
semilogy(tiempo, max(errorFTCS, eps), 'LineWidth', 1.6); hold on;
semilogy(tiempo, max(errorBTCS, eps), 'LineWidth', 1.6);
semilogy(tiempo, max(errorCN, eps), 'LineWidth', 1.8);
semilogy(tiempo, max(errorKeller, eps), 'LineWidth', 1.5);
grid on; box on;
xlabel('Tiempo'); ylabel('Error L_2 discreto');
title('Evolucion del error respecto de la solucion exacta');
legend('FTCS', 'BTCS', 'Crank-Nicolson', 'Keller-Box', ...
    'Location', 'best');
if guardarFiguras
    print(fig2, fullfile(carpetaSalida, '02_error_calor.png'), ...
        '-dpng', '-r200');
end

%% PRUEBA DE ESTABILIDAD DE FTCS
% Se agrega un modo de frecuencia alta muy pequeno. Si r > 1/2, ese modo
% crece y hace visible la inestabilidad del metodo explicito.

numeroPasos = 150;
[maxEstable, pasoEstable] = pruebaEstabilidadFTCS( ...
    0.45, Nx, L, alpha, numeroPasos);
[maxInestable, pasoInestable] = pruebaEstabilidadFTCS( ...
    0.60, Nx, L, alpha, numeroPasos);

fig3 = figure('Color', 'w');
semilogy(pasoEstable, maxEstable, 'LineWidth', 1.7); hold on;
semilogy(pasoInestable, maxInestable, 'LineWidth', 1.7);
grid on; box on;
xlabel('Numero de paso temporal');
ylabel('max |u|');
title('Prueba de estabilidad de FTCS');
legend('r = 0.45 (estable)', 'r = 0.60 (inestable)', ...
    'Location', 'best');
if guardarFiguras
    print(fig3, fullfile(carpetaSalida, '03_estabilidad_ftcs.png'), ...
        '-dpng', '-r200');
end

%% PARTE B. PROBLEMA DE STEFAN REDUCIDO A UNA EDO
% Se considera una capa liquida 0 < x < s(t), calentada en x = 0.
% T(0,t) = Th y T(s(t),t) = Tm.
%
% Aproximacion elegida:
%   T(x,t) = Th - (Th-Tm)*x/s(t)
%
% Un balance de energia conduce a
%   ds/dt = k*(Th-Tm)/(rho*Lefectivo*s),
% donde Lefectivo = Lf + cp*(Th-Tm)/2 si se incluye el calor sensible
% promedio de la capa. Si se desea la aproximacion cuasi-estacionaria mas
% simple, se puede cambiar usarCalorSensible a false.

k = 0.60;                % Conductividad termica [W/(m K)]
rho = 1000;              % Densidad [kg/m^3]
Lf = 334000;             % Calor latente [J/kg]
cp = 4180;               % Calor especifico [J/(kg K)]
Th = 20;                 % Temperatura en el borde caliente [deg C]
Tm = 0;                  % Temperatura de fusion [deg C]
deltaT = Th - Tm;
usarCalorSensible = true;

if usarCalorSensible
    Lefectivo = Lf + cp*deltaT/2;
    etiquetaModelo = 'ode45 con calor sensible';
else
    Lefectivo = Lf;
    etiquetaModelo = 'ode45 cuasi-estacionario';
end

A = k*deltaT/(rho*Lefectivo);       % Constante con unidades m^2/s
s0 = 1.0e-4;                        % Interfaz inicial [m]
tFinalStefan = 6*3600;              % Seis horas [s]
edoStefan = @(t,s) A/max(s, 1.0e-12);
opciones = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
[tStefan, sNumerica] = ode45(edoStefan, [0 tFinalStefan], s0, opciones);
sAnalitica = sqrt(s0^2 + 2*A*tStefan);
errorFrontera = max(abs(sNumerica - sAnalitica));

fprintf('\nPARTE B: PROBLEMA DE STEFAN REDUCIDO\n');
fprintf('A = %.6e m^2/s\n', A);
fprintf('s(%.2f h) = %.4f mm\n', tFinalStefan/3600, 1000*sNumerica(end));
fprintf('Error maximo ode45 vs. formula = %.3e m\n', errorFrontera);

% Comparacion opcional con el modelo que solo toma calor latente.
Acuasi = k*deltaT/(rho*Lf);
sCuasi = sqrt(s0^2 + 2*Acuasi*tStefan);

fig4 = figure('Color', 'w');
plot(tStefan/3600, 1000*sNumerica, 'LineWidth', 2.0); hold on;
plot(tStefan/3600, 1000*sAnalitica, '--', 'LineWidth', 1.5);
plot(tStefan/3600, 1000*sCuasi, ':', 'LineWidth', 1.8);
grid on; box on;
xlabel('Tiempo [h]'); ylabel('Posicion de la frontera s(t) [mm]');
title('Movimiento de la frontera en el modelo reducido de Stefan');
legend(etiquetaModelo, 'Formula del modelo reducido', ...
    'Aproximacion cuasi-estacionaria', 'Location', 'best');
if guardarFiguras
    print(fig4, fullfile(carpetaSalida, '04_frontera_stefan.png'), ...
        '-dpng', '-r200');
end

% Perfiles de temperatura aproximados dentro de la zona liquida.
horasPerfil = [0.5, 2, 4, 6];
tiemposPerfil = horasPerfil*3600;
sPerfil = interp1(tStefan, sNumerica, tiemposPerfil, 'pchip');
xTemperatura = linspace(0, 1.05*max(sPerfil), 250);

fig5 = figure('Color', 'w'); hold on;
for j = 1:numel(tiemposPerfil)
    temperatura = nan(size(xTemperatura));
    dentro = xTemperatura <= sPerfil(j);
    temperatura(dentro) = Th - deltaT*xTemperatura(dentro)/sPerfil(j);
    plot(1000*xTemperatura, temperatura, 'LineWidth', 1.7, ...
        'DisplayName', sprintf('t = %.1f h', horasPerfil(j)));
end
grid on; box on;
xlabel('Posicion x [mm]'); ylabel('Temperatura [deg C]');
title('Perfiles lineales usados en la reduccion a una EDO');
legend('Location', 'best');
if guardarFiguras
    print(fig5, fullfile(carpetaSalida, '05_perfiles_stefan.png'), ...
        '-dpng', '-r200');
end

%% PARTE C. CONDICION DE FRONTERA EXPONENCIAL
% Se reemplaza la temperatura constante del borde por
%   T(0,t) = Tm + deltaT*(1-exp(-t/tau)).
% En la aproximacion cuasi-estacionaria de una fase se obtiene
%   ds/dt = C*(1-exp(-t/tau))/s,
% donde C = k*deltaT/(rho*Lf). La integracion directa conduce a
%   s(t)^2 = s0^2 + 2*C*(t-tau*(1-exp(-t/tau))).

tau = 1*3600;                         % Tiempo caracteristico [s]
Cexp = k*deltaT/(rho*Lf);             % Constante [m^2/s]
edoExponencial = @(t,s) Cexp*(1-exp(-t/tau))/max(s, 1.0e-12);
[tExpNum, sExpNum] = ode45(edoExponencial, ...
    [0 tFinalStefan], s0, opciones);
sExpExacta = sqrt(s0^2 + 2*Cexp*(tExpNum ...
    - tau*(1-exp(-tExpNum/tau))));
errorExp = max(abs(sExpNum-sExpExacta));

tExp = linspace(0, tFinalStefan, 500)';
Texp = Tm + deltaT*(1-exp(-tExp/tau));
sExp = sqrt(s0^2 + 2*Cexp*(tExp-tau*(1-exp(-tExp/tau))));
sBordeConstante = sqrt(s0^2 + 2*Cexp*tExp);

fprintf('\nPARTE C: BORDE EXPONENCIAL\n');
fprintf('tau = %.2f h\n', tau/3600);
fprintf('s exponencial (%.2f h) = %.4f mm\n', ...
    tFinalStefan/3600, 1000*sExp(end));
fprintf('Error maximo ode45 vs. formula exponencial = %.3e m\n', ...
    errorExp);

fig6 = figure('Color', 'w');
subplot(2,1,1);
plot(tExp/3600, Texp, 'LineWidth', 1.8);
grid on; box on;
ylabel('Temperatura [deg C]');
title('Condicion de frontera exponencial');
subplot(2,1,2);
plot(tExp/3600, 1000*sBordeConstante, '--', 'LineWidth', 1.5); hold on;
plot(tExp/3600, 1000*sExp, 'LineWidth', 1.8);
grid on; box on;
xlabel('Tiempo [h]'); ylabel('Posicion s(t) [mm]');
legend('Borde constante', 'Borde exponencial', 'Location', 'best');
if guardarFiguras
    print(fig6, fullfile(carpetaSalida, '06_frontera_exponencial.png'), ...
        '-dpng', '-r200');
end

if guardarFiguras
    fprintf('\nLas figuras se guardaron en la carpeta: %s\n', carpetaSalida);
end

%% FUNCIONES LOCALES

function [U, errorL2, tiempo] = resolverCalor( ...
    metodo, x, dt, Nt, alpha, u0, uExacta)

%RESOLVERCALOR Resuelve u_t = alpha*u_xx con bordes homogeneos.

    Nx = numel(x) - 1;
    dx = x(2) - x(1);
    r = alpha*dt/dx^2;
    tiempo = (0:Nt)'*dt;
    U = zeros(Nx + 1, Nt + 1);
    U(:,1) = u0(x);
    U(1,1) = 0;
    U(end,1) = 0;
    errorL2 = zeros(Nt + 1, 1);
    errorL2(1) = sqrt(dx*sum((U(:,1) - uExacta(x,0)).^2));

    Ni = Nx - 1;
    e = ones(Ni,1);
    D = spdiags([e, -2*e, e], -1:1, Ni, Ni);
    I = speye(Ni);

    switch upper(metodo)
        case 'BTCS'
            Aimp = I - r*D;
        case 'CN'
            Acn = I - 0.5*r*D;
            Bcn = I + 0.5*r*D;
    end

    for n = 1:Nt
        u = U(:,n);
        unueva = zeros(size(u));

        switch upper(metodo)
            case 'FTCS'
                unueva(2:end-1) = u(2:end-1) + r*( ...
                    u(3:end) - 2*u(2:end-1) + u(1:end-2));
            case 'BTCS'
                unueva(2:end-1) = Aimp\u(2:end-1);
            case 'CN'
                unueva(2:end-1) = Acn\(Bcn*u(2:end-1));
            otherwise
                error('Metodo desconocido: %s', metodo);
        end

        U(:,n+1) = unueva;
        exacta = uExacta(x, tiempo(n+1));
        errorL2(n+1) = sqrt(dx*sum((unueva - exacta).^2));
    end
end

function [historialMaximo, pasos] = pruebaEstabilidadFTCS( ...
    r, Nx, L, alpha, numeroPasos)
%PRUEBAESTABILIDADFTCS Muestra el efecto de usar r menor o mayor que 1/2.

    x = linspace(0, L, Nx + 1)';
    dx = L/Nx;
    dt = r*dx^2/alpha; %#ok<NASGU>
    u = sin(pi*x/L) + 1e-8*sin((Nx - 1)*pi*x/L);
    u(1) = 0;
    u(end) = 0;
    historialMaximo = zeros(numeroPasos + 1, 1);
    historialMaximo(1) = max(abs(u));

    for n = 1:numeroPasos
        unueva = zeros(size(u));
        unueva(2:end-1) = u(2:end-1) + r*( ...
            u(3:end) - 2*u(2:end-1) + u(1:end-2));
        u = unueva;
        historialMaximo(n+1) = max(abs(u));
    end
    pasos = (0:numeroPasos)';
end

function [U, Q, errorL2, tiempo] = resolverKellerBox( ...
    x, dt, Nt, alpha, u0, q0, uExacta)
%RESOLVERKELLERBOX Aplica cajas centradas a u_x=q y q_x=u_t/alpha.

    Nx = numel(x)-1;
    Nn = Nx+1;
    dx = x(2)-x(1);
    tiempo = (0:Nt)'*dt;
    U = zeros(Nn, Nt+1);
    Q = zeros(Nn, Nt+1);
    U(:,1) = u0(x);
    Q(:,1) = q0(x);
    U(1,1) = 0;
    U(end,1) = 0;
    errorL2 = zeros(Nt+1,1);
    errorL2(1) = sqrt(dx*sum((U(:,1)-uExacta(x,0)).^2));

    % Dos ecuaciones de borde y dos ecuaciones por cada caja espacial.
    A = spalloc(2*Nn, 2*Nn, 8*Nx+2);
    A(1,1) = 1;
    A(2,Nn) = 1;
    fila = 3;
    for j = 2:Nn
        % Promedio centrado de u_x=q en (x_{j-1/2},t_{n+1/2}).
        A(fila,j) = 1/dx;
        A(fila,j-1) = -1/dx;
        A(fila,Nn+j) = -0.5;
        A(fila,Nn+j-1) = -0.5;
        fila = fila+1;

        % Promedio centrado de q_x=u_t/alpha en la misma caja.
        A(fila,Nn+j) = 1/dx;
        A(fila,Nn+j-1) = -1/dx;
        A(fila,j) = -1/(alpha*dt);
        A(fila,j-1) = -1/(alpha*dt);
        fila = fila+1;
    end

    for n = 1:Nt
        u = U(:,n);
        q = Q(:,n);
        b = zeros(2*Nn,1);
        fila = 3;
        for j = 2:Nn
            b(fila) = -(u(j)-u(j-1))/dx + 0.5*(q(j)+q(j-1));
            fila = fila+1;
            b(fila) = -(q(j)-q(j-1))/dx ...
                - (u(j)+u(j-1))/(alpha*dt);
            fila = fila+1;
        end

        solucion = A\b;
        U(:,n+1) = solucion(1:Nn);
        Q(:,n+1) = solucion(Nn+1:end);
        exacta = uExacta(x, tiempo(n+1));
        errorL2(n+1) = sqrt(dx*sum((U(:,n+1)-exacta).^2));
    end
end
