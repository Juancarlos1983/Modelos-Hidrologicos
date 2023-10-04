% =======================================================================
% -----------------------------------------------------------------------
% -----------------------   JUAN CARLOS TICONA  -------------------------
% ----------- INSTITUTO DE PESQUISAS HIDRAULICAS (IPH) UFRGS  -----------
% ---------------------------- HyMOD MOdel ------------------------------
% -------------------------- Dezembro do 2023 ---------------------------
% -----------------------------------------------------------------------
% =======================================================================

function [Q,QO] = HYMOD(X)
% =======================================================================
% O modelo HyMOD combina uma rotina de umidade do solo do tipo PDM (exemplo, 
% Moore (2007)) com uma cascata Nash de três reservatórios lineares que 
% simulam fluxo rápido e um único reservatório linear destinado a simular 
% fluxo lento (Wagener et al., 2001; Boyle, 2001). Possui 5 parâmetros e 5 reservatórios.

% Descrição e unidades dos parametros:
% HI1, HI2, HI3 e HI4 altura de cada armazenamento inicial *[mm]*
% HA1, HA2, HB1 e HC1 altura de cada um dos orifícios laterais *[mm]*
% a1, a2, b1, c1, d1 seus respectivos coeficientes de escoamento *[dia-1]*
% a0, b0 e c0 coeficiente de infiltração para os tanques *[dia-1]*
% Parametros fixos:
%     Sm0 (armazenamento inicial do reservatório 1)               : X1
%     F10 (armazenamento inicial do 1o reservatório rapido)       : X2
%     F20 (armazenamento inicial do 2o reservatório rapido)       : X3
%     F30 (armazenamento inicial do 3o reservatório rapido)       : X4
%     SL0 (armazenamento inicial do reservatório lento)            : X5
% Parametros calibraveis:
%     Cmax (Armazenamento máximo de umidade do solo)              : X1
%     b (Parâmetro de forma da curva da área contribuinte)        : X2
%     a (Fração da precipitação efetiva que é de fluxo rápido)    : X3
%     kf (Coeficiente de escoamento dos reservatórios rápidos)    : X4
%     ks (Coeficiente de escoamento do reservatório lento)        : X5
%     X = [Cmax, b, a, kf, ks]

%% =======================================================================
% Começa a ler os dados de entrada
QO = textread('vaz_goias_c.txt','%f')';         % vazão observado em m3/s
P = textread('prec_goias_c.txt','%f');          % Precitacão em mm/dia
ETR  = textread('evap_goias_c.txt','%f');       % Evapotranspiração mm/dia
NT   = length(QO);
% Bacia Ijui
% Area = 5414; % km^2
% Area = 86.4/Area;  % conversão vazão em unidades de m^3/s
% % Bacia Canoas
% Area = 989; % km^2
% Area = 86.4/Area;  % conversão vazão em unidades de m^3/s
% % Bacia Goias
Area = 1817; % km^2
Area = 86.4/Area;  % conversão vazão em unidades de m^3/s

%%Inicializar reservatorios do modelo
% load ('storeinitial_hymod_goias.prn')   % comeca a ler os dados de entrada, armazenamento inicial dos reservatorios
% S(1)  = storeinitial_hymod_goias(1);
% SL(1) = storeinitial_hymod_goias(2);
% F1(1) = storeinitial_hymod_goias(3);
% F2(1) = storeinitial_hymod_goias(4);
% F3(1) = storeinitial_hymod_goias(5);
S(1)  = X(1);
SL(1) = X(2);
F1(1) = X(3);
F2(1) = X(4);
F3(1) = X(5);
Cmax = X(6);
b    = X(7);
a    = X(8);
kf   = X(9);
ks   = X(10);


    for t = 1:NT-1
            %% Calculo da pecipitação efetiva
	    So = S(t);
            c1 = Cmax*(1 - (1 -(b + 1)*S(t)/Cmax)^(1/(b + 1)));

	    er1 = max(P(t) - Cmax + c1, 0);

            P(t) = P(t) - er1;

	    d = min((c1 + P(t))/Cmax, 0);
            S(t) = (Cmax/(b + 1))*(1 - (1 - d)^(b + 1));

            er2 = max(P(t) - (S(t) - So), 0);

            evap = (1 - ( ( Cmax/(b + 1) - S(t) )/( Cmax/(b + 1) )))*Ep(t);
            S(t + 1) = max(S(t) - evap, 0);

            Pe = er1 + er2;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Calcula fluxo do Tanque 1 para o reservatório lento
            Ps = (1 - a)*Pe;

            %% Calcula armazenamento no reservatório lento 
            if SL(t) + Ps <0
                SL(t + 1) = 0;
            else
                SL(t + 1) = (1 - ks)*(SL(t) + Ps);
            end

            %% Calcula fluxo do reservatório lento
            Qs(t + 1) = (ks/(1 - ks))*SL(t + 1);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Calcula fluxo do Tanque 1 para o 1o reservatório rápido
            Pf = (a)*Pe;
            
            %% Calcula fluxo que sai do 1o reservatório rápido
            Ff1 = kf*F1(t);
            %% Calcula armazenamento no 1o reservatório rápido 
            if F1(t) - Ff1 + Pf <0
                F1(t + 1) = 0;
            else
                F1(t + 1) = F1(t) - Ff1 + Pf;
            end

            %% Calcula fluxo que sai do 2o reservatório rápido
            Ff2 = kf*F2(t);
            %% Calcula armazenamento no 2o reservatório rápido 
            if F2(t) - Ff2 + Ff1 <0
                F2(t + 1) = 0;
            else
                F2(t + 1) = F2(t) - Ff2 + Ff1;
            end

            %% Calcula fluxo que sai do 3o reservatório rápido
            Ff3 = kf*F3(t);

            %% Calcula armazenamento no 3o reservatório rápido 
            if F3(t) -Ff3 + Ff2<0
                F3(t + 1) = 0;
            else
                F3(t + 1) = F3(t) -Ff3 + Ff2;
            end

            %% Escoamento total que chega ao exutório
            Q(t + 1) = Qs(t + 1) + Ff3;
    end
    Q = Q/Area;
end