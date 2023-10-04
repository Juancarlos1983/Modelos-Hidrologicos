% =====================================================================
% ---------------------   JUAN CARLOS TICONA  -------------------------
% ---------- INSTITUTO DE PESQUISAS HIDRAULICAS (IPH) UFRGS  ----------
% -------------------------- OUTUBRO DE 2023 --------------------------    
% --------------------------------------------------------------------- 
% =====================================================================

function [Q,QO] = TANK4hidro(x)
% =======================================================================
% O modelo de tanque foi desenvolvido originalmente para uso em solos 
% constantemente saturados, no Japão (Sugawara, 1979, 1995). Possui 4
% tanques e 12 parametros.

% Descrição e unidades dos parametros:
% HI1, HI2, HI3 e HI4 altura de cada armazenamento inicial *[mm]*
% HA1, HA2, HB1 e HC1 altura de cada um dos orifícios laterais *[mm]*
% a1, a2, b1, c1, d1 seus respectivos coeficientes de escoamento *[dia-1]*
% a0, b0 e c0 coeficiente de infiltração para os tanques *[dia-1]*

%     HI1 (armazenamento inicial do tanque 1)                    : X1
%     HI2 (armazenamento inicial do tanque 2)                    : X2
%     HI3 (armazenamento inicial do tanque 3)                    : X3
%     HI4 (armazenamento inicial do tanque 4)                    : X4
%     HA1 (altura do orificio a1 do tanque 1)                    : X5
%     HA2 (altura do orificio a2 do tanque 1)                    : X6
%     HB1 (altura do orificio do tanque 2)                       : X7
%     HC1 (altura do orificio do tanque 3)                       : X8
%     a1  (coeficiente de escoamento 1, do tanque 1)             : X9
%     a2  (coeficiente de escoamento 2, do tanque 1)             : X10
%     b1  (coeficiente de escoamento do tanque 2)                : X11
%     c1 (coeficiente de escoamento do tanque 3                  : X12
%     d1 (coeficiente de escoamento do tanque 4)                 : X13
%     a0 (coeficiente de infiltração do tanque 1 para o tanque 2): X14
%     b0 (coeficiente de infiltração do tanque 2 para o tanque 3): X15
%     C0 (coeficiente de infiltração do tanque 3 para o tanque 4): X16
% % X = [HI1, HI2, HI3, HI4, HA1, HA2, HB1, HC1, a1, a2, b1, c1, d1, a0, b0, c0]
%% =======================================================================
% Começa a ler os dados de entrada
QO = textread('vaz_goias_v.txt','%f')';                % vazão observado em m3/s
P = textread('prec_goias_v.txt','%f');          % Precitacão em mm/dia
ETR  = textread('evap_goias_v.txt','%f');  % Evapotranspiração mm/dia
NT=length(QO);
PAR = x; % Le os parametros em uma coluna

% Areas das bacias hidrograficas aplicadas
% Bacia Ijui
% Area = 5414; % km^2
% Area = 86.4/Area;  % conversão vazão em unidades de m^3/s
% % Bacia Canoas
% Area = 989; % km^2
% Area = 86.4/Area;  % conversão vazão em unidades de m^3/s
% % Bacia Goias
Area = 1817; % km^2
Area = 86.4/Area;  % conversão vazão em unidades de m^3/s

    for I = 1:NT

        X = PAR;

        %% CALCULOS DO PRIMEIRO DIA
        if I == 1
            %% Calcula armazenamento no Tanque 1 e se há evapotranspiração no 2o tanque
            if X(1) + P(I) - ETR(I) <0
                S1(I) = 0;
                ETR2(I) = abs(X(1) + P(I) - ETR(I));
            else
                S1(I) = X(1) + P(I) - ETR(I);
                ETR2(I) = 0;
            end

            %% Calcula infiltração do 1o para o 2o tanque
            qi1(I) = X(14)*S1(I);

            %% Calcula armazenamento no Tanque 1 e se há evapotranspiração no 2o tanque
            if X(2) + qi1(I) - ETR2(I) < 0
                S2(I) = 0;
                ETR3(I) = abs(X(2) + qi1(I) - ETR2(I));
            else
                S2(I) = X(2) + qi1(I) - ETR2(I);
                ETR3(I) = 0;
            end

            %% Calcula infiltração do 2o para o 3o tanque
            qi2(I) = X(15)*S2(I);

            %% Calcula armazenamento no Tanque 3 e se há evapotranspiração no 4o tanque
            if X(3) + qi2(I) - ETR3(I) < 0
                S3(I) = 0;
                ETR4(I) = abs(X(3) + qi2(I) - ETR3(I));
            else
                S3(I) = X(3) + qi2(I) - ETR3(I);
                ETR4(I) = 0;
            end

            %% Calcula infiltração do 3o para o 4o tanque
            qi3(I) = X(16)*S3(I);

            %% Calcula armazenamento do primeiro dia tanque 4
            if X(4) + qi3(I) - ETR4(I) < 0
                S4(I) = 0;
            else
                S4(I) = X(4) + qi3(I) - ETR4(I);
            end

            %% Calcula escoamento do primeiro orificio do 1o tanque
            if  S1(I) > X(6)
                qs1(I) = X(10)*(S1(I) - X(6));
            else
                qs1(I) = 0;
            end
            %% Calcula escoamento do segundo orificio do 1o tanque
            if  S1(I) > X(5)
                qs2(I) = X(9)*(S1(I) - X(5));
            else
                qs2(I) = 0;
            end
            %% Calculo do escoamento dos orificios do 1o tanque
            qst1(I) = qs1(I) + qs2(I);

            %% Calcula o escoamento do orificio do 2o tanque
            if  S2(I) > X(7)
                qs3(I) = X(11)*(S2(I) - X(7));
            else
                qs3(I) = 0;
            end
            qst2(I) = qs3(I);

            %% Calcula o escoamento do orificio do 3o tanque
            if S3(I) > X(8)
                qs4(I) = X(12)*(S3(I) - X(8));
            else
                qs4(I) = 0;
            end

            %% Calcula o escoamento do orificio do 4o tanque
            qs5(I) = X(13)*(S4(I));

            %% Calculo do armazenamento final dos tanques
            S1f(I) = S1(I) - qst1(I) - qi1(I);
            S2f(I) = S2(I) - qst2(I) - qi2;
            S3f(I) = S3(I) - qs4(I);
            S4f(I) = S4(I) - qs5(I);

            %% Escoamento total que chega ao exutório
            Q(I) = qs1(I) + qs2(I) + qs3(I) + qs4(I) + qs5(I);
        end

        %% CALCULOS DOS OUTROS DIAS
        if I>1
            %% Calcula armazenamento do tanque1 e se há evapotranspiração no 2o tanque
            if (S1f(I - 1) + P(I) - ETR(I)) < 0
                S1(I) = 0;
                ETR2(I) = abs(S1f(I - 1) + P(I) - ETR(I));
            else
                S1(I) = S1f(I - 1) + P(I) - ETR(I);
                ETR2(I) = 0;
            end

            %% Calcula infiltração do 1o para o 2o tanque
            qi1(I) = X(14)*S1(I);

            %% Calcula armazenamento do tanque2 e se há evapotranspiração no 3o tanque
            if (S2f(I - 1) + qi1(I) - ETR2(I)) < 0
                S2(I) = 0;
                ETR3(I) = abs(S2f(I - 1) + qi1(I) - ETR2(I));
            else
                S2(I) = S2f(I - 1) + qi1(I) - ETR2(I);
                ETR3(I) = 0;
            end

            %% Calcula infiltração do 2o para o 3o tanque
            qi2(I) = X(15)*S2(I);

            %% Calcula armazenamento do tanque3 e se há evapotranspiração no 4o tanque
            if S3f(I - 1) + qi2(I) - ETR3(I) < 0
                S3(I) = 0;
                ETR4(I) = abs(abs(S3f(I - 1) + qi2(I) - ETR3(I))) ;
            else
                S3(I) = S3f(I - 1) + qi2(I) - ETR3(I);
                ETR4(I) = 0;
            end

            %% Calcula infiltração do 3o para o 4o tanque
            qi3(I) = X(16)*S3(I);

            %% Calcula armazenamento dos demais dias tanque4
            if S4f(I - 1) + qi3(I) - ETR4(I) < 0
                S4(I) = 0;
            else
                S4(I) = S4f(I - 1) + qi3(I) - ETR4(I) ;
            end

            %% Calcula escoamento do primeiro orificio do 1o tanque
            if  S1(I) > X(6)
                qs1(I) = X(10)*(S1(I) - X(6));
            else
                qs1(I) = 0;
            end
            %% calcula escoamento do segundo orificio do 1o tanque
            if  S1(I) > X(5)
                qs2(I) = X(9)*(S1(I) - X(5));
            else
                qs2(I) = 0;
            end
            %% Calculo do escoamento dos orificios do 1o tanque
            qst1(I) = qs1(I) + qs2(I);

            %% calcula o escoamento do orificio do 2o tanque
            if  S2(I) > X(7)
                qs3(I) = X(11)*(S2(I) - X(7));
            else
                qs3(I) = 0;
            end
            qst2(I) = qs3(I);

            %% calcula o escoamento do orificio do 3o tanque
            if S3(I) > X(8)
                qs4(I) = X(12)*(S3(I) - X(8)); 
            else
                qs4(I) = 0;
            end

            %% Calcula o escoamento do orificio do 4o tanque
            qs5(I) = X(13)*S4(I);

            %% Calculo do armazenamento final dos tanques
            S1f(I) = S1(I) - qst1(I) - qi1(I);
            S2f(I) = S2(I) - qst2(I) - qi2(I);
            S3f(I) = S3(I) - qs4(I);
            S4f(I) = S4(I) - qs5(I);

            %% Escoamento total que chega ao exutório
            Q(I) = qs1(I) + qs2(I) + qs3(I) + qs4(I) + qs5(I);
        end
    end
    Q = Q/Area;
end