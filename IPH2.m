% =====================================================================
% ---------------------------------------------------------------------
% ---------------------   JUAN CARLOS TICONA  -------------------------
% ---------- INSTITUTO DE PESQUISAS HIDRAULICAS (IPH) UFRGS  ----------
% ---------------------------- IPH II Model ---------------------------
% -------------------------- OUTUBRO DE 2023 --------------------------    
% --------------------------------------------------------------------- 
% =====================================================================

function [Q,QO] = IPH2(x)
% Modelo hidrológico conceitual: IPH II
%
% Copyright (C) 2023 Juan Carlos Ticona Gutierrez
% O modelo IPH II ([6, 80, 81]) é um modelo conceitual concentrado chuva-vazão 
% que simula a vazão do rio por meio de dados de precipitação e evaporação como 
% dados de entrada e usando quatro rotinas diferentes. O modelo é amplamente 
% utilizado no Brasil para modelar as respostas hidrológicas de pequenas bacias 
% e normalmente funciona em intervalos de tempo diários ou horários.Possui
% 7 parametros detalhada em:
% BRAVO, J. M. et al. Avaliação visual e numérica da calibração do modelo hidrológico 
% IPH II com fins educacionais. In: XVII Simpósio Brasileiro de Recursos Hídricos, 2007,
% São Paulo. Anais do XVII Simpósio Brasileiro de Recursos Hídricos. Porto Alegre: 
% Associação Brasileira de Recursos Hídricos, v. 1. 2007.
% TUCCI, C. E. M.; CLARKE, R. T. Adaptative forecasting with a conceptual rainfall-runoff 
% model. Hydrological forecasting Proceedings of the Oxford Symposium IAHS, n. 129, p. 425-454, 1980.
% TUCCI, C. E. M.; ORDONEZ, J. S.; SIMÕES LO-PES, M. Modelo Matemático Precipi-tação-Vazão
% IPH II Alguns Resultados. AnaisIV Simpósio Brasileiro de Recursos Hídricos. Fortaleza: [s.n.]. 1981.

% Começa a ler os dados de entrada
QO = textread('vaz_goias_v.txt','%f')';              % vazao observado em m3/s
PREC = textread('prec_goias_v.txt','%f');            % precipitacao em mm/dia
EVAP = textread('evap_goias_v.txt','%f');            % precipitacao mm/dia
NT=length(PREC);

% Comeca a ler os dados de entrada inicial
load ('storeinitial_iphii_goias.prn')    
% ------------------------------------------------------------------------
NH  = storeinitial_iphii_goias(1);  % Tempo de concentração (dias)
XN  = storeinitial_iphii_goias(2);  % fator de forma da bacia (n)
%---------- dados condicao inicial e da corrida -------
T1   = storeinitial_iphii_goias(3)';     % Percolacao(m3/s)          
QT1  = storeinitial_iphii_goias(4)';     % Q subterraneo(m3/s)        
QS1  = storeinitial_iphii_goias(5)';     % Q superficial(m3/s)       
AREA = storeinitial_iphii_goias(6)';     % Area bacia(km2)          
AAT  = storeinitial_iphii_goias(7)';     % Delta de tempo(min)
AINP = storeinitial_iphii_goias(8)';     % Area impermeavel(%)
% ------------------------------------------------------------------------
%-------PARAMETROS--------------------------------------------------------
%x=load ('PARAMETROS.prn');  %parametros da corrida 
%PAR=parametros(:); % Armazem dos parametros em uma coluna
PAR = x; % LE OS parametros em uma coluna
%------------------------------------------------------------------------
% CALCULO DO HISTOGRAMA TEMPO-AREA SINTETICO
N=0;
%
C=0.5/((NH/2.0)^XN);% constante "a" do histograma tempo area sintetico
L=NH/2.0;
XA=0.0;
HIST=[];
HIST(1)=C;
for I=1:L
  XA1=C*(I^XN);     % formula do histograma tempo area sintetico até 50% da area da bacia e utilizada para TConc/2
  HIST(I)=XA1-XA;
  XA=XA1;
  L1=NH-I+1;
  HIST(L1)=HIST(I);
end  
%------------------------------------------------------
% TRANSFORMAÇÃO DE UNIDADES ------------
FATOR=AREA*100/AAT/6.0; %fator de transformacao de unidades
TC=NH*AAT/60.0;
MT=NT;
AUXIT1=T1;
AUXIQT1=QT1;
AUXIQS1=QS1;
T1=T1/FATOR;
QT1=QT1/FATOR;
QS1=QS1/FATOR;
%
% CHAMA AO IPH2 PARA CALCULO SERIE CONTINUA
% CALL OBJEC(X,NF,VFO,KM,NPAR,QC,SFNS,SQ,QP,NSF)
% DENTRO DO ANTERIOR OBJECT
X=PAR;    % parametros%
%
R=0.0;
N=0;
BI=X(1)/log(X(3))/(X(1)-X(2));  
AI=-X(1)*BI;
AT=0.0;
BT=-X(1)/X(2)/log(X(3));
AIL=-AI/BI;
BIL=1/BI;
ATL=0.0;
BTL=1/BT;
SMAX=-X(1)/log(X(3));
X(4)=exp(-1/X(4));
X(5)=exp(-1/X(5));
%
for KT=1:NH
   PV(KT)= 0.0;         
end
S=AT+BT*T1;
RI=AIL+BIL*S;
QT=QT1;      
QS=QS1;
ALF=X(7);
VES=0.0;   
%
Q=[];
% CALCULO DE SIMULAÇÃO DE SERIE CONTINUA
for J=1:MT           
   P=PREC(J);
   E=EVAP(J);
   RIB= X(2);
   H= X(3);
   XK1=X(4);
   XK2=X(5);
   RMAX=X(6);
   %----------------------------------------------------------------------
   %   CALL MODEL(PAUX,EAUX,RI,SX,QTT,QST,QMM,X(2),X(3),X(4),X(5),X(6),
   %	Esta subrotina é o modelo IPH2 propriamente dito.	
   %----------------------------------------------------------------------
   bal=P-E;
   if bal < 0.0
      EP=E-P;
      P=0.0;
      res= EP-R;
      if res <= 0.0
         R=R-EP;
      else
         EP=EP-R;
         R=0.0;
         ER=EP*S/SMAX;
         S=S-ER;
         if S < 0.0
            ER=ER+S;
            S=0.0;
         else
            RI=AIL+BIL*S;
            T=ATL+BTL*S;
            ER=ER+P;
         end
      end  
   else
      P=P-E;
      ER=E;
      RD=RMAX-R;
      sob=P-RD;
      if sob <= 0.0
         R=R+P;
         P=0.0;
      else
         P=P-RD;
         R=RMAX;  
      end   
   end
   AT1=1.0;
   par=P;
   cas=P-RI;
   if cas < 0.0 
      CR=(P/RI)^2/((P/RI)+ALF);
      VES=VES+P*CR;
      P=P*(1-CR);
      S1=(S*(2.0-1.0/BT)+2.0*P)/(2.0+1.0/BT);
      RI1=AIL+BIL*S1;
      sac= P-RI1;
      if sac < 0.0
         T=ATL+BTL*S1; 
         VE=0.0;
         VI=P;
      else
         SX=AI+BI*P;
         ATX=2.0*BT*(SX-S)/(2.0*P*BT+2.0*AT-SX-S);
         A=AT1-ATX;
         RAUX=P;
         VAUX=P*ATX;
         RI1=RIB+(RAUX-RIB)*H^AT1;
         S1=AI+BI*RI1;
         T=ATL+BTL*S1;
         VI=RIB*AT1+(RAUX-RIB)*(H^AT1-1.0)/log(H)+VAUX;
         VE=P*AT1-VI+VAUX;
      end
   else
      RAUX=RI;
      VAUX=0.0;
      RI1=RIB+(RAUX-RIB)*H^AT1;
      S1=AI+BI*RI1;
      T=ATL+BTL*S1;
      VI=RIB*AT1+(RAUX-RIB)*(H^AT1-1)/log(H)+VAUX;
      VE=P*AT1-VI+VAUX;
   end          
   VP=S-S1+VI;
   VE= VE*(1-AINP/100.0) + par*AINP/100.0;   
   for KT=1:NH
      PV(KT)=PV(KT)+VE*HIST(KT);
   end
   VE=PV(1);
   LL=NH-1;
   for KT=1:LL
      PV(KT)=PV(KT+1);
   end
   PV(NH)=0.0;
   QS = QS*XK1+VE*(1.-XK1);
   QT = QT*XK2+VP*(1.-XK2);
   Q(J)=QS+QT;
   S=S1;
   RI=RI1;          
   %   RETURN
   %--------------  END da subroutina MODEL  -------------------------
   %--------------  CONTINUA OBJECT  ---------------------------------
   %
   Q(J)=Q(J)*FATOR;   % Q calculado
   SUP(J)=QS*FATOR; % SUPERF
   TER(J)=QT*FATOR; % SUBT
end
T1=AUXIT1;
QT1=AUXIQT1;
QS1=AUXIQS1;