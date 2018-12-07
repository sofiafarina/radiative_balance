clear % close all; 
clc
% close all chiude tutte le figure aperte 
code = "RBL";

% *** input parameters ***

ztoa=50; % TOA height in km
wp1=1; % mixing ratio profile type
H1=5; % and scale height (km) for IR-active gas
wp2=3; % mixing ratio profile type
H2=5; % and scale height (km) for SW-active gas
d1t=1; % integral from ztoa to surface of density of gas-1 (IR)
d2t=1; % integral from ztoa to surface of density of gas-2 (SW)
k1a=0; % IR absorption cross section per unit mass
k2a=0; % SW absorption cross section per unit mass
icl=0; % icl=0 if no cloud layer is present and 1 if it's present 
nube=[0,0]; % =[base, top] cloud extends from base (km) to top (km)
kp1=0; % IR cloud absorption coefficient
kp2=0; % SW cloud absorption coefficient


% parameters varying from computation to computation

% default variables 
nlay = 101; 
wp1 = 1; H1 = 5; k1a = 0; d1t = 1;
wp2 = 1; H2 = 5; k2a = 0; d2t = 1;
icl = 0; nube = [0,0]; kp1 = 1; kp2 = 1; 

%nlay = 6; k1a = 5; k2a = 2;  %JOB 0 
nlay = 1; %JOB 1
nlay = 2; ztoa = 50; k1a = 0.00001; %JOB 2 
nlay = 2; k1a = 1; % JOB 3 
nlay = 6; % JOB 4
nlay = 21; % JOB 5 


% basic variables common to all computations

ap=0.3;             % Planetary albedo
sun=(1-ap)*1370/4;  % Averaged solar irradiance at TOA
mudif=3/5;          % mu* value to compute diffuse tranmittance
sigma=5.6704e-8;    % [W/(m^2k^4)] Stefan­Boltzmann constant

% create vector profiles from scalar input quantities 

k1(1:nlay) = k1a;
k2(1:nlay) = k2a;


% scaling of input quantities (SI units)

ztoa=ztoa*1d3; % ztoa [m]
H1=H1*1d3; % H1 [m]
H2=H2*1d3; % H2 [m]
nube = nube*1.e3; 

% define the geometry and layering 
% gli strati sono tutti uguali (non è il modo più sensato per scriverli, si
% dovrebbe mantenere costante lo spessore ottico o la quantità di materia)

if (nlay == 1)
    z = [0]; % altezza
    dz = 0; % spessore
    zm = ztoa; 
else
    dzs = ztoa/(nlay-1);
    dz(1:nlay-1) = dzs; % tutti gli strati da 1 a nlay-1 hanno spessore dzs
    dz(nlay) = 0; % questo invece ha spessore 0
    z = [ztoa:-dzs:dzs 0]; % level heights [m], partendo da ztoa si toglie dzs fino ad arrivare a dzs e poi l'ultimo strato è a zero
    zm = [ztoa - dzs*0.5:-dzs:dzs*0.5 0]; % mean layer heights [m]
end

% compute also zk and zmk, same as zm and z but in KM
zk = z/1000;
zmk = zm/1000;

% define the atmospheric density profile using the scale height parameter H

do = 1.225;
H = 101325/(9.8*do); % air density scale height (m)
d = do*exp(-z/H); % density at level z 

% define the mixing ratio profile for gas-1 (active in the IR) usign wp1
% ad esempio il secondo è un modo per descrivere il rapporto di
% mescolanza del vapor acqueo che decresce con la quota in modo exp
if (wp1 == 1)
    w1(1:nlay) = 1;  % constant mixing ratio (like CO2)
elseif (wp1 == 2)
    w1 = exp(-z/H1); % exponantially decreasing mixing ratio (like H2O)
end

% define mixing ratio profile for gas-2 (active in SW)

if (wp2 == 1)
      w2(1:nlay)=1;  % constant mix.ratio (like CO2)
   elseif (wp2 == 2)
      w2=exp(-z/H2); % exponentially decreasing mix.ratio (like H2O)
   elseif (wp2 == 3)   % constant mix.ratio for a number of layers (like O3)
      zkb = 20; zkt = 50;
      lay = find(zk > zkb & zk < zkt);
      zo = (zkt+zkb)*0.5;
      sig = (zkt - zkb)/10;
      w2(1:nlay) = 1.e-10;
      for k = lay
          w2(k) = exp(-(zk(k)-zo)^2./(2*sig^2));
      end
      % gaussiano centrato intorno ad un certo intervallo kinda rapporto di
      % mescolanza
      %figure, plot(w2,zk)
end

% computation of the density of absorber profile and normalization of the
% integral from ztoa to atmosphere 

if (wp2==1)
      w2(1:nlay)=1;  % constant mix.ratio (like CO2)
   elseif (wp2==2)
      w2=exp(-z/H2); % exponentially decreasing mix.ratio (like H2O)
   elseif (wp2==3)   % constant mix.ratio for a number of layers (like O3)
      w2(1:nlay)=0.0001; 
      lay=find(zk>20 & zk<50); 
      w2(lay)=1;
end

% computation of density of the absorber profile and normalization 

da1 = d.*w1;         % density of absorber at each level (ricorda che d e w1 sono vettori, senza il punto ci viene fuori una matrice!)
da2 = d.*w2;         % density of absorber at each level

% compute integral of dax and normalise
tot1 = 0; 
tot2 = 0;
ix = 1:nlay-1;


% calcolo dell'integrale con il metodo dei trapezi 

for j = ix
    tot1 = tot1+(da1(j) + da1(j+1))*dz(j)*0.5;
    tot2 = tot2+(da2(j) + da2(j+1))*dz(j)*0.5;
end

% scale da using integral values da1 and da2:
da1 = da1*d1t/tot1;
da2 = da2*d2t/tot2;


% compute the optical depth profile vector and layer transmissivities

% computation of the optical dpet of each atm layer in IR and SW 
% delta Chi = k_1 rho_a \Deltaz 

chir(1:nlay) = 0; % IR layer (dove dz è un vettore)  
chsw(1:nlay) = 0; % SW layer 

for j = 1:nlay-1
    chir(j) = dz(j)*0.5*(k1(j)*da1(j) + k1(j+1)*da1(j+1))/mudif;
    chsw(j) = dz(j)*0.5*(k2(j)*da2(j) + k2(j+1)*da2(j+1))/mudif;
end 

% supponiamo che ci sia una nube -> dobbiamo calcolare lo spessore ottico
% della nube e poi aggiungerlo a quello del gas 
% si usano dei coefficienti al posto delle sezioni d'urto in unità di massa
% quello che cambia sono le unità di misura 
% \Delta chi_n = k_p*\Deltaz dove k_p = sezionedurto*densità

if icl == 1 % se c'è la nube
    index = find(zm >= nube(1) & zm <= nube(2)); % così ho l'indice degli strati dove c'è la nube 
    chir(index) = chir(index) + kp1*dz(index)/mudif; 
    chsw(index) = chsw(index) + kp2*dz(index)/mudif; 
end

figure, semilogx(chir, zk, chsw, zk), title("chir")

% si vede bene se si trasforma il grafico in semilogaritmico  
    
% compute the integrated total OD from TOA to each level (tchir e tchsw)
tchir(1:nlay) = 0.;  % settiamo il primo a zero 
tchsw(1:nlay) = 0.; 
tchir(2:nlay-1) = cumsum(chir(2:nlay-1)); % lui ci mette reverse 
tchsw(2:nlay-1) = cumsum(chsw(2:nlay-1));
% attenzione agli indici, deve essere fino a nlay-1 perchè non dobbiamo
% sommare anche quello della terra che è zero

figure, plot(tchir, zk, tchsw, zk), title("tchir")

% compute the layer transmissivity profile vector in the IR and SW
% remember to set the trasmissivity of the surface layer to zero (perchè è
% quella di un corpo nero la sua trasmissività)
% trasmissività dello strato = exp(-\Delta Chi)


tir = exp(-chir);
tir(nlay) = 0; 
tsw = exp(-chsw); 
tsw(nlay) = 0;

figure, plot(tir, zk, tsw, zk), title("tir")

% computation of layer absorptivity and emissivity

air = 1 - tir; 
eir = air; 
asw = 1 - tsw; 

% compute the total transmittance matrix TRAL 

tral(1:nlay, 1:nlay) = 1; 
for j = 1:nlay-2
    for k = j+2:nlay
        tral(j,k) = tral(j,k-1)*tir(k-1);
        tral(k,j) = tral(j,k); 
    end 
end


ttsw(1:nlay) = 1.; 
for j = 2:nlay
    ttsw(j) = ttsw(j-1)*tsw(j-1);
end 

% compute matrix M 
% M = E * TRAL * A 

M = eir.*tral.*air' 

for i = 1:nlay-1
    M(i,i) = eir(i)*(-2)
end

M(nlay, nlay) =  eir(nlay)*(-1);


% compute vector S 
Sv = -sun.*ttsw.*asw; 

% computation of system solution 

Yv=linsolve(M,Sv');
Tv=(Yv/sigma).^(0.25); % mean layer temperature at equilibrium 
Tv 
























    
    
    
    
    
    
    
    
    
