format short g

% v0 - brzina aviona (const) i pocetna brzina projektila u m/s
% v_eksp - brzina sirenja oblaka eksplozije u m/s
v0=250;
v_eksp=350;

% r - precnik zaokreta aviona u m
r=1000;

% g - koeficijent sile zemljine teze

% t0 - pocetni trenutak (pocetak zaokreta)
% t1 - trenutak izbacivanja projektila
% t2 - trenutak zavrsavanja zaokreta (predjeni put (l)/brzina (v0))
% t_eksp - trenutak eksplozije (pad projektila)
% t_kraj - trenutak sretanja oblaka eksplozije sa avionom
t0=0;
t2=r*pi/v0;

%%%%%%%%%%%%%%%%%%%%%%%
%%% KRETANJE AVIONA %%%
%%%%%%%%%%%%%%%%%%%%%%%

% u intervalu [t0,t2], avion se krece u zaokretu
teta = @(t) pi*t/t2;
x_1 = @(t) r*sin(teta(t));
y_1 = @(t) 11000 - r*cos(teta(t));

% u intervalu [t2, inf] avion se krece po pravoj liniji
x_2 = @(t) -v0*(t-t2);
y_2 = @(t) 12000 + t*0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KRETANJE PROJEKTILA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pozicija projektila zavisice od trenutka t, ali i od trenutka
% ispaljivanja projektila, t1, jer on odredjuje ugao pod kojim
% je projektil ispaljen
x_proj = @(t, t1) v0*cos(teta(t1))*t + r*sin(teta(t1));
y_proj = @(t, t1) -g*t.^2/2 + v0*sin(teta(t1))*t + 11000 - r*cos(teta(t1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODREDJIVANJE t_eksp I t_susreta %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vektor razmatranih trenutaka ispaljivanja
t1_array=[];
% vektor odgovarajucih trenutaka eksplozije
t_eksp_array=[];
% vektor odgovarajucih x koordinata pada eksploziva
x_eksp_array=[];
% vektor odgovarajucih trenutaka susreta udarnog talasa sa avionom
t_susreta_array=[];
% vektor odgovarajucih poluprecnika udarnog talasa 
r_eksplozije_array=[];

% ispitujemo samo trenutke od 0 do t2/2, jer inace avion krece u zaokret
% unazad i ispaljuje projektil u nepovoljnom smeru. Racunacemo za svako t1
% sa korakom 0.001, smatramo da je to dovoljno precizno.
for t1=0:0.001:6
    % ISCRTAVANJE LEPIH GRAFIKA
    % Za dobijanje raznih grafika kao onih iz rada, potrebno je izmeniti 
    % opseg t1 i otkomentarisati naredne redove i definisati times_proj
    % ili ga prebaciti iznad ovoga odozgo.
    % ylim([0,15000]);
    % plot(x_proj(times_proj, t1), y_proj(times_proj,t1), 'DisplayName', ['t1=' num2str(t1)]);
    % hold on
    
    % nulu funkcije nalazimo funkcijom fzero u okolini 100, eksperimentalno
    % utvrdjeno da je to okej za sve slucajeve
    t_eksp=fzero(@(t) y_proj(t,t1), 100);
    
    % lako odredjujemo x koordinatu projektila kada imamo t_eksp
    x_eksp=x_proj(t_eksp, t1);
    
    % cuvamo obe vrednosti u vektore
    t_eksp_array=[t_eksp_array t_eksp];
    x_eksp_array=[x_eksp_array x_eksp];

    % t_susreta eksplozije sa avionom desava se u trenutku kada je 
    % poluprecnik lopte, R(t) = v_eksp*(t-t_eksp), jednak distanci aviona
    % od mesta pada projektila.
    % t-t_eksp je vreme proteklo od eksplozije
    % trazicemo nulu sledece funkcije:
    f= @(t) v_eksp*(t-t_eksp) - sqrt((-x_2(t)+x_eksp).^2 + 12000.^2);
    t_susreta = fzero(f, 300);
    t_susreta_array = [t_susreta_array t_susreta];
    
    % na kraju je jos potrebno odrediti i poluprecnih udanog talasa
    % u trenutku t_susreta
    r_eksplozije = sqrt((-x_2(t_susreta)+x_eksp).^2 + 12000.^2);
    r_eksplozije_array = [r_eksplozije_array r_eksplozije];
    
    % dodajemo i razmatrano t1 u vektor
    t1_array = [t1_array t1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODREDJIVANJE OPTIMALNOG t1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r_max, i_max] = max(r_eksplozije_array);
disp('Najveci moguci poluprecnik lopte eksplozije je:');
disp(r_max);
disp('i to se desava ako je trenutak ispustanja t1=');
t1_opt=t1_array(i_max);
disp(t1_opt);

disp('x_eksp=');
disp(x_eksp_array(i_max));
disp('t_eksp=');
disp(t_eksp_array(i_max));
disp('t_susreta=');
disp(t_susreta_array(i_max));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ISCRTAVANJE FINALNOG STANJA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ISCRTAVANJE KRETANJA AVIONA
% razmatramo vreme podeljeno na dva dela, zbog dva razlicita tipa
% kretanja aviona
% potrebno je otkomentarisati sve linije ispod

 times_1 = linspace(t0,t2,100);
 times_2 = linspace(t2, t_susreta_array(i_max), 100);
 times_plane = [times_1 times_2];
 times_proj = linspace(0, t_eksp_array(i_max), 200);

 x = [x_1(times_1) x_2(times_2)];
 y = [y_1(times_1) y_2(times_2)];
 plot(x,y,'DisplayName','avion')
 hold on

% ISCRTAVANJE KRETANJA NAJBOLJEG PROJEKTILA
plot(x_proj(times_proj, t1_opt), y_proj(times_proj,t1_opt), 'DisplayName', 'projektil');
hold on
 
 legend


