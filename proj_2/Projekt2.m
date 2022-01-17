clear;

% gwiazda Spica z gwiazdozbioru Panny
rektascensja = 13 + 25/60 + 11.579/3600; %%13h 25m 12s  
deklinacja = -(11 + 9/60 + 40.75/3600); %%-11° 9′ 41″ 

%  Norylsk  
phi = 69.33629228312138;
lambda = 88.18223614146471;

%% Quito
%  phi = -0.17972737242763911;
% lambda = -78.46556813496137; 

%% Porto Alegre 
%phi = -30.035112361544094;
%lambda = -51.214916894632374;

h = transpose([0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]);

% Wartości kątowe dla każdej godziny
kot_godzinowy = katgodz(2021, 12, 7, h, lambda, rektascensja);

% obliczenia kąta dla każdej z warunkiem 
for i = 1:24
    if kot_godzinowy(i) > 360
        kot_godzinowy(i) = kot_godzinowy(i) - 360;
    end
end

% odleglość zenitalna dla każdej godziny
cos_z = sind(phi).*sind(deklinacja) + cosd(phi).*cosd(deklinacja).*cosd(kot_godzinowy);
z = acosd(cos_z);
wysokosc = 90 - z;

% Azymut gwiazdy dla każdej godziny
licznik = -cosd(deklinacja).*sind(kot_godzinowy);
mianownik = cosd(phi).*sind(deklinacja) - sind(phi).*cosd(deklinacja).*cosd(kot_godzinowy);
Az = atan2d(licznik, mianownik);
for i=1:24
   if Az(i) < 0 
       Az(i) = Az(i) + 360; 
   elseif Az(i) > 360
       Az(i) = Az(i) - 360;
   end
end 

% Obliczenie wsp do układu prostokątnego
x = 1.*sind(z).*cosd(Az);
y = 1.*sind(z).*sind(Az);
z = 1.*cosd(z);

% Wizualizacja półkuli 
%[X,Y,Z] = sphere(24);
%X = X(13:end,:);
%Y = Y(13:end,:);
%Z = Z(13:end,:);
%surf(X,Y,Z,'FaceColor','green','FaceAlpha',0.5)
%axis equal, 
%hold on;
%% Wizualizacja gwiazdy
%scatter3(x,y,z, 500, 'red', '.')

plot(h, wysokosc)
%plot(h, z)


function [t] = katgodz(y, m, d, h, lambda, alfa)
    jd = juliandate(datetime(y, m, d)); % dni
    g = GMST(jd); % stopnie
    UT1 = h * 1.002737909350795; % godziny

    % obliczenie czasu gwiazdowego (w stopniach)
    S = UT1*15 + lambda + g; 
    % obliczenie kąta godzinowego (w stopniach)
    t = S - alfa*15;
end

function g = GMST(jd)
    T = (jd - 2451545) / 36525;
    g = 280.46061837 + 360.98564736629  * (jd - 2451545.0) + 0.000387933*T.^2 - T.^3/38710000;
    g = mod(g, 360);

end