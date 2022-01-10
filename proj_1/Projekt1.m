clear;
a=6378137;
e2=0.00669437999013;
macierzDane= load('dane2.txt');
phiB = 41.297445;
lambdaB = 2.0832941;
hB = 4;

phi=macierzDane(:,1); % wszystkie el z kolumny nr1
lambda=macierzDane(:,2);% wszystkie el z kolumny nr2
h=macierzDane(:,3);% wszystkie el z kolumny nr3

% współrzędne samolotu

N_lotnisko = (a./sqrt(1-e2 .* sind(phi) .* sind(phi)));
x = (N_lotnisko+h).*cosd(phi).*cosd(lambda);
y = (N_lotnisko+h).*cosd(phi).*sind(lambda);
z = (N_lotnisko.*(1-e2)+h).*sind(phi);

% lotniska współrzędne
N_lotnisko = (a./sqrt(1-e2*sind(phiB)));
XB = (N_lotnisko+hB).*cosd(phiB).*cosd(lambdaB);
YB = (N_lotnisko+hB).*cosd(phiB).*sind(lambdaB);
ZB = (N_lotnisko.*(1-e2)+hB).*sind(phiB);

n =  ((-sind(phi)) .* cosd(lambda)) .* (x - XB) 
+(-sind(phi) .* sind(lambda)) .* (y - YB) 
+ cosd(phi) .* (z - ZB);

e = (-sind(lambda)) .* (x - XB) + cosd(lambda .* y - YB);

u = (cosd(phi) .* cosd(lambda)) .* (x - XB) 
+ (cosd(phi) .* sind(lambda)) .* (y - YB)
+ sind(phi) .* (z - ZB);

% tworzymy macierz 
xyz = [x y z];
neu = [n e u];

s = sqrt(n.*n + e.*e + u.*u);

pomocniczyZenit = u ./ sqrt(n.*n + e.*e + u.*u);
zenit = acosd(pomocniczyZenit);

pomocniczyAzymut = e ./ n;
azymut = atand(pomocniczyAzymut);



geoscatter(phi,lambda,5,'ro');

plot3(x,y,z);

plot3(n,e,h);

for i = 1:size(neu, 1)
    if ((n(i, 1) < 0) && (e(i, 1) > 0)) || ((n(i, 1) < 0) && (e(i, 1) < 0))
    azymut(i, 1) = azymut(i, 1) + 180
    elseif ((n(i, 1) > 0) && (e(i, 1) < 0))
     azymut(i, 1) = azymut(i, 1) + 360
    end

end



    