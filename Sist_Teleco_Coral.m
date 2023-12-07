clear 
clc
load ("Practica_Sist_Tec_Teleco_2324.mat")
plot(bt(:,1),bt(:,2),'o')
hold on
plot(xp(:,1),xp(:,2),'x')
title ('Distribución de estaciones base')
grid on


C1 = 0; %usuarios a los que cubre
C2 = 0; %Coste total
Radius = 1.75;
Matrix = nan * ones(100,100);

for i = 1:length(bt(:,1))
    contador = 1;
    for j = 1:length(xp(:,1))
        %Si el usuario está cubierto, es el usuario j
        distance = sqrt((bt(i,1)-xp(j,1)).^2+(bt(i,2)-xp(j,2)).^2);
        if (distance<=Radius)
            Matrix(i,contador) = j;
            contador = contador+1;
        end
    end
end

Posicion_antenas = randperm(100,25)';


%% Calculamos C1 %%



%% Calculamos C2 %%

%g(x) = Alpha Gentecubre(x) + Betta 1/Costetotal(x) 