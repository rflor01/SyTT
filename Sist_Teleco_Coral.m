clear 
close all
clc

load ("Practica_Sist_Tec_Teleco_2324.mat")

figure(1)
plot(bt(:,1),bt(:,2), 'o', 'Color','red')
hold on
plot(xp(:,1),xp(:,2), 'x', 'Color','blue')
viscircles(bt,1.5)
xlabel('Distancia [km]')
ylabel('Distancia [km]')
title('Distribución de Estaciones Base y Usuarios')

C1 = 0; %usuarios a los que cubre
C2 = 0; %Coste total
Radius = 1.75;
Matrix = nan * ones(100,100);

tic
for i = 1:length(bt(:,1))
    contador = 1;
    for j = 1:length(xp(:,1))
        %Si el usuario está cubierto, es el usuario j
        distance = norm(bt(i,:)-xp(j,:));
        if (distance<=Radius)
            Matrix(i,contador) = j;
            contador = contador+1;
        end
    end
end
run_time_norm = toc;

C1 = 0; %usuarios a los que cubre
C2 = 0; %Coste total
Radius = 1.75;
Matrix = nan * ones(100,100);

tic
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
run_time_dist = toc;
Posicion_antenas = randperm(100,25)';

N_BTS = 25;
temp = [ones(1,N_BTS) zeros(1,height(bt)-N_BTS)];
rand_temp = temp(randperm(height(bt),height(bt)));
ones_pos = find(rand_temp == 1);




%% Calculamos C1 %%



%% Calculamos C2 %%

%g(x) = Alpha Gentecubre(x) + Betta 1/Costetotal(x) 