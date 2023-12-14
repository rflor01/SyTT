clear 
close all

%%%%% INICIA VALORES
load ("Practica_Sist_Tec_Teleco_2324.mat")
N_BTS = 25;
num_individuos_iniciales = 80;
Tamanyo_matriz_coral = 100;
%%%%%%%%%%
figure(1)
plot(bt(:,1),bt(:,2), 'o', 'Color','red')
hold on
plot(xp(:,1),xp(:,2), 'x', 'Color','blue')
%viscircles(bt,1.5)
xlabel('Distancia [km]')
ylabel('Distancia [km]')
title('Distribución de Estaciones Base y Usuarios')

Radius = 1.75;

Matriz_coral = nan * ones(Tamanyo_matriz_coral,100);

Personas = obtain_personas(bt,xp,Radius);

poblacion_inicial = init_poblacion(num_individuos_iniciales,bt,N_BTS);


posiciones = randperm(100,num_individuos_iniciales);

for i = 1:num_individuos_iniciales

    Matriz_coral(posiciones(i),:) = poblacion_inicial(i,:);
end
hijos = cruce(100,Matriz_coral,N_BTS);
N_personas = obtain_alcance(Matrix(1,:),Personas,0);
N_personas_unicas = obtain_alcance(Matrix(1,:),Personas,1);
N_personas
N_personas_unicas


%% Calculamos C1 %%



%% Calculamos C2 %%

%g(x) = Alpha Gentecubre(x) + Betta 1/Costetotal(x) 


function [hijos] = cruce(num_hijos, Matrix,num_BTS_Total)
    hijos = zeros(num_hijos,100);
    padres = find(~isnan(Matrix(:,1)));
    num_padres = length(padres);
    for i = 1:num_hijos
        progenitores_indexes = randperm(num_padres,2);
        dos_padres = Matrix(padres(progenitores_indexes),:);
        corte = randperm(99,1);
        hijo = [dos_padres(1,1:corte) dos_padres(2,corte+1:end)];
        num_BTS = sum(hijo);
        diferencia = num_BTS_Total-num_BTS;
        if diferencia > 0
            potenciales = find(hijo == 0);
            nuevas_idx = randperm(length(potenciales),diferencia);
            hijo(potenciales(nuevas_idx)) = 1;
        elseif diferencia < 0
            potenciales = find(hijo == 1);
            viejas_idx = randperm(length(potenciales),abs(diferencia));
            hijo(potenciales(viejas_idx)) = 0;
        end
        hijos(i,:)= hijo;
        %Ya hemos corregido los hijos, falta mutarlos
    end
end
function [Cost] = obtain_cost(muestra,costes)
    Cost = muestra*costes;
end

function [N_personas] = obtain_alcance(muestra,personas,modo)
%En el modo 0 una persona cuenta como 2 si le cubren 2 antenas
%En el modo 1 una persona cuenta como 1 si le cubren 2 antenas
pos=find(muestra==1);
    if modo == 0
        fila = personas(pos,:);
        N_personas = sum(sum(~isnan(fila)));
    else
        recorridos = [];
        for i = 1:length(pos)
            idx_nan = find(isnan(personas(pos(i),:)));
            temp_recorridos = personas(pos(i),1:(idx_nan(1)-1));
            recorridos = unique([recorridos temp_recorridos]);
        end
        N_personas = length(recorridos);

    end
end

function [Matrix] = init_poblacion(num_vectores,bt,N_BTS)
    Matrix = zeros(num_vectores,height(bt));
    
    for i = 1:num_vectores
        temp = zeros(1,height(bt));
        flag = 1;
        while flag == 1
            flag = 0;
            idx_rand = randperm(height(bt),N_BTS);
    
            temp(idx_rand) = 1;
            for j = 1:i
                if(all(temp == Matrix(j,:)))
                    flag = 1;
                end
    
            end
        end
        Matrix(i,:) = temp;
    end

end


function [Personas] = obtain_personas(bt,xp,Radius)
    Personas = nan * ones(100,24);
    for i = 1:length(bt(:,1))
        contador = 1;
        for j = 1:length(xp(:,1))
            %Si el usuario está cubierto, es el usuario j
            distance = sqrt((bt(i,1)-xp(j,1)).^2+(bt(i,2)-xp(j,2)).^2);
            if (distance<=Radius)
                Personas(i,contador) = j;
                contador = contador+1;
            end
        end
    end
end