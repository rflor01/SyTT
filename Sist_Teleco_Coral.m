clear 
close all

rng(0,"threefry")

%%%%% INICIA VALORES
load ("Practica_Sist_Tec_Teleco_2324.mat")
N_BTS = 25;
num_individuos_iniciales = 80;
Tamanyo_matriz_coral = 100;
num_generaciones = 100;
max_vida = 3;
alpha = 1;
betta = 1;
modo = 0;
Radius = 1.75;
%%%%%%%%%%




Matriz_coral = nan * ones(Tamanyo_matriz_coral,100);

Personas = obtain_personas(bt,xp,Radius);

poblacion_inicial = init_poblacion(num_individuos_iniciales,bt,N_BTS);

coste_min = obtain_min_cost(C);
alcance_max = obtain_max_alcance(Personas);



posiciones = randperm(100,num_individuos_iniciales);

for i = 1:num_individuos_iniciales

    Matriz_coral(posiciones(i),:) = poblacion_inicial(i,:);
end
vector_temporal = zeros(1,num_generaciones);
for i = 1:num_generaciones
    max_funcion_obj = 0;
    for j = 1:height(Matriz_coral)
        if ~isnan(Matriz_coral(j,1))
        objetivo_actual = function_objetivo(Matriz_coral(j,:),alcance_max,coste_min, alpha, betta,modo,Personas,C);
            if objetivo_actual > max_funcion_obj
              max_funcion_obj = objetivo_actual;
              mejor = j;
            end
        end

    end
    vector_temporal(i)=max_funcion_obj;
    figure(1)
    plot(1:i,vector_temporal(1:i),':.')
    grid minor
    xlabel('Número de iteraciones')
    ylabel('Resultado de la Función de Coste g(x)')
    %title('Resultado de la Función de Coste en relación a las iteraciones sobre el Algoritmo de Coral')
    hijos = cruce(100,Matriz_coral,N_BTS);
    hijos = mutar_hijos(hijos);
    Matriz_coral = BATALLA(Matriz_coral,hijos,max_vida,alcance_max,coste_min,alpha,betta,modo,Personas,C);
end
max_funcion_obj = 0;
for i = 1:height(Matriz_coral)
    objetivo_actual = function_objetivo(Matriz_coral(i,:),alcance_max,coste_min, alpha, betta,modo,Personas,C);
    if objetivo_actual > max_funcion_obj
        max_funcion_obj = objetivo_actual;
        mejor = i;
    end

end

fprintf("La mejor solución en índice: <strong>%d</strong>\n",mejor);
fprintf("El resultado de la función objetivo es: <strong>%.20f</strong>\n", max_funcion_obj);
bts_usadas = find(Matriz_coral(mejor,:));
bts_no_usadas = find(~Matriz_coral(mejor,:));
bt_sol_idx = sprintf('%d, ', bts_usadas);
bt_sol_idx = bt_sol_idx(1:end-2);
fprintf("Las BT empleadas son: <strong>%s</strong>\n", bt_sol_idx);

figure(2)
plot(bt(bts_usadas,1),bt(bts_usadas,2), 'o', 'Color','red')
hold on
plot(xp(:,1),xp(:,2), 'x', 'Color','blue')
hold on 
viscircles(bt(bts_usadas,:),Radius*ones(25,1));
hold on
plot(bt(bts_no_usadas,1),bt(bts_no_usadas,2), 'o', 'Color','#77AC30')
hold on
xlabel('Distancia [km]')
ylabel('Distancia [km]')
title('Distribución de Estaciones Base y Usuarios')

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

function [hijos_mutados] = mutar_hijos(hijos)
    hijos_mutados = zeros(length(hijos),100);
    for i = 1:length(hijos)
        aleatorio = randperm(100,1);
        hijos_mutados(i,:) = hijos(i,:);
        if aleatorio == 1
            unos = find(hijos(i,:));
            n_zeros = find(~hijos(i,:));
            indices_unos = randperm(25,10);
            indices_zeros = randperm(75,10);
            hijos_mutados(i,unos(indices_unos)) = 0;
            hijos_mutados(i,n_zeros(indices_zeros)) = 1;
        end
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

function [max_personas] = obtain_max_alcance(Personas)
    num_personas_por_estacion = zeros(1,100);
    for i = 1:height(Personas)

        num_personas_por_estacion(i) = length(find(~isnan(Personas(i,:))));

    end
    ordenada = sort(num_personas_por_estacion);
    max_personas = sum(ordenada(76:end));
end
function [min_cost] = obtain_min_cost(C)
    ordenada = sort(C);
    min_cost = sum(ordenada(1:25));
end


function [valor] = function_objetivo(muestra,alcance_max,coste_min, alpha, betta,modo,Personas,coste)
    valor = alpha*(obtain_alcance(muestra,Personas,modo)/alcance_max) + betta*(coste_min/obtain_cost(muestra,coste));
end




function [Matriz_coral] = BATALLA(Matriz_coral, hijos,max_vida, alcance_max, coste_min, alpha, betta,modo,Personas,coste)
    
    for i = 1:height(hijos)
        vidas = max_vida;
        potencia_hijo = function_objetivo(hijos(i,:),alcance_max,coste_min,alpha,betta,modo,Personas,coste);
        while(vidas)
            pos_contrincante = randperm(height(Matriz_coral),1);
            if(isnan(Matriz_coral(pos_contrincante,1)))
                Matriz_coral(pos_contrincante,:) = hijos(i,:);
                break;
            else
                potencia_padre = function_objetivo(Matriz_coral(pos_contrincante,:),alcance_max,coste_min,alpha,betta,modo,Personas,coste);
                if (potencia_padre>potencia_hijo)
                    vidas = vidas-1;
                else
                    Matriz_coral(pos_contrincante,:) = hijos(i,:);
                    break;
                end
            end
        end
    end
end