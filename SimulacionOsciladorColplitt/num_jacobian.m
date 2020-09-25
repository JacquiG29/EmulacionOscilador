% =========================================================================
% IE3011 - TAREA 2: linealización numérica de sistemas no lineales
% -------------------------------------------------------------------------
% Complete la función linloc, encargada de efectuar una linealización local
% de un sistema no lineal con campos vectoriales f(x,u) y h(x,u) alrededor
% de un punto de equilibrio (xe,ue), basándose en el algoritmo descrito en
% el problema 5 del documento de la tarea. Los parámetros f_handle y
% h_handle son punteros que llaman a las funciones que apuntan. Por
% ejemplo, para este ejemplo particular usted debe usar la sintaxis
% >> [A,B,C,D] = linloc(@dinamica,@salida,xe,ue)
% ya sea en la línea de comando o dentro de un script para aplicar la
% función de linealización a las funciones que completó previamente.
% =========================================================================
function J = num_jacobian(f_handle, x,G)
    %% Inicialización
    % Se obtienen los tamaños de cada uno de los vectores del sistema
    n = length(x);
    delta = 1e-12; 
    J = zeros(size(G'));
    % Desviación empleada para la aproximaxión numérica de los jacobianos
    f = f_handle(x);
    for k = 1:n
        dx = zeros(n,1);
        dx(k) = delta;
        df = f_handle(x+dx);
        J(:,k) = (df-f)/delta;
    end

        
end