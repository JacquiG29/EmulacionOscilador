% =========================================================================
% IE3011 - TAREA 2: linealizaci�n num�rica de sistemas no lineales
% -------------------------------------------------------------------------
% Complete la funci�n linloc, encargada de efectuar una linealizaci�n local
% de un sistema no lineal con campos vectoriales f(x,u) y h(x,u) alrededor
% de un punto de equilibrio (xe,ue), bas�ndose en el algoritmo descrito en
% el problema 5 del documento de la tarea. Los par�metros f_handle y
% h_handle son punteros que llaman a las funciones que apuntan. Por
% ejemplo, para este ejemplo particular usted debe usar la sintaxis
% >> [A,B,C,D] = linloc(@dinamica,@salida,xe,ue)
% ya sea en la l�nea de comando o dentro de un script para aplicar la
% funci�n de linealizaci�n a las funciones que complet� previamente.
% =========================================================================
function J = num_jacobian(f_handle, x,G)
    %% Inicializaci�n
    % Se obtienen los tama�os de cada uno de los vectores del sistema
    n = length(x);
    delta = 1e-12; 
    J = zeros(size(G'));
    % Desviaci�n empleada para la aproximaxi�n num�rica de los jacobianos
    f = f_handle(x);
    for k = 1:n
        dx = zeros(n,1);
        dx(k) = delta;
        df = f_handle(x+dx);
        J(:,k) = (df-f)/delta;
    end

        
end