function window = kaiser3D(sizes,beta);
% creates a 3D Kaiser Window
window = ones(sizes);


    ventana = kaiser(sizes(1),beta);

    ventana = repmat(ventana,[1 sizes(2) sizes(3)]);
    window = window.*ventana;
    
    ventana = kaiser(sizes(2),beta);

    ventana = repmat(ventana',[sizes(1) 1 sizes(3)]);
    

    window = window.*ventana;

    ventana = kaiser(sizes(3),beta);
    ventana = reshape(ventana,[1 1 sizes(3)]);
    
    ventana = repmat(ventana,[sizes(1) sizes(2) 1]);
    
    window = window.*ventana;
    