       m = 4.48e-26;     % Mass of Electron [kg]
        h = 6.62607e-34;%Planck's Constant [m^2 kg/s]
        epsilon = 15.76; %Ionization Energy wrt. Atomic Ground State [eV]
        T = 1;
        P=1;
        P = P*760*133.322; %torr to Pa [kg/ m s^2]
        %{
        syms term2 real
        term1 = ((2*((2*pi*m)^(3/2))*((T)^(5/2)))/(P*h^3));
        
        term3 = exp(-epsilon/(T));
        eqn1 = term1*term2*term3;
        sol = vpa(solve(eqn1==(1/3), term2))
 %}
        term2=3.9939318527271464835152172363226e-54
        term1 = ((2*((2*pi*m)^(3/2))*((T)^(5/2)))/(P*h^3));
        
        term3 = exp(-epsilon/(T));

        syms alpha real

        eqn1 = ((alpha^2)/(1-alpha^2)) == term1*term2*term3;
        sol = vpa(solve(eqn1, alpha));
        alpha = sol(2)