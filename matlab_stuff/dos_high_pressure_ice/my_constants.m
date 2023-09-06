
%constants
k = 1.38066e-23;    %Boltzmann constant in J/K
h = 6.62618e-34;    %Planck constant in J s
c = 299792458*100;  %speed of light in cm/s
N = 6.0221409e+23;  %Avogadro's number

%system-specific information
mw = 18.01528;      %molar mass in g/mol
%number of atom types in our system
num_atom_types = 2;
%number of atoms of each type and color scheme to use in plotting
for i = 1:num_atom_types
    switch i
        case 1
            color(i) = 'r';
            num_atom(i) = 1;
        case 2
            color(i) = 'b';
            num_atom(i) = 2;
        case 3
            color(i) = 'g';
    end
end