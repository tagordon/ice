
function [Evib] = calculate_Evib(sp,nu_start,nu_end,T)
%input: a DOS in b-spline format, start and end points for the frequency (cm^-1),
%and the temperature (K)
%output: the internal energy due to ionic motion in J/kg

%read in constants and system-specific information
my_constants;

%convert Planck constant from J s to J cm (since frequency is in cm^-1)
h = h*c;

%function to integrate representing the energy of each state and its
%occupational probability multiplied by the DOS
fun = @(nu) h*nu.*((exp((h.*nu)/k/T)-1)).^-1.*fnval(sp,nu);

%integral of the function over the range of frequencies
Evib = integral(fun,nu_start,nu_end);

%convert from J/particle to J/kg
Evib = Evib*N/mw*1000;

 