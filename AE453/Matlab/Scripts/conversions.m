function conversions
%{
%File containing current International Astronomical Union (IAU) 2000A Earth orientation data
aeroReadIERSData	
%Convert from acceleration units to desired acceleration units
convacc
%Convert from angle units to desired angle units
convang	
%Convert from angular acceleration units to desired angular acceleration units
convangacc		
%Convert from angular velocity units to desired angular velocity units
convangvel	
%Convert from density units to desired density units
convdensity	
%Convert from force units to desired force units
convforce	
%Convert from length units to desired length units
convlength	
%Convert from mass units to desired mass units
convmass	
%Convert from pressure units to desired pressure units
convpres	
%Convert from temperature units to desired temperature units
convtemp	
%Convert from velocity units to desired velocity units
convvel
%}
%{
exa     E	1000000000000000000	   10e18
peta	P	1000000000000000	10e15
tera	T	1000000000000	10e12
giga	G	1000000000	10e9
mega	M	1000000	  10e6
kilo	k	1000	10e3
hecto	h	100	  10e2
deca	da	10	10e1
(none)	1	100
deci	d	0.1	10e?1
centi	c	0.01   10e?2
milli	m	0.001	10e?3
micro	?	0.000001	10e?6
nano	n	0.000000001	   10e?9
pico	p	0.000000000001	10e?12
femto	f	0.000000000000001	10e?15
atto	a	0.000000000000000001	10e?18
%}
end
