% Calculate analytical solution

[x,y,z]= meshgrid(0:0.05:1,0:0.01:1,0:0.2:1);


Panal = cosh(2*pi.*x)./cosh(2*pi).*cos(2*pi.*y);

