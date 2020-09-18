% Plots symbol of wavelet filter
% with two real numbers.

close all
clear all

x=linspace(0,1,1000); %%%% create an array for x from 0 to 1 with 1000 points

% Plot real & imaginary parts.
plot(x, real(symbol(x)));
title('Real Part');
figure();
plot(x, imag(symbol(x)));
title('Imaginary Part');
figure();
plot(x, abs(symbol(x)));
title('Absolute Value');


% Creates the symbol as a function.
function dum = symbol(x)
dum = 1/2 + (1/2)*exp(-2*pi*i*x);
end
