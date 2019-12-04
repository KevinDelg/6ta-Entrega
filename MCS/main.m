clc;
clear all;
clearvars;
warning off all;


% Se escoge la archivo a analisar
[filename,pathname] = uigetfile({'*.txt';'*.csv'});
s = MCS(strcat(pathname,filename));
s.find_opt(0.35);

n = length(s.x);
minx = min(s.x);
maxx = max(s.x);
dx = 0.33*(maxx-minx)/(n-1);

xplt = [];
for i = minx : dx : maxx + dx
    xplt = [xplt i];
end

disp('El tiempo del algoritmo es O(n^3)');
disp('Se muestran las ecuaciones de la recta del mejor ajuste para determinados puntos');
s.plot_fit(xplt); %tiempo dado en O(n^2)