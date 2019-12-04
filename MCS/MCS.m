classdef MCS < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x;
        y;
        err_arr;
        a_arr;
        b_arr;
        posi;
        posj;
        v;
        opt_coefs;
    end
    
    methods
        % Constructor
        function obj = MCS (filename)
            obj = obj.getData(filename);
        end
        
        function obj = getData(obj, filename)
            M = csvread(filename);
            
            obj.x = M(:,1);
            obj.y = M(:,2);
            obj.v = std(obj.y).^2;
            
            n = length(obj.x);
            obj.err_arr = zeros(n,n);
            obj.a_arr = zeros(n,n);
            obj.b_arr = zeros(n,n);
            
            for j=1:n             
                for i=1:j
                    fprintf('j = %d, i = %d\n', j,i);
                    [a, b, e2] = obj.lscoef(obj.x(i:j),obj.y(i:j));
                    obj.err_arr(i,j) = e2;
                    obj.a_arr(i,j) = a;
                    obj.b_arr(i,j) = b;
                end
            end
        end
        
        function [a,b,e2] = lscoef(~,x,y)
            n = length(x);
            
            if (n == 1)
               a = 0.0;
               b = y(1);
               e2 = 0.0;
            else
               sx = sum(x);
               sy = sum(y);
               sx2 = sum(x.^2);
               sxy = sum(x.*y);

               % calculo del error minimo
               a = (n*sxy-sx*sy)/(n*sx2-sx*sx);
               b = (sy-a*sx)/n;
               e2 = sum((y-a*x-b).^2);
            end
        end
        
        function find_opt(obj, penalty_factor)
            penalty = obj.v .* penalty_factor;
            n = size(obj.err_arr);
            n = n(1);

            opt_arr = zeros(n,1);
            for j=1:n
                tmp_opt = zeros(j,1);
                tmp_opt(1) = obj.err_arr(1,j) + penalty;
                for i=2:j
                    tmp_opt(i) = opt_arr(i-1) + obj.err_arr(i,j) + penalty;
                end
                opt_arr(j) = min(tmp_opt);
            end
            
            opt_coefs = [];
            posi = zeros(n,1);
            posj = zeros(n,1);
            j = n;
            while j >= 1
                tmp_opt = zeros(j,1);
                tmp_opt(1) = obj.err_arr(1,j) + penalty;
                for i=2:j
                    tmp_opt(i) = opt_arr(i-1) + obj.err_arr(i,j) + penalty;
                end
                [~,i_opt] = min(tmp_opt);
                a_opt = obj.a_arr(i_opt,j);
                b_opt = obj.b_arr(i_opt,j);
                posi(i_opt) = j;
                posj(j) = i_opt;
                obj.posi = posi;
                obj.posj = posj;
                
                if (i_opt <= 1)
                    xmin = -Inf;
                else
                    xmin = (obj.x(i_opt-1) + obj.x(i_opt))/2;
                end
                
                if (j >= n)
                    xmax = Inf;
                else
                    xmax = (obj.x(j) + obj.x(j+1))/2;
                end
                
                opt_coefs = [[xmin,xmax,a_opt,b_opt] ; opt_coefs];

                j = i_opt-1;
                
            end
            
            obj.opt_coefs = opt_coefs;
                
        end
        
        function yfit = get_fit(obj, x)
            f = 1;
            n = length(x);

            obj.posi = sort(obj.posi);
            obj.posj = sort(obj.posj);
            yfit = zeros(n,1);
            
            for j = 1 : length(obj.opt_coefs)
                opt_coef = obj.opt_coefs(j,:);
                ind = [];
                for i = 1 : length(x)
                    if (x(i) >= opt_coef(1) && x(i) <= opt_coef(2))
                        ind = [ind i];
                    end
                end
                
                for i = 1 : length(ind)
                    yfit(ind(i)) = x(ind(i)) * opt_coef(3) + opt_coef(4);
                end
                
                fprintf('y = %f * x + %f\n', opt_coef(3), opt_coef(4));
                %disp('Para los puntos '), disp(obj.posj(f)+1), disp(' a '), disp(obj.posi(f)+1);
                
                f = f+1;
            end
        end
        
        function plot_fit(obj, xplt)
            yfit = obj.get_fit(xplt);
            
            figure(1);
            plot(obj.x, obj.y, 'r.', xplt, yfit,'g');
            xlabel('x');
            ylabel('y');
            title('Minimos Cuadrados Segmentados (SLS)');
            %legend(('Puntos ingresados', 'L?nea del mejor ajuste'),loc='upper right');
            %plt.show()
        end
        
    end
    
end

