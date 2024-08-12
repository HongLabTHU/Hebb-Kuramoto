classdef Traveling_Mode
    
    properties
        Prop
    end
    
    methods
        function obj = Traveling_Mode(val)
            % ¹¹Ôìº¯Êý
            if nargin > 0
                obj.Prop = val;
            end
        end
        
        function grids = set_grids(obj,n_grid)
            % Define grid positions
            n_vec = n_grid ^ 2;
            grids = zeros([n_vec, 2]);
            k = 1;
            for i = 1:n_grid
                for j = 1:n_grid
                    grids(k ,:) = [(i-1) / n_grid, (j-1) / n_grid];
                    k = k + 1;
                end
            end
            grids = grids + 0.5 / n_grid - 0.5;
        end
        
        function v = translational(obj, grids, dirs, wl)
            n_dirs = length(dirs);
            n_vec = size(grids, 1);
            v = zeros(n_dirs, n_vec);
            for i = 1:n_dirs
                proj_vec = [cos(dirs(i)); sin(dirs(i))];
                v(i,:) = exp(1j * grids * proj_vec * pi * 2 / wl)';
            end
        end
        
        function v = rotational(obj, grids, centers, polarity, wn)
            % polarity=1, wn=1
            n_centers = size(centers, 1);
            n_vec = size(grids, 1);
            v = zeros([n_centers, n_vec]);
            for i = 1:n_centers
                dir_vec = grids - repmat(centers(i, :),n_vec,1);
                v(i,:) = exp(1j * angle(dir_vec(:, 1) + 1j * dir_vec(:, 2)) * wn * polarity)';
            end
        end
        
        function v = singular(obj, grids, centers, polarity, wl)
            % polarity=1, wl=1;
            n_centers = size(centers, 1);
            n_vec = size(grids, 1);
            v = zeros([n_centers, n_vec]);
            for i = 1:n_centers
                dir_vec = grids - centers(i, :);
                v(i,:) = exp(1j * sum(dir_vec.^2,2).^0.5 * pi * 2 / wl * polarity)';
            end
        end
        
        function v = saddle(obj, grids, centers, orientations, wl)
            n_centers = size(centers, 1);
            n_vec = length(grids);
            v = zeros(n_centers, n_vec);
            for i = 1:n_centers
                dir_vec = grids - centers(i,:);
                ang_vec = angle(dir_vec(:,1) + 1i * dir_vec(:,2));
                v(i,:) = exp(-1i * sqrt(sum(dir_vec.^2, 2)) .* cos((ang_vec - orientations(i)) * 2 * pi * 2 / wl));
            end
        end
        
        function r = purify(obj, v)
            v = v / sqrt(v * v'); % Normalize v
            c = conj(v) * v'; % Calculate the complex conjugate dot product
            t = (-1 + sqrt(max(0, 1 - real(c * conj(c))))) / c; % Calculate t

            r = v + t * conj(v); % Calculate the purified vector r
            r_norm = r * r'; % Calculate the complex conjugate dot product of r

            if r_norm > 1e-6
                r = r / sqrt(r_norm); % Normalize r if its norm is above the threshold
            else
                r = []; % Return an empty array if the norm of r is below the threshold
            end
        end


    end
end