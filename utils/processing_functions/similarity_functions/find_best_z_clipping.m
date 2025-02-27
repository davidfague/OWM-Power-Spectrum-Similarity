% function find_best_z_clipping(WI, BI)
%     % find the z-value that will be just outside the non-inf data
%     rho = 0.999; % start with rho = 0.999 and add a 9 to the end until it meets criteria.
%     z = 3.8; % fisher-z of starting rho
%     while ~(max(BI(~isinf(BI))) < z && max(WI(~isinf(WI))) < z && min(WI(~isinf(WI))) > -z && min(BI(~isinf(BI))) > -z)
%         z = z + 
%     end
% end

function best_z = find_best_z_clipping(WI, BI)
    % find_best_z_clipping finds the smallest Fisher z-value 
    % (derived from an increasing rho) that exceeds the finite data in WI and BI.
    %
    % WI, BI: matrices of Fisher z-values (or any similar measures)
    % best_z: the z-value that is just outside the finite values in WI and BI.
    
    % Start with rho = 0.999 (i.e., a very high correlation)
    initial_rho = 0.999;
    rho = initial_rho;
    % Compute its Fisher Z-transform
    z = 0.5 * log((1 + rho) / (1 - rho));
    
    % Continue increasing rho (by appending a '9') until the condition is met:
    % All finite values in both BI and WI must be within (-z, z)
    while ~( max(BI(~isinf(BI))) < z && max(WI(~isinf(WI))) < z && ...
             min(WI(~isinf(WI))) > -z && min(BI(~isinf(BI))) > -z )
         
         % Append a '9' to the end of rho's string representation
         rho_str = num2str(rho);
         rho_str = [rho_str '9'];
         fprintf("trying rho %s\n", rho_str)
         rho = str2double(rho_str);

         if rho == 1
             error("finding best z did not work using str2double method")
         end
         
         % Update z according to the Fisher Z-transform
         z = 0.5 * log((1 + rho) / (1 - rho));
    end
    
    % Return the best z value found
    best_z = z;
    fprintf("clipping rho to %s\n", rho)
end

% function best_z = find_best_z_clipping(WI, BI)
%     % find_best_z_clipping finds the smallest Fisher z-value that is just
%     % outside the range of the finite data in WI and BI.
%     %
%     % WI, BI: matrices containing Fisher z-values (or similar measures)
%     % best_z: the smallest z such that all finite values in WI and BI satisfy:
%     %         -best_z < value < best_z.
% 
%     % Extract finite values from both matrices
%     finite_WI = WI(isfinite(WI));
%     finite_BI = BI(isfinite(BI));
% 
%     % Check if there are any finite values
%     if isempty(finite_WI) && isempty(finite_BI)
%         error('Both WI and BI contain no finite values.');
%     end
% 
%     % Find the maximum absolute finite value from both matrices
%     max_val = max([max(finite_WI(:)), max(finite_BI(:)), ...
%                    -min(finite_WI(:)), -min(finite_BI(:))]);
% 
%     % Add a small tolerance so that best_z is strictly larger than max_val
%     tol = 1e-10;
%     best_z = max_val + tol;
% 
%     % Compute the corresponding correlation value using the inverse Fisher transform
%     best_rho = tanh(best_z);
%     fprintf("clipping rho to %g\n", best_rho);
% end
