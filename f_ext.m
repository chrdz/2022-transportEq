function res = f_ext(t)
% function M_1
% to define the external force

    if t < 2.5
        res = 200/2.5*t;
    elseif t < 5
        res = 400 - 200/2.5*t;
    else
        res = 0;
    end

end


% function res = f_ext(t)
% % function M_1
% % to define the external force
% 
%     tmax = 1;
%     x0 = 0.5*tmax;   % center
%     a = 0.2*tmax;    % width
%     c0 = 0.1;       % magnitude
%     f = @(tt) c0 * exp( -1/a^2*((tt-x0).^2) ) ; % bump function
%     res = f(t);
%     
% 
% end


