classdef solver_param
    %SOLVER_PARAM undefined
    %   undefined

    properties
        %Armijo backtracking parameters
        backtracking_max {mustBeNumeric}
        beta {mustBeNumeric} 
        alpha {mustBeNumeric}
        sigma {mustBeNumeric}
        allowed_failed_backtracks {mustBeNumeric}

        iter_max {mustBeNumeric}
    end
end    