classdef matrix_handle
    properties (Access = private)
        multFun		% function handle for y = A*x
        dim			% dimension (square matrix)		
    end
    
    methods
        function obj = matrix_handle(A, n)
            % assume that A is a function handle, n integer
            obj.multFun = A;
            obj.dim = n;
        end
        
        function n = size(obj, idx)
            if nargin == 1
                n = [obj.dim, obj.dim];
            else
                if (idx == 1 || idx == 2)
                    n = obj.dim;
                else
                    n = 1;	% default
                end
            end
        end
        
        function y = mtimes(obj, x)
            y = obj.multFun(x);
        end
    end
end
