function HFun = createFunHandle(H, funType)
%creatFunHandle: Create a function handle from a matrix H, which is for H*x if funType is 'mtimes' or H\x if funType is 'mldivide'.
% If H is already a function handle, then do nothing but return H itself.
%HFun = solver.createFunHandle(H, 'mtimes')
%HFun = solver.createFunHandle(H, 'mldivide')

if isa(H, 'function_handle')
    HFun = H;
elseif isnumeric(H)
    switch lower(funType)
        case 'mtimes'
            HFun = @(x) H*x;        
        case 'mldivide'
            HFun = @(x) H\x;
        otherwise
            % not implemented other function types, how about add @(x) H'*x in future etc.
            error('Input function type %s is not implemented!', funType); 
    end
else        
    error('Input need to be either a numeric array or a function handle');
end
end