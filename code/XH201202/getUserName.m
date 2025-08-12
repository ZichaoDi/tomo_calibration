function username = getUserName()
% get user name, works under Mac, Linux and Windows 10.
    if isunix() 
        username = getenv('USER'); 
    else 
        username = getenv('username'); 
    end
end