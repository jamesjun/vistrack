function prefix = getSpike2Prefix(S)
    prefix = fields(S); 
    prefix = prefix{1}; 
    k = strfind(prefix, '_Ch');
    k=k(end);
    prefix = prefix(1:k-1);    
end