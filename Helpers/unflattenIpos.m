function IMap = unflattenIpos(Ipos,okpix,dims)
%
%
%

%%
    IMap = double(okpix);
    IMap(okpix) = Ipos;
    iszero = IMap==0;
    isNaN = isnan(IMap);
    
    %Flip zeros and nans.
    IMap(iszero) = nan;
    IMap(isNaN) = 0;
    
    IMap = reshape(IMap,dims);
    
end