function IMap = unflattenIpos(Ipos,okpix,dims)
%
%
%

%%
    temp = double(okpix);
    temp(okpix) = Ipos;
    iszero = temp==0;
    isNaN = isnan(temp);
    
    %Flip zeros and nans.
    temp(iszero) = nan;
    temp(isNaN) = 0;
    
    IMap = reshape(temp,dims);
    
end