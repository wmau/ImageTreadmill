function R = flipR(R,src,snk)
%R = flipR(R,src,snk)
%
%

%%
    R{src,snk}.trials =         fliplr(R{snk,src}.trials);
    R{src,snk}.curve =          fliplr(R{snk,src}.curve); 
    R{src,snk}.CI =             fliplr(R{snk,src}.CI);
    R{src,snk}.tshff.mu =       fliplr(R{snk,src}.tshff.mu);
    R{src,snk}.tshff.upper =    fliplr(R{snk,src}.tshff.upper);
    R{src,snk}.tshff.lower =    fliplr(R{snk,src}.tshff.upper);
    R{src,snk}.sigt =           fliplr(R{snk,src}.sigt);
    R{src,snk}.trlshff.mu =     fliplr(R{snk,src}.trlshff.mu);
    R{src,snk}.trlshff.upper =  fliplr(R{snk,src}.trlshff.upper);
    R{src,snk}.trlshff.lower =  fliplr(R{snk,src}.trlshff.upper);
    R{src,snk}.sigtrls =        fliplr(R{snk,src}.sigtrls);
    R{src,snk}.sig =            fliplr(R{snk,src}.sig);
    
end