function pCorrect = DayDecodeError(decodedDay,realDay)
%
%
%

%%
    correct = (decodedDay-realDay) == 0;
    pCorrect = (sum(correct)/length(realDay))*100; 
    
end