function [S] = unpackFrame(frameData)
    
    N = numel(frameData);
    
    %vectorized
    S = zeros(N*2,2);
    S([[1:2:N*2, 2:2:N*2] N*2+[1:2:N*2, 2:2:N*2]]) =  typecast(uint16([bitand(frameData,65535) ; bitshift(bitand(frameData,4294901760), -16) ; ...
        bitshift(bitand(frameData,281470681743360), -32) ; bitshift(bitand(frameData,18446462598732840960), -48)]), 'int16');
    
%     signalData = zeros(N*2,2);
%     for i = 1:N
%         dat = typecast(uint16([bitand(frameData(i),65535)...
%             bitshift(bitand(frameData(i),4294901760), -16)...
%             bitshift(bitand(frameData(i),281470681743360), -32)...
%             bitshift(bitand(frameData(i),18446462598732840960), -48)]),'int16');
%         signalData([i*2-1 i*2],:) = [dat(1) dat(3); dat(2) dat(4);];
%     end

end

