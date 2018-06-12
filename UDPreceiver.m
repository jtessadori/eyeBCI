clear
close all
clc

inPort=7561;
outPort=7560;

net=load('lastNet.mat','net');

udps=dsp.UDPSender('RemoteIPPort',outPort);
udpr=dsp.UDPReceiver('LocalIPPort',inPort,'ReceiveBufferSize',26000,'MaximumMessageLength',26000);

tic;
while true
    inString=[];
    while isempty(inString)
        inString=udpr.step;
        pause(.01);
    end
    if strcmp(char(inString)','close')
        break
    else
        outValue=double(classify(net.net,eval(char(inString))'));
    end
    switch outValue
        case 1
            outString=uint8(mat2str(2));
        case 2
            outString=uint8(mat2str(0.5));
    end
    fprintf('%s\n',outString)
    udps.step(outString);
end

udpr.release