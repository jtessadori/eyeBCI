clear
close all
clc

inPort=7561;
outPort=7560;

net=load('lastNet.mat','net');
% bestCoeffs=load('lastNet.mat','bestCoeffs');

udps=dsp.UDPSender('RemoteIPPort',outPort);
udpr=dsp.UDPReceiver('LocalIPPort',inPort,'ReceiveBufferSize',26000,'MaximumMessageLength',26000);

tic;
sigmoid=@(coeffs,x)1./(1+exp(-coeffs(1)*(x-(coeffs(2)))));
while true
    inString=[];
    while isempty(inString)
        inString=udpr.step;
        pause(.01);
    end
    if strcmp(char(inString)','close')
        break
    else
%         prediction=predict(net.net,eval(char(inString))');
%         outValue=0.5+1.5*sigmoid(bestCoeffs.bestCoeffs,prediction);
        outValue=double(classify(net.net,eval(char(inString'))'));
    end
    outString=uint8(mat2str(outValue));
    fprintf('%s\n',outString)
    udps.step(outString);
end

udpr.release