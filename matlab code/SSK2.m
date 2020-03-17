tic
clear; clc
M = 16;  %the number of mudualtion constellation 
Mt = 4;% the nubmer of transmit antenna
Mr = 2;% the number of receive antenna
SNR = 0: 2 :30;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt    = log2(Mt);            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = log2(Mt); 
Nbits = 50000;
%generate bit stream
bit_T = randi([0 1],Nbits,log2(Mt));
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;
%ssk mapping
x=zeros(Nbits,Mt);

for i=1:Nbits
   x(i,active_ant(i))=1;
end
%tranmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 

for ii=1:L_SNR 
   for j = 1 : size(x,1) 
    channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))*sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*x(j,:)'+noise;      
   end
end
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,2));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
       for jj=1:Mt
         noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1));  
       x_detection(jj)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj))'* channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+jj));
       end
       r=max(x_detection); 
       x_index(ii,j)=find(x_detection==r);
   end
end
%dem llr
 bit_Tdem=zeros(size(x,1),Nt*L_SNR);
for ii=1:L_SNR 
    for j = 1 : size(x,1) 
       bit_Tdem(j,ii*Nt-(Nt-1):ii*Nt)=  de2bi(x_index(ii,j)-1,log2(Mt),'left-msb');
    end
end

 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*Nt-(Nt-1):ii*Nt)-bit_T)~=0));
    ber(ii)=a(ii)/(Nt*Nbits);
end
semilogy(SNR,ber,'-co','LineWidth',1); hold on
clear; clc
M = 16;  %the number of mudualtion constellation 
Mt = 8;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
SNR = 0: 2 :30;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt    = log2(Mt);            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = log2(Mt); 
Nbits = 50000;
%generate bit stream
bit_T = randi([0 1],Nbits,log2(Mt));
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;
%ssk mapping
x=zeros(Nbits,Mt);

for i=1:Nbits
   x(i,active_ant(i))=1;
end
%tranmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 

for ii=1:L_SNR 
   for j = 1 : size(x,1) 
    channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))*sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*x(j,:)'+noise;      
   end
end
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,2));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
       for jj=1:Mt
         noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1));  
       x_detection(jj)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj))'* channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+jj));
       end
       r=max(x_detection); 
       x_index(ii,j)=find(x_detection==r);
   end
end
%dem llr
 bit_Tdem=zeros(size(x,1),Nt*L_SNR);
for ii=1:L_SNR 
    for j = 1 : size(x,1) 
       bit_Tdem(j,ii*Nt-(Nt-1):ii*Nt)=  de2bi(x_index(ii,j)-1,log2(Mt),'left-msb');
    end
end

 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*Nt-(Nt-1):ii*Nt)-bit_T)~=0));
    ber(ii)=a(ii)/(Nt*Nbits);
end
semilogy(SNR,ber,'-bv','LineWidth',1); hold on
clear; clc
M = 16;  %the number of mudualtion constellation 
Mt = 16;% the nubmer of transmit antenna
Mr = 8;% the number of receive antenna
SNR = 0: 2 :30;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt    = log2(Mt);            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = log2(Mt); 
Nbits = 50000;
%generate bit stream
bit_T = randi([0 1],Nbits,log2(Mt));
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;
%ssk mapping
x=zeros(Nbits,Mt);

for i=1:Nbits
   x(i,active_ant(i))=1;
end
%tranmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 

for ii=1:L_SNR 
   for j = 1 : size(x,1) 
    channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))*sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*x(j,:)'+noise;      
   end
end
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,2));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
       for jj=1:Mt
         noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1));  
       x_detection(jj)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj))'* channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+jj));
       end
       r=max(x_detection); 
       x_index(ii,j)=find(x_detection==r);
   end
end
%dem llr
 bit_Tdem=zeros(size(x,1),Nt*L_SNR);
for ii=1:L_SNR 
    for j = 1 : size(x,1) 
       bit_Tdem(j,ii*Nt-(Nt-1):ii*Nt)=  de2bi(x_index(ii,j)-1,log2(Mt),'left-msb');
    end
end

 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*Nt-(Nt-1):ii*Nt)-bit_T)~=0));
    ber(ii)=a(ii)/(Nt*Nbits);
end
semilogy(SNR,ber,'-kd','LineWidth',1); hold on
clear ;   
clc;
M = 4;
Mt = 1; %for number of transmit antenna
Mr = 2;%for number of receive antenna
Nbits =50000;
SNR =0: 2 :30;                  % signal-to-noise ratio in dB    
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
bit_SMsym = log2(M)*Mt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M); 
bit_T=randi([0 1],Nbits*Mt,Numberofbit);
bit_T=bit_T.';
x=zeros(Nbits*Mt,1);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType', 'Full precision' );
No=10.^(-SNR/10);               % noise variance
for i=1:Nbits*Mt
Hmod(i) = step(hMod,bit_T(1:end,i));
x(i)=Hmod(i);
end
x=reshape(x,Mt,Nbits);
x=x.';
bit_T=bit_T.';


%transmission
y=zeros(L_SNR*Mr,size(x,1));
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end


%detection
x_detected=zeros(L_SNR,size(x,1));
for ii=1:L_SNR 
   for j = 1 : Nbits
       h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt);
       x_detected(ii*Mt-(Mt-1):ii*Mt,j)=(h'*h+No(ii)*eye(Mt,Mt))\h'*y(ii*Mr-(Mr-1):ii*Mr,j);
   end
   x_detect(ii,:)=reshape(x_detected(ii*Mt-(Mt-1):ii*Mt,:),1,Nbits*Mt);
end

%demodulation
bit_Tdem=zeros(size(x,1),L_SNR*Mt);
for ii=1:L_SNR 
   for j = 1 : Nbits*Mt   
       bit_Tdem(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)=step(hdeMod,(x_detect(ii,j)));
   end
end
 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)-bit_T)~=0));
    ber(ii)=a(ii)/(bit_SMsym*Nbits);
end 
semilogy(SNR,ber,'--co','LineWidth',1);hold on
M = 8;
Mt = 1; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nbits =50000;
SNR =0: 2 :30;                  % signal-to-noise ratio in dB    
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
bit_SMsym = log2(M)*Mt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M); 
bit_T=randi([0 1],Nbits*Mt,Numberofbit);
bit_T=bit_T.';
x=zeros(Nbits*Mt,1);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType', 'Full precision' );
No=10.^(-SNR/10);               % noise variance
for i=1:Nbits*Mt
Hmod(i) = step(hMod,bit_T(1:end,i));
x(i)=Hmod(i);
end
x=reshape(x,Mt,Nbits);
x=x.';
bit_T=bit_T.';


%transmission
y=zeros(L_SNR*Mr,size(x,1));
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end


%detection
x_detected=zeros(L_SNR,size(x,1));
for ii=1:L_SNR 
   for j = 1 : Nbits
       h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt);
       x_detected(ii*Mt-(Mt-1):ii*Mt,j)=(h'*h+No(ii)*eye(Mt,Mt))\h'*y(ii*Mr-(Mr-1):ii*Mr,j);
   end
   x_detect(ii,:)=reshape(x_detected(ii*Mt-(Mt-1):ii*Mt,:),1,Nbits*Mt);
end

%demodulation
bit_Tdem=zeros(size(x,1),L_SNR*Mt);
for ii=1:L_SNR 
   for j = 1 : Nbits*Mt   
       bit_Tdem(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)=step(hdeMod,(x_detect(ii,j)));
   end
end
 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)-bit_T)~=0));
    ber(ii)=a(ii)/(bit_SMsym*Nbits);
end 
semilogy(SNR,ber,'--bv','LineWidth',1);hold on
M = 16;
Mt = 1; %for number of transmit antenna
Mr = 8;%for number of receive antenna
Nbits =50000;
SNR =0: 2 :30;                  % signal-to-noise ratio in dB    
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
bit_SMsym = log2(M)*Mt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M); 
bit_T=randi([0 1],Nbits*Mt,Numberofbit);
bit_T=bit_T.';
x=zeros(Nbits*Mt,1);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType', 'Full precision' );
No=10.^(-SNR/10);               % noise variance
for i=1:Nbits*Mt
Hmod(i) = step(hMod,bit_T(1:end,i));
x(i)=Hmod(i);
end
x=reshape(x,Mt,Nbits);
x=x.';
bit_T=bit_T.';


%transmission
y=zeros(L_SNR*Mr,size(x,1));
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end


%detection
x_detected=zeros(L_SNR,size(x,1));
for ii=1:L_SNR 
   for j = 1 : Nbits
       h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt);
       x_detected(ii*Mt-(Mt-1):ii*Mt,j)=(h'*h+No(ii)*eye(Mt,Mt))\h'*y(ii*Mr-(Mr-1):ii*Mr,j);
   end
   x_detect(ii,:)=reshape(x_detected(ii*Mt-(Mt-1):ii*Mt,:),1,Nbits*Mt);
end

%demodulation
bit_Tdem=zeros(size(x,1),L_SNR*Mt);
for ii=1:L_SNR 
   for j = 1 : Nbits*Mt   
       bit_Tdem(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)=step(hdeMod,(x_detect(ii,j)));
   end
end
 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)-bit_T)~=0));
    ber(ii)=a(ii)/(bit_SMsym*Nbits);
end 
semilogy(SNR,ber,'--kd','LineWidth',1);hold on
grid on;
xlabel('$$SNR$$','Interpreter','latex')
ylabel('BER','Interpreter','latex')
title('Bit Error Rate of SM using ML detector with receive antennas')
legend('8X4 SSK','4X2 SSK','16X8 SSK','MMSE-SIMO  4QAM 1X2 ' , 'MMSE-SIMO  8QAM 1X4 ','MMSE-SIMO  16QAM 1X8 ', 'Location','NorthEast')
toc