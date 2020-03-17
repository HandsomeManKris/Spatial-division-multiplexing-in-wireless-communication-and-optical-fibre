clear;clc
M = 16;
Mt = 2; %for number of transmit antenna
Mr = 8;%for number of receive antenna
Nbits =100;
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
       x_detected(ii*Mt-(Mt-1):ii*Mt,j)=inv(h'*h+No(ii)*eye(Mt,Mt))*h'*y(ii*Mr-(Mr-1):ii*Mr,j);
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
semilogy(SNR,ber,'-g+','LineWidth',1.5);hold on

M = 16;
Mt = 2; %for number of transmit antenna
Mr = 8;%for number of receive antenna
Nbits =100;
SNR =0: 2 :30;                  % signal-to-noise ratio in dB    
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
bit_SMsym = log2(M)*Mt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M); 
bit_T=randi([0 1],Nbits*Mt,Numberofbit);
bit_T=bit_T.';
x=zeros(Nbits*Mt,1);
hMod = comm.PSKModulator('ModulationOrder',M,'BitInput',true,'SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.PSKDemodulator('ModulationOrder',M,'BitOutput',true,'SymbolMapping','Gray','OutputDataType', 'Full precision' );
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
       x_detected(ii*Mt-(Mt-1):ii*Mt,j)=inv(h'*h+No(ii)*eye(Mt,Mt))*h'*y(ii*Mr-(Mr-1):ii*Mr,j);
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
semilogy(SNR,ber,'--ro','LineWidth',1.5);hold on