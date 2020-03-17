tic
clear;   
clc;
M = 16;
Mt = 16;                         %for number of transmit antenna
Mr = 4;                         %for number of receive antenna
bit_SMsym = log2(M*Mt);         % number of bit per spatial modulation sysmbol
Nbits =20000;
SNR = 0 : 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = log2(Mt);                  %bits for transmit antenna
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
bit_T=bit_T.';
x=zeros(Nbits,Mt);
hmodem =   modem.qammod('M',M,  'SymbolOrder', 'Gray','InputType', 'bit');
AP = (mean(hmodem.Constellation .* conj(hmodem.Constellation)));
No=10.^(-SNR/10);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType', 'Full precision' );
for i=1:Nbits
Hmod(i) = step(hMod,bit_T(Nt+1:end,i));
end
% noise variance
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits
    x(i,active_ant(i))=Hmod(i);
end
bit_T=bit_T.';
%transmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
   channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise;      
   end
end
%detection
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt
            
            x_detection(jj)=(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)'*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
        end
         noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1))* sqrt(No(ii));
        r=max(abs(x_detection));
        x_index(ii,j)=find(abs(x_detection)==r); 
        p=channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j))'*channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j));
        x_detected(ii,j)=channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j))'*y(ii*Mr-Mr+1:ii*Mr,j)/p;
        %x_detected(ii,j)=mean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=harmmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=trimmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)),0.25);
   end
end
%demudulation
bit_Tdem1=zeros(size(x,1),L_SNR*Nt);
bit_Tdem2=zeros(size(x,1),L_SNR*Numberofbit);
bit_Tdem=zeros(size(x,1),L_SNR*bit_SMsym);
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
       noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1))* sqrt(No(ii));
     bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt)= de2bi(x_index(ii,j)-1,log2(Mt),'left-msb');
    
     bit_Tdem2(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)=step(hdeMod,(x_detected(ii,j)));
   end
end

for ii=1:L_SNR 
    for j = 1 : size(x,1) 
       bit_Tdem(j,ii*bit_SMsym-bit_SMsym+1:ii*bit_SMsym-bit_SMsym+Nt)= bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt);
       bit_Tdem(j,ii*bit_SMsym-bit_SMsym+Nt+1:ii*bit_SMsym)=bit_Tdem2(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit);
    end
end
       
 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*bit_SMsym-(bit_SMsym-1):ii*bit_SMsym)-bit_T)~=0));
    ber(ii)=a(ii)/(bit_SMsym*Nbits);
    b(ii)=length(find((bit_Tdem2(:,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)-bit_T2)~=0));
    c(ii)=length(find((bit_Tdem1(:,ii*Nt-(Nt-1):ii*Nt)-bit_T1)~=0));
end 

 semilogy(SNR,ber,'-','LineWidth',1.5);hold on

clear;   
clc;
M = 32;
Mt = 8;                         %for number of transmit antenna
Mr = 4;                         %for number of receive antenna
bit_SMsym = log2(M*Mt);         % number of bit per spatial modulation sysmbol
Nbits =20000;
SNR = 0 : 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = log2(Mt);                  %bits for transmit antenna
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
bit_T=bit_T.';
x=zeros(Nbits,Mt);
hmodem =   modem.qammod('M',M,  'SymbolOrder', 'Gray','InputType', 'bit');
AP = (mean(hmodem.Constellation .* conj(hmodem.Constellation)));
No=10.^(-SNR/10);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType', 'Full precision' );
for i=1:Nbits
Hmod(i) = step(hMod,bit_T(Nt+1:end,i));
end
% noise variance
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits
    x(i,active_ant(i))=Hmod(i);
end
bit_T=bit_T.';
%transmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
   channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise;      
   end
end
%detection
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt
            
            x_detection(jj)=(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)'*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
        end
         noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1))* sqrt(No(ii));
        r=max(abs(x_detection));
        x_index(ii,j)=find(abs(x_detection)==r); 
        p=channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j))'*channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j));
        x_detected(ii,j)=channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j))'*y(ii*Mr-Mr+1:ii*Mr,j)/p;
        %x_detected(ii,j)=mean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=harmmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=trimmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)),0.25);
   end
end
%demudulation
bit_Tdem1=zeros(size(x,1),L_SNR*Nt);
bit_Tdem2=zeros(size(x,1),L_SNR*Numberofbit);
bit_Tdem=zeros(size(x,1),L_SNR*bit_SMsym);
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
       noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1))* sqrt(No(ii));
     bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt)= de2bi(x_index(ii,j)-1,log2(Mt),'left-msb');
    
     bit_Tdem2(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)=step(hdeMod,(x_detected(ii,j)));
   end
end

for ii=1:L_SNR 
    for j = 1 : size(x,1) 
       bit_Tdem(j,ii*bit_SMsym-bit_SMsym+1:ii*bit_SMsym-bit_SMsym+Nt)= bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt);
       bit_Tdem(j,ii*bit_SMsym-bit_SMsym+Nt+1:ii*bit_SMsym)=bit_Tdem2(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit);
    end
end
       
 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*bit_SMsym-(bit_SMsym-1):ii*bit_SMsym)-bit_T)~=0));
    ber(ii)=a(ii)/(bit_SMsym*Nbits);
    b(ii)=length(find((bit_Tdem2(:,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)-bit_T2)~=0));
    c(ii)=length(find((bit_Tdem1(:,ii*Nt-(Nt-1):ii*Nt)-bit_T1)~=0));
end 

 semilogy(SNR,ber,'-','LineWidth',1.5);hold on

clear;   
clc;
M = 64;
Mt = 4;                         %for number of transmit antenna
Mr = 4;                         %for number of receive antenna
bit_SMsym = log2(M*Mt);         % number of bit per spatial modulation sysmbol
Nbits =20000;
SNR = 0 : 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = log2(Mt);                  %bits for transmit antenna
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
bit_T=bit_T.';
x=zeros(Nbits,Mt);
hmodem =   modem.qammod('M',M,  'SymbolOrder', 'Gray','InputType', 'bit');
AP = (mean(hmodem.Constellation .* conj(hmodem.Constellation)));
No=10.^(-SNR/10);
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType', 'Full precision' );
for i=1:Nbits
Hmod(i) = step(hMod,bit_T(Nt+1:end,i));
end
% noise variance
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits
    x(i,active_ant(i))=Hmod(i);
end
bit_T=bit_T.';
%transmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
   channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise;      
   end
end
%detection
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt
            
            x_detection(jj)=(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)'*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
        end
         noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1))* sqrt(No(ii));
        r=max(abs(x_detection));
        x_index(ii,j)=find(abs(x_detection)==r); 
        p=channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j))'*channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j));
        x_detected(ii,j)=channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j))'*y(ii*Mr-Mr+1:ii*Mr,j)/p;
        %x_detected(ii,j)=mean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=harmmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=trimmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)),0.25);
   end
end
%demudulation
bit_Tdem1=zeros(size(x,1),L_SNR*Nt);
bit_Tdem2=zeros(size(x,1),L_SNR*Numberofbit);
bit_Tdem=zeros(size(x,1),L_SNR*bit_SMsym);
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
       noise = sqrt(.5)*(randn(1,1) + 1i*randn(1,1))* sqrt(No(ii));
     bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt)= de2bi(x_index(ii,j)-1,log2(Mt),'left-msb');
    
     bit_Tdem2(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)=step(hdeMod,(x_detected(ii,j)));
   end
end

for ii=1:L_SNR 
    for j = 1 : size(x,1) 
       bit_Tdem(j,ii*bit_SMsym-bit_SMsym+1:ii*bit_SMsym-bit_SMsym+Nt)= bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt);
       bit_Tdem(j,ii*bit_SMsym-bit_SMsym+Nt+1:ii*bit_SMsym)=bit_Tdem2(j,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit);
    end
end
       
 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*bit_SMsym-(bit_SMsym-1):ii*bit_SMsym)-bit_T)~=0));
    ber(ii)=a(ii)/(bit_SMsym*Nbits);
    b(ii)=length(find((bit_Tdem2(:,ii*Numberofbit-(Numberofbit-1):ii*Numberofbit)-bit_T2)~=0));
    c(ii)=length(find((bit_Tdem1(:,ii*Nt-(Nt-1):ii*Nt)-bit_T1)~=0));
end 

 semilogy(SNR,ber,'-','LineWidth',1.5);hold on
 clear; clc
M = 16;  %the number of mudualtion constellation 
Mt = 256;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
SNR = 0: 2 :24;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt    = log2(Mt);            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = log2(Mt); 
Nbits = 10000;
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
semilogy(SNR,ber,'-','LineWidth',1.5); hold on
clear ;   
clc;
M = 16;
Mt = 4; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nbits =20000;
SNR =0: 2 :24;                  % signal-to-noise ratio in dB    
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
semilogy(SNR,ber,'-','LineWidth',1.5);hold on

clear; clc
M = 8;  %the number of mudualtion constellation 
Mt = 11;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
Nu=5;   % the number of activated antennas
SNR = 0: 2 :24;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt = floor(log2(nchoosek(Mt,Nu)));            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = Nt; 
Nbits = 10000;
%generate bit stream
bit_T = randi([0 1],Nbits,Nt);
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;
%ssk mapping
x=zeros(Nbits,Mt);
l=(Mt-Nu+1);
%genrate mapping table
t=1;
for i=1:Mt-4
    for j=(i+1):Mt-3
        for k=(j+1):Mt-2
            for o=(k+1):Mt-1
                for u=(o+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            ta(t,4)=o;
            ta(t,5)=u;
            t=t+1;
                end
            end
        end
    end
end
%spatial modulation
for i=1:Nbits    
    x(i,ta(active_ant(i),:))=1/sqrt(Nu);   
end
%tranmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 

for ii=1:L_SNR 
   for j = 1 : size(x,1) 
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,2));
 x_detected=zeros(L_SNR,size(x,2));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt-4
           for kk=(jj+1):Mt-3
              for mm=(kk+1):Mt-2
                  for nn=(mm+1):Mt-1
                      for oo=(nn+1):Mt
            %x_detection(jj,kk,mm,nn,oo)=((channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)'+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
           %x_detection(jj,kk,mm,nn,oo)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk))*z(ll)).^2);
           h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+nn)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+oo);
            x_detection(jj,kk,mm,nn,oo)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*(1/sqrt(Nu))*h)'*(1/sqrt(Nu))*(h));
                      end
           end
        end
        end
        end
        x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(x_detection(:));
        [d,f,g,w,v] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g w v]
                 x_index(ii,j)=i-1;
                 break
             end
        end     
        %x_detection(x_detection==0)=inf;
        %[x_min,ind] = min(x_detection(:));
        %[d,f,g] = ind2sub(size(x_detection),ind);
       % x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        %p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f));
        %x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f))'*y(ii*Mr-Mr+1:ii*Mr,j)/p;
        %x_detected(ii,j)=mean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=harmmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)));
        %x_detected(ii,j)=trimmean(y(ii*Mr-Mr+1:ii*Mr,j)./channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+x_index(ii,j)),0.25);
     end   
end
%dem llr
 bit_Tdem=zeros(size(x,1),Nt*L_SNR);
for ii=1:L_SNR 
    for j = 1 : size(x,1) 
       bit_Tdem(j,ii*Nt-(Nt-1):ii*Nt)=  de2bi(x_index(ii,j),Nt,'left-msb');
    end
end

 %calculate ber
for ii=1:L_SNR
    a(ii)=length(find((bit_Tdem(:,ii*Nt-(Nt-1):ii*Nt)-bit_T)~=0));
    ber(ii)=a(ii)/(Nt*Nbits);
end
semilogy(SNR,ber,'-','LineWidth',1.5); hold on
 grid on;
ylim([10^(-5) 10^(0)])
xlabel('$$SNR$$','Interpreter','latex')
ylabel('BER','Interpreter','latex')
title('8 bits transmission BER of SM using i-MRC detector')
legend('SM 16QAM 16X4channel' ,'SM 32QAM 8X4channel' ,'SM 64QAM 4X4channel' ,'SSK  256X4channel' ,'MMSE-MIMO 4X4channel' ,'GSSK (11,5)  11X4channel' ,  'Location','NorthEast')    
toc