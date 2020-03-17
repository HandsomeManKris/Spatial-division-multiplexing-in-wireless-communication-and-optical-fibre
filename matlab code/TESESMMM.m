 clear;   
clc;
M = 16;
Mt = 16;                         %for number of transmit antenna
Mr = 4;                         %for number of receive antenna
bit_SMsym = log2(M*Mt);         % number of bit per spatial modulation sysmbol
Nbits =100;
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
z=constellation(hMod);
aa=reshape(z,M,1);
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
             for ll=1:M
            %x_detection(jj)=(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)'*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj);
            x_detection(jj,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);
             end
        end
        [x_max,ind] = max(x_detection(:));
        [d,f] = ind2sub(size(x_detection),ind);
       % r=max(abs(x_detection));
        x_index(ii,j)=d; 
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
