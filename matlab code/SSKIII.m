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