tic
clear; clc
M = 8;  %the number of mudualtion constellation 
Mt = 12;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
Nu=2;   % the number of activated antennas
SNR = 0: 2 :24;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt = floor(log2(nchoosek(Mt,Nu)));            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = Nt; 
Nbits = 100;
%generate bit stream
bit_T = randi([0 1],Nbits,Nt);
active_ant=bi2de(bit_T(:,1:Nt),'left-msb');
%ssk mapping
x=zeros(Nbits,Mt);
l=(Mt-Nu+1);
for i=1:Nbits
    a=0;
    t=1;
    for ii=l:-1:1
    a=a+ii;
        if active_ant(i)<a
        x(i,t)=1/sqrt(Nu);
        x(i,t+active_ant(i)-(l+l-t+2)*(t-1)/2+1)=1/sqrt(Nu);
        break
        end
         t=t+1;
    end
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
 y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*x(j,:).'+noise;      
   end
end
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,2));
 x_detected=zeros(L_SNR,size(x,2));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:10
           for kk=(jj+1):(Mt-floor(jj/10))
              h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk);
            %x_detection(jj,kk)=(h')*y(ii*Mr-(Mr-1):ii*Mr,j)/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j)));
           %x_detection(jj,kk)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk))*z(ll)).^2);
            x_detection(jj,kk)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*(1/sqrt(Nu))*(h))'*(1/sqrt(Nu))*(h));
            
           end
        end
   
        x_detection(x_detection==0)=-inf;
        %[x_max,ind] = max(abs((x_detection(:))));
        [x_max,ind] = max(x_detection(:));
        [d,f] = ind2sub(size(x_detection),ind);
        x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2; 
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
semilogy(SNR,ber,'-','LineWidth',1); hold on
clear; clc
M = 8;  %the number of mudualtion constellation 
Mt = 12;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
Nu=2;   % the number of activated antennas
SNR = 0: 2 :24;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt = floor(log2(nchoosek(Mt,Nu)));            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = Nt; 
Nbits = 100;
%generate bit stream
bit_T = randi([0 1],Nbits,Nt);
active_ant=bi2de(bit_T(:,1:Nt),'left-msb');
%ssk mapping
x=zeros(Nbits,Mt);
l=(Mt-Nu+1);
for i=1:Nbits
    a=0;
    t=1;
    for ii=l:-1:1
    a=a+ii;
        if active_ant(i)<a
        x(i,t)=1/sqrt(Nu);
        x(i,t+active_ant(i)-(l+l-t+2)*(t-1)/2+1)=1/sqrt(Nu);
        break
        end
         t=t+1;
    end
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
 y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*x(j,:).'+noise;      
   end
end
x_detection=zeros(Mr,1);  
x_index=zeros(L_SNR,size(x,2));
 x_detected=zeros(L_SNR,size(x,2));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:10
           for kk=(jj+1):(Mt-floor(jj/10))
              h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk);
            x_detection(jj,kk)=(h')*y(ii*Mr-(Mr-1):ii*Mr,j)/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j)));
           %x_detection(jj,kk)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk))*z(ll)).^2);
           % x_detection(jj,kk)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*(1/sqrt(Nu))*(h))'*(1/sqrt(Nu))*(h));
            
           end
        end
   
        %x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(abs((x_detection(:))));
        %[x_max,ind] = max(x_detection(:));
        [d,f] = ind2sub(size(x_detection),ind);
        x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2; 
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
semilogy(SNR,ber,'--','LineWidth',1); hold on
grid on;
xlabel('$$SNR$$','Interpreter','latex')
ylabel('BER','Interpreter','latex')
title('Bit Error Rate of GSSK using ML and i-MRC detector 2 active antennas')
legend( '(5,2) ML 3bits','(5,2)  MRC 3bits','(7,2)  ML 4bits','(7,2)  MRC 4bits', '(12,2)  ML 6bits','(12,2)  MCR 6bits','MMSE-MIMO 3X4 M=2 3bits','MMSE-MIMO 2X4 M=4 4bits', 'MMSE-MIMO 3X4 M=4','SSK 3bits','SSK 5bits','SSK 6bits','Location','NorthEast')
toc
     
