        
tic
clear all;   
clc;
M = 4;
Mt = 4; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=2;
Nk=10;
Nbits =50000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb');  %generate active antenna number
bit_T=bit_T.';
x=zeros(Nbits,Mt);
hmodem =   modem.qammod('M',M,  'SymbolOrder', 'Gray','InputType', 'bit');
AP = (mean(hmodem.Constellation .* conj(hmodem.Constellation)));
No=10.^(-SNR/10);               % noise variance
hMod = comm.RectangularQAMModulator('ModulationOrder',M,'BitInput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.RectangularQAMDemodulator('ModulationOrder',M,'BitOutput',true,'NormalizationMethod','Average Power','SymbolMapping','Gray','OutputDataType', 'Full precision' );
z=constellation(hMod);
aa=reshape(z,M,1);
l=(Mt-Nu+1);
for i=1:Nbits
Hmod(i) = step(hMod,bit_T(Nt+1:end,i));
end

%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits
    a=0;
    t=1;
    for ii=l:-1:1
    a=a+ii;
        if active_ant(i)<a
        x(i,t)=Hmod(i);
        x(i,t+active_ant(i)-(l+l-t+2)*(t-1)/2+1)=Hmod(i);
        break
        end
         t=t+1;
    end
end
bit_T=bit_T.';
%transmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
z=[randn(Mt)+1i*randn(Mt)]/sqrt(2);
[Q,R] = qr(z);
d = diag(R);
ph = diag(d./abs(d));
Q=Q*ph*Q;
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) =Q;
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii))*sqrt(Nk);
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end
%detection
x_detection=zeros(Mt,Mt,M);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Nt
           for kk=(jj+1):(Mt-floor(jj/Nt))
             %  for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk);
           x_detection(jj,kk)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           %x_detection(jj,kk,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
           
            %   end
           end
        end
       % x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(abs(x_detection(:)));
        %[d,f,g] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        %[x_max,ind] = max(x_detection(:));
        [d,f] = ind2sub(size(x_detection),ind);
        x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
       % p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)')*y(ii*Mr-Mr+1:ii*Mr,j)/2;
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
     bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt)= de2bi(x_index(ii,j),Nt,'left-msb');
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
clear all;   
clc;
M = 2;
Mt = 5; %for number of transmit antenna
Mr = 5;%for number of receive antenna
Nu=2;
Nk=10;
Nbits =50000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb');  %generate active antenna number
bit_T=bit_T.';
x=zeros(Nbits,Mt);
hmodem =   modem.qammod('M',M,  'SymbolOrder', 'Gray','InputType', 'bit');
AP = (mean(hmodem.Constellation .* conj(hmodem.Constellation)));
No=10.^(-SNR/10);               % noise variance
hMod = comm.PSKModulator('ModulationOrder',M,'BitInput',true,'SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.PSKDemodulator('ModulationOrder',M,'BitOutput',true,'SymbolMapping','Gray','OutputDataType', 'Full precision' );
z=constellation(hMod);
aa=reshape(z,M,1);
l=(Mt-Nu+1);
for i=1:Nbits
Hmod(i) = step(hMod,bit_T(Nt+1:end,i));
end

%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits
    a=0;
    t=1;
    for ii=l:-1:1
    a=a+ii;
        if active_ant(i)<a
        x(i,t)=Hmod(i);
        x(i,t+active_ant(i)-(l+l-t+2)*(t-1)/2+1)=Hmod(i);
        break
        end
         t=t+1;
    end
end
bit_T=bit_T.';
%transmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
z=[randn(Mt)+1i*randn(Mt)]/sqrt(2);
[Q,R] = qr(z);
d = diag(R);
ph = diag(d./abs(d));
Q=Q*ph*Q;
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) =Q;
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii))*sqrt(Nk);
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end
%detection
x_detection=zeros(Mt,Mt,M);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Nt
           for kk=(jj+1):(Mt-floor(jj/Nt))
             %  for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk);
           x_detection(jj,kk)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           %x_detection(jj,kk,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
           
            %   end
           end
        end
       % x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(abs(x_detection(:)));
        %[d,f,g] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        %[x_max,ind] = max(x_detection(:));
        [d,f] = ind2sub(size(x_detection),ind);
        x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
       % p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)')*y(ii*Mr-Mr+1:ii*Mr,j)/2;
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
     bit_Tdem1(j,ii*Nt-(Nt-1):ii*Nt)= de2bi(x_index(ii,j),Nt,'left-msb');
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
M = 8;
Mt = 2;                         %for number of transmit antenna
Mr = 2;     
Nk=10;%for number of receive antenna
bit_SMsym = log2(M*Mt);         % number of bit per spatial modulation sysmbol
Nbits =50000;
SNR = 0 : 2 :26;                % signal-to-noise ratio in dB                          
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
hMod = comm.PSKModulator('ModulationOrder',M,'BitInput',true,'SymbolMapping','Gray','OutputDataType','double');
hdeMod=comm.PSKDemodulator('ModulationOrder',M,'BitOutput',true,'SymbolMapping','Gray','OutputDataType', 'Full precision' );
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
    z=[randn(Mt)+1i*randn(Mt)]/sqrt(2);
[Q,R] = qr(z);
d = diag(R);
ph = diag(d./abs(d));
Q=Q*ph*Q;
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) =Q;
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii))*sqrt(Nk);
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mr-1):j*Mr)*(x(j,:).')+noise;      
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
            
            x_detection(jj)=(channel(ii*Mt-(Mt-1):ii*Mt,j*Mt-Mt+jj)'*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mt-(Mt-1):ii*Mt,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
        end
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
M = 4;
Mt = 4;                         %for number of transmit antenna
Mr = 4;     
Nk=10;%for number of receive antenna
bit_SMsym = log2(M*Mt);         % number of bit per spatial modulation sysmbol
Nbits =50000;
SNR = 0 : 2 :26;                % signal-to-noise ratio in dB                          
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
    z=[randn(Mt)+1i*randn(Mt)]/sqrt(2);
[Q,R] = qr(z);
d = diag(R);
ph = diag(d./abs(d));
Q=Q*ph*Q;
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) =Q;
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii))*sqrt(Nk);
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mr-1):j*Mr)*(x(j,:).')+noise;      
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
            
            x_detection(jj)=(channel(ii*Mt-(Mt-1):ii*Mt,j*Mt-Mt+jj)'*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mt-(Mt-1):ii*Mt,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
        end
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
 
 clear ;   
clc;
M = 4;
Nk=10;
Mt = 2; %for number of transmit antenna
Mr = 2;%for number of receive antenna
Nbits =50000;
SNR =0: 2 :26;                  % signal-to-noise ratio in dB    
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
z=[randn(Mt)+1i*randn(Mt)]/sqrt(2);
[Q,R] = qr(z);
d = diag(R);
ph = diag(d./abs(d));
Q=Q*ph*Q;
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) =Q;
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii))*sqrt(Nk);
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end


%detection
x_detected=zeros(L_SNR,size(x,1));
for ii=1:L_SNR 
   for j = 1 : Nbits
       x_detected(ii*Mr-(Mr-1):ii*Mr,j)=(eye(Mr,Mt)+No(ii)*eye(Mr,Mt))\channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)'*y(ii*Mr-(Mr-1):ii*Mr,j);
   end
   x_detect(ii,:)=reshape(x_detected(ii*Mr-(Mr-1):ii*Mr,:),1,Nbits*Mt);
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
semilogy(SNR,ber,'--','LineWidth',1.5);hold on
clear;clc
M = 4;
Mt = 2; %for number of transmit antenna
Mr = 2;%for number of receive antenna
Nbits =50000;
Nk=10;
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
z=[randn(Mt)+1i*randn(Mt)]/sqrt(2);
[Q,R] = qr(z);
d = diag(R);
ph = diag(d./abs(d));
Q=Q*ph*Q;
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) =Q;
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii))*sqrt(Nk);
      y(ii*Mr-(Mr-1):ii*Mr,j)=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)*(x(j,:).')+noise; 
   end
end


%detection
x_detected=zeros(L_SNR,size(x,1));
for ii=1:L_SNR 
   for j = 1 : Nbits
       x_detected(ii*Mr-(Mr-1):ii*Mr,j)=(eye(Mr,Mt))\channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt)'*y(ii*Mr-(Mr-1):ii*Mr,j);
   end
   x_detect(ii,:)=reshape(x_detected(ii*Mr-(Mr-1):ii*Mr,:),1,Nbits*Mt);
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
semilogy(SNR,ber,'-.','LineWidth',1.5);hold on
grid on;
xlabel('$$SNR$$','Interpreter','latex')
ylabel('BER','Interpreter','latex')
title('Receiver complexity 4bits transmission  in MMF Nk=10')
legend('MMF GSM 4QAM (4,2) ','MMF GSM BPSK (5,2) ','MMF SM 8PSK 2x2 ','MMF SM 4QAM 4X4 ' ,'MMF MMSE-MIMO 4QAM 4X4 ','MMF ZF-MIMO 16QAM 2X2 ','Location','NorthEast')
toc