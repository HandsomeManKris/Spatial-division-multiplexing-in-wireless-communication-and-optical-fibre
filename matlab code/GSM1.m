tic
clear;   
clc;
M = 16;
Mt = 5; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=2;

Nbits =10000;
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
 channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
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
               for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk);
           %x_detection(jj,kk,ll)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           x_detection(jj,kk,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
           
               end
           end
        end
        x_detection(x_detection==0)=-inf;
        %[x_max,ind] = max(abs(x_detection(:)));
        %[d,f,g] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        [x_max,ind] = max(x_detection(:));
        [d,f,g] = ind2sub(size(x_detection),ind);
        x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'-ko','LineWidth',1.5);hold on

 clear ;   
clc;
M = 16;
Mt = 5; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=3;

Nbits =10000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
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
%genrate mapping table
t=1;
for i=1:Mt-2
    for j=(i+1):Mt-1
        for k=(j+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            t=t+1;
        end
    end
end
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits    
    x(i,ta(active_ant(i),:))=Hmod(i);   
end
bit_T=bit_T.';
%transmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt) + 1i*randn(Mr,Mt));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
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
        for jj=1:Mt-2
           for kk=(jj+1):Mt-1
               for mm=(kk+1):Mt
               for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm);
           %x_detection(jj,kk,mm,ll)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,mm,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           x_detection(jj,kk,mm,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
           
               end
           end
           end
        end
        x_detection(x_detection==0)=-inf;
        %[x_max,ind] = max(abs(x_detection(:)));
        %[d,f,g] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        [x_max,ind] = max(x_detection(:));
        [d,f,g,b] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g]
                 x_index(ii,j)=i-1;
                 break
             end
        end    
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'-gs','LineWidth',1.5);hold on
 clear ;   
clc;
M = 4;
Mt = 7; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=4;

Nbits =10000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
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
%genrate mapping table
t=1;
for i=1:Mt-3
    for j=(i+1):Mt-2
        for k=(j+1):Mt-1
              for p= (k+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            ta(t,4)=p;
            t=t+1;
              end
        end
    end
end
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits    
    x(i,ta(active_ant(i),:))=Hmod(i);   
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
x_detection=zeros(Mt,Mt,M);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt-3
           for kk=(jj+1):Mt-2
               for mm=(kk+1):Mt-1
                   for nn=(mm+1):Mt
               for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+nn);
           %x_detection(jj,kk,mm,nn,ll)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,mm,nn,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           x_detection(jj,kk,mm,nn,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
               end
               end
           end
           end
        end
        x_detection(x_detection==0)=-inf;
        %[x_max,ind] = max(abs(x_detection(:)));
        %[d,f,g] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        [x_max,ind] = max(x_detection(:));
        [d,f,g,u,b] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g u]
                 x_index(ii,j)=i-1;
                 break
             end
        end    
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'-bd','LineWidth',1.5);hold on
clear ;   
clc;
M = 4;
Mt = 8; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=5;

Nbits =10000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
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
%genrate mapping table
t=1;
for i=1:Mt-4
    for j=(i+1):Mt-3
        for k=(j+1):Mt-2
              for p= (k+1):Mt-1
                  for o=(p+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            ta(t,4)=p;
            ta(t,5)=o;
            t=t+1;
                  end
              end
        end
    end
end
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits    
    x(i,ta(active_ant(i),:))=Hmod(i);   
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
x_detection=zeros(Mt,Mt,M);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt-4
           for kk=(jj+1):Mt-3
               for mm=(kk+1):Mt-2
                   for nn=(mm+1):Mt-1
                       for oo=(nn+1):Mt
               for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+nn)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+oo);
           %x_detection(jj,kk,mm,nn,oo,ll)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,mm,nn,oo,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           x_detection(jj,kk,mm,nn,oo,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
               end
                       end
               end
           end
           end
        end
        x_detection(x_detection==0)=-inf;
        %[x_max,ind] = max(abs(x_detection(:)));
        %[d,f,g] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        [x_max,ind] = max(x_detection(:));
        [d,f,g,u,r,b] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g u r]
                 x_index(ii,j)=i-1;
                 break
             end
        end    
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+r))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+r));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+r)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'-mv','LineWidth',1.5);hold on
 clear;   
clc;
M = 16;
Mt = 5; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=2;

Nbits =10000;
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
 channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt,1) + 1i*randn(Mr,Mt,1));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
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
              % for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk);
           x_detection(jj,kk)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           %x_detection(jj,kk,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
           
             %  end
           end
        end
        %x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(abs(x_detection(:)));
        [d,f] = ind2sub(size(x_detection),ind);
        x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
        %[x_max,ind] = max(x_detection(:));
       % [d,f,g] = ind2sub(size(x_detection),ind);
       % x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'--ko','LineWidth',1.5);hold on

clear ;   
clc;
M = 16;
Mt = 5; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=3;

Nbits =10000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
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
%genrate mapping table
t=1;
for i=1:Mt-2
    for j=(i+1):Mt-1
        for k=(j+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            t=t+1;
        end
    end
end
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits    
    x(i,ta(active_ant(i),:))=Hmod(i);   
end
bit_T=bit_T.';
%transmission
L_SNR=length(SNR);
ber= zeros (L_SNR,1);
y=zeros(L_SNR*Mr,size(x,1));
No= 10.^(-SNR/10); 
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-(Mt-1):j*Mt) = sqrt(.5)*( randn(Mr,Mt) + 1i*randn(Mr,Mt));
     noise = sqrt(.5)*(randn(Mr,1) + 1i*randn(Mr,1))* sqrt(No(ii));
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
        for jj=1:Mt-2
           for kk=(jj+1):Mt-1
               for mm=(kk+1):Mt
             %  for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm);
           x_detection(jj,kk,mm)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,mm,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           %x_detection(jj,kk,mm,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
           
             %  end
           end
           end
        end
        %x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(abs(x_detection(:)));
        [d,f,g] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
       % [x_max,ind] = max(x_detection(:));
       %[d,f,g,u] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g]
                 x_index(ii,j)=i-1;
                 break
             end
        end    
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'--gs','LineWidth',1.5);hold on
 clear ;   
clc;
M = 4;
Mt = 7; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=4;

Nbits =10000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
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
%genrate mapping table
t=1;
for i=1:Mt-3
    for j=(i+1):Mt-2
        for k=(j+1):Mt-1
              for p= (k+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            ta(t,4)=p;
            t=t+1;
              end
        end
    end
end
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits    
    x(i,ta(active_ant(i),:))=Hmod(i);   
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
x_detection=zeros(Mt,Mt,M);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt-3
           for kk=(jj+1):Mt-2
               for mm=(kk+1):Mt-1
                   for nn=(mm+1):Mt
             %  for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+nn);
           x_detection(jj,kk,mm,nn)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,mm,nn,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           %x_detection(jj,kk,mm,nn,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
             %  end
               end
           end
           end
        end
       % x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(abs(x_detection(:)));
        [d,f,g,u] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
       % [x_max,ind] = max(x_detection(:));
       % [d,f,g,u,b] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g u]
                 x_index(ii,j)=i-1;
                 break
             end
        end    
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'--bd','LineWidth',1.5);hold on
clear ;   
clc;
M = 4;
Mt = 8; %for number of transmit antenna
Mr = 4;%for number of receive antenna
Nu=5;

Nbits =10000;
SNR =0: 2 :24;                % signal-to-noise ratio in dB                          
L_SNR=length(SNR);              %lentgh of ber picture
ber= zeros (L_SNR,1);           %for plotting ber
Nt = floor(log2(nchoosek(Mt,Nu))); %bits for transmit antenna
bit_SMsym = log2(M)+Nt;         % number of bit per spatial modulation sysmbol
Numberofbit = log2(M);          %bits for diffierentkindsof modulation  
bit_T1 = randi([0 1],Nbits,Nt); %randomly choose transmita antenna
bit_T2 = randi([0 1],Nbits,Numberofbit);  %transmit bits
bit_T=[bit_T1 bit_T2];%get encoded sysmbol for transmiision
active_ant=bi2de(bit_T(:,1:Nt),'left-msb')+1;  %generate active antenna number
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
%genrate mapping table
t=1;
for i=1:Mt-4
    for j=(i+1):Mt-3
        for k=(j+1):Mt-2
              for p= (k+1):Mt-1
                  for o=(p+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            ta(t,4)=p;
            ta(t,5)=o;
            t=t+1;
                  end
              end
        end
    end
end
%spatial modulation
x=zeros(Nbits,Mt);
for i=1:Nbits    
    x(i,ta(active_ant(i),:))=Hmod(i);   
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
x_detection=zeros(Mt,Mt,M);  
x_index=zeros(L_SNR,size(x,1));
x_detected=zeros(L_SNR,size(x,1));
%detection
for ii=1:L_SNR 
   for j = 1 : size(x,1) 
        for jj=1:Mt-4
           for kk=(jj+1):Mt-3
               for mm=(kk+1):Mt-2
                   for nn=(mm+1):Mt-1
                       for oo=(nn+1):Mt
              % for ll=1:M
            h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+nn)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+oo);
           x_detection(jj,kk,mm,nn,oo)=((h')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(h)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))); % i-MRC detector
           %x_detection(jj,kk,mm,nn,oo,ll)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-h*z(ll)).^2);                            %MLdetector
           %x_detection(jj,kk,mm,nn,oo,ll)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*aa(ll)*h)'*aa(ll)*h);                  %MLdetector
             %  end
                       end
               end
           end
           end
        end
        %x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(abs(x_detection(:)));
        [d,f,g,u,r] = ind2sub(size(x_detection),ind);
        %x_index(ii,j)=f-d-1+(l+l-d+2)*(d-1)/2;
        %x_detected(ii,j)=z(g);
       % [x_max,ind] = max(x_detection(:));
       % [d,f,g,u,r,b] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g u r]
                 x_index(ii,j)=i-1;
                 break
             end
        end    
        p=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+r))'*(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+r));
        %x_detected(ii,j)=z(g);
        x_detected(ii,j)=(channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+d)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+f)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+g)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+u)'+channel(ii*Mr-Mr+1:ii*Mr,j*Mt-Mt+r)')*y(ii*Mr-Mr+1:ii*Mr,j)/p;
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

 semilogy(SNR,ber,'--mv','LineWidth',1.5);hold on
grid on;
ylim([10^(-5) 10^(0)])
xlabel('$$SNR$$','Interpreter','latex')
ylabel('BER','Interpreter','latex')
title('7 bits transmission BER of GSM using ML detector and i-MCR detector')
legend('GSM (5,2) ML' ,'GSM (5,3) ML', 'GSM (7,4) ML','GSM (8,5)ML', 'GSM (5,2) i-MRC' ,'GSM (5,3) i-MRC', 'GSM (7,4) i-MRC','GSM (8,5) i-MRC','Location','NorthEast')       
 toc       
     
     
