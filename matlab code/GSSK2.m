tic
clear; clc
M = 4;  %the number of mudualtion constellation 
Mt = 24;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
Nu=2;   % the number of activated antennas
SNR = 0: 2 :24;          %the snr requirement
No= 10.^(-SNR/10);     % noise variance  for the requirement of snr     
L_SNR=length(SNR);
Nt = floor(log2(nchoosek(Mt,Nu)));            %number of bit of tranmit antenna
Nobit = log2(M);             %number of bit of tranmsit symbol
bit_SMsym = Nt; 
Nbits = 10000;
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
        for jj=1:18
           for kk=(jj+1):(Mt-floor(jj/18)*5)
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
semilogy(SNR,ber,'-','LineWidth',1.5); hold on

clear; clc
M = 8;  %the number of mudualtion constellation 
Mt = 13;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
Nu=3;   % the number of activated antennas
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
        for jj=1:Mt-2
           for kk=(jj+1):Mt-1
              for mm=(kk+1):Mt
            %x_detection(jj,kk,mm)=((channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)'+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
           %x_detection(jj,kk,mm)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk))*z(ll)).^2);
           h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm);
            x_detection(jj,kk,mm)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*(1/sqrt(Nu))*h)'*(1/sqrt(Nu))*(h));
            
           end
        end
        end
   
        x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(x_detection(:));
        [d,f,g] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g]
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
clear; clc
M = 8;  %the number of mudualtion constellation 
Mt = 11;% the nubmer of transmit antenna
Mr = 4;% the number of receive antenna
Nu=4;   % the number of activated antennas
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
for i=1:Mt-3
    for j=(i+1):Mt-2
        for k=(j+1):Mt-1
            for o=(k+1):Mt
            ta(t,1)=i;
            ta(t,2)=j;
            ta(t,3)=k;
            ta(t,4)=o;
            t=t+1;
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
        for jj=1:Mt-3
           for kk=(jj+1):Mt-2
              for mm=(kk+1):Mt-1
                  for nn=(mm+1):Mt
            %x_detection(jj,kk,mm,nn)=((channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)'+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)')*y(ii*Mr-(Mr-1):ii*Mr,j))/(norm(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)*norm(y(ii*Mr-(Mr-1):ii*Mr,j))));
           %x_detection(jj,kk,mm,nn)=sum(abs(y(ii*Mr-(Mr-1):ii*Mr,j)-(channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk))*z(ll)).^2);
           h=channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+jj)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+kk)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+mm)+channel(ii*Mr-(Mr-1):ii*Mr,j*Mt-Mt+nn);
            x_detection(jj,kk,mm,nn)=real((y(ii*Mr-(Mr-1):ii*Mr,j)-(1/2)*(1/sqrt(Nu))*h)'*(1/sqrt(Nu))*(h));
            
           end
        end
        end
        end
        x_detection(x_detection==0)=-inf;
        [x_max,ind] = max(x_detection(:));
        [d,f,g,w] = ind2sub(size(x_detection),ind);
        for i=1:2^Nt
             if ta(i,:)==[d f g w]
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
clear;clc
M = 4;
Mt = 4; %for number of transmit antenna
Mr =4;%for number of receive antenna
Nbits =50000;
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
semilogy(SNR,ber,'--','LineWidth',1.5);hold on
clear;clc
M = 16;
Mt = 2; %for number of transmit antenna
Mr =4;%for number of receive antenna
Nbits =10000;
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
semilogy(SNR,ber,'--','LineWidth',1.5);hold on
grid on;
xlabel('$$SNR$$','Interpreter','latex')
ylabel('BER','Interpreter','latex')
title('8bits transmission BER of GSSK using ML detector with different active antennas ')
legend('GSSK (24,2)' , 'GSSK (13,3)' ,'GSSK (11,4)' ,'GSSK (11,5)' , 'MMSE-MIMO M=4 4X4' ,'MMSE-MIMO M=16 2X4' ,'Location','NorthEast')
toc