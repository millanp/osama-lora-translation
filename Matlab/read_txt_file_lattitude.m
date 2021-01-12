clear all
close all
clc;
SF = 8;
N = 2^SF;
DC = conj(sym_to_data_ang([1],N));
load('example.txt')
Rx_Buffer = example(:,1) + i*example(:,2);
Rx_Buffer = Rx_Buffer.';
stft(Rx_Buffer,N,DC(1:N),0,1);

%% write symbols to txt file, to be fed to RPP0
clear all
close all
clc
num_nodes = 15;
path = ['C:\temp\USRP_data\1min_15nodes_lab_2'];
for k = floor(num_nodes.*[0.5:0.5:5.5])%[1:15]%
    for iter = [1]
        filename = [path '\lamda_' num2str(k) '_' num2str(iter) '.mat'];
%         filename = [path '\' num2str(floor(k)) 'tx_lamda' num2str(iter) '.mat'];
        load(filename);
        fileID_1 = fopen([path '\sym_' num2str(k) '_' num2str(iter) '.txt'],'w');
%         fileID_1 = fopen(['C:\temp\USRP_data\sym_txt\sym_' num2str(k) 'tx_lm' num2str(iter) '.txt'],'w');
        %fwrite(fileID,demod_sym,'int16');
        demod_sym_stack(find(demod_sym_stack < 0)) = 0;
        for i = 1:size(demod_sym_stack,1)
            fprintf(fileID_1,'%d\n',-1);
            fprintf(fileID_1,'%d\n',demod_sym_stack(i,:));
        end
    end
    fclose(fileID_1);
end

%%  Compute BER + PER
clear all
clc
file_dur = 60;
BER_arr = [];
PER_arr = [];
% fileID = fopen('SF8_new/SF8_bits','r');
path = ['C:\temp\USRP_data\meeting-5-1-binout'];
load('C:\temp\USRP_data\meeting-5-1-binout\bits.mat')
% text = fileread('SF8_new/bin_out');
for k = 15*[0.5:0.5:5.5]%[1:15]%[10:5:60]%
for iter = 1
    filename = [path '\bin_' num2str(floor(k)) '_' num2str(iter) '.txt'];
%     filename = [path '\bin_' num2str(floor(k)) 'tx_lm' num2str(iter) '.txt'];
    fileID_1 = fopen(filename);

    n = 1;
    tline_1 = fgetl(fileID_1);
    out_bits = {};
    while ischar(tline_1)
        temp = [];
        for i = 1:length(tline_1)
            temp(i) = str2num(tline_1(i));
        end
        out_bits{n} = temp;
    %     out=char(bin2dec(num2str(reshape(temp,8,[])).'))'
        tline_1 = fgetl(fileID_1);
        n = n+1;
    end

    BER = 0;
    pkt = 0;
    for i = 1:size(out_bits,2)
        clear temp;
        temp = out_bits{i};
        if(length(temp) >= length(bits))
            if( sum(bits == temp(1:length(bits))) == length(bits) )
                pkt = pkt + 1;
            end
        end
        if(length(temp) >= length(bits))
            BER = BER + sum(temp(1:length(bits)) == bits);
        elseif(length(temp) < length(bits))
            BER = BER + sum(temp == bits(1:length(temp)));
        end
    end
%     PER_arr = [PER_arr (1 - (pkt/(k * file_dur)))];
    PER_arr = [PER_arr (1 - (pkt/(k * file_dur * iter)))];
    BER_arr = [BER_arr (1 - BER/(k * file_dur * length(bits)) )];
end
end
% BER = BER/size(out_bits,2);

% out=char(bin2dec(num2str(reshape(out_bits{3},8,[])).'))';
% out=char(bin2dec(num2str(reshape(bits,8,[])).'))';
%[0 0.0034 0.0109 0.0552 0.0385]
%%  CRC Calculation
clc
% pol_hex =  0x2A;%;
% pol_bin = hexToBinaryVector(0x2A);
hexStr = '0x1021';%'0xFFFF';%'0x8005';%'0x1D0F';%
pol_bin = [1 hexToBinaryVector(hexStr,16)]
% pol_bin = [1 0 1 1];
N = length(pol_bin);
bin_data_buf = bits;
% temp = [1 1 0 1 0 0 1 1 1 0 1 1 0 0 1     0     0];
temp = [bin_data_buf]
% temp = [bin_data_buf zeros(1,16)]
temp
for i = 1:length(temp) - N
    if(temp(i) == 0)
        continue;
    end
    temp(i : i+N-1) = xor(temp(i : i+N-1),pol_bin);
end
temp(end - N + 2:end)
%%  Calculating ser for poisson experiments
clear all
clc
ser_arr = [];
file_dur = 60;
num_sym = 28;
num_nodes = 15;
path = 'C:\temp\USRP_data\1min_15nodes_lab_2';
for i = num_nodes.*[0.5:0.5:5.5]%[1:15]%[10:5:60]%
    for iter = 1
    %     for i = [5:5:60]
        demod = [];
        load([path '\lamda_' num2str(floor(i)) '_' num2str(iter) '.mat']);
%         load([path '\' num2str(floor(i)) 'tx_lamda' num2str(iter) '.mat']);
        demod = [demod; demod_sym_stack];
    %     load([stream '/lamda_' num2str(i) '_2.mat']);
    %     demod = [demod; demod_sym_stack];
    %     load([stream '/lamda_' num2str(i) '_3.mat']);
    %     demod = [demod; demod_sym_stack];
    %     load([stream '/lamda_' num2str(i) '_4.mat']);
    %     demod = [demod; demod_sym_stack];
    %     load([stream '/lamda_' num2str(i) '_5.mat']);
    %     demod = [demod; demod_sym_stack];
        load([path '\sym.mat']);
        sym = sym(1:num_sym);
        sym = repmat(sym,size(demod,1),1);
%         tot = (((i)*file_dur*iter)*num_sym);
        tot = (((i)*file_dur)*num_sym);
        ser_arr = [ser_arr (1 - sum(sum((sym == demod)))/tot)];
    %     tot - sum(sum((sym == demod_sym_stack)))
    end
end
%%  Calculating Tx located %age for poisson experiments
clear all
clc
pkt_arr = [];
file_dur = 60;
num_nodes = 15;
path = 'C:\temp\USRP_data\1min_15nodes_lab_2';
% for i = [5:5:60]
 for i = num_nodes.*[0.5:0.5:5.5]%[1:15]%[10:5:60]%
         demod = [];
 for iter = 1
    load([path '\lamda_' num2str(floor(i)) '_' num2str(iter) '.mat']);
%     load([path '\' num2str(floor(i)) 'tx_lamda' num2str(iter) '.mat']);
    demod = [demod; demod_sym_stack];
%     load([stream '/lamda_' num2str(i) '_2.mat']);
%     demod = [demod; demod_sym_stack];
%     load([stream '/lamda_' num2str(i) '_3.mat']);
%     demod = [demod; demod_sym_stack];
%     load([stream '/lamda_' num2str(i) '_4.mat']);
%     demod = [demod; demod_sym_stack];
%     load([stream '/lamda_' num2str(i) '_5.mat']);
%     demod = [demod; demod_sym_stack];
%     tot_pkts = ((i)*file_dur*iter);
    tot_pkts = ((i)*file_dur);
    pkt_arr = [pkt_arr (size(demod,1)/tot_pkts)*100];
%     tot - sum(sum((sym == demod_sym_stack)))
 end
end
%%  Calculating PER from symbols for poisson experiments
clear all
clc
per_arr = [];
file_dur = 60;
num_data_sym = 28;
path = 'C:\temp\USRP_data\1min_SF8_tx';
% for i = [30:10:100]
for i = [1:15]%15*[0.5:0.5:6]%
        demod = [];
for iter = 6       
%     load([path '\lamda_' num2str(floor(i)) '_' num2str(iter) '.mat']);
    load([path '\' num2str(floor(i)) 'tx_lamda' num2str(iter) '.mat']);
    demod = [demod; demod_sym_stack];
    load([path '\sym.mat']);
    sym = repmat(sym,size(demod,1),1);
    tot = ((i*file_dur * iter));
%     tot = ((i*file_dur));
%     tot = ((((i/20) * 15)*file_dur)*28);
    v = sum((demod == sym),2);
    per_arr = [per_arr (1 - sum(v == num_data_sym)/tot)];
%     tot - sum(sum((sym == demod_sym_stack)))
end
end
%%  read symbols from RPP0 and save
clear all
clc
num_sym = 28;

path = ['C:\Osama\Matlab_lt\sym_bin\1min_15nodes_lab_sym_RPP0'];
for k = 15*[0.5:0.5:6]%[10:5:60]%
for iter = 1
    
    filename = [path '\sym_' num2str(floor(k)) '_' num2str(iter) 'lm.txt'];
    fileID_1 = fopen(filename);
    
    demod_sym_stack = [];
    n = 1;
    tline_1 = fgetl(fileID_1);
    temp = [];
    % A = fscanf(fileID,'%u\n')
    filename1 = [path '\lamda_' num2str(floor(k)) '_' num2str(iter) '.mat'];
    
    n = 1;
    while ischar(tline_1)
        tline_1 = fgetl(fileID_1);
        if(tline_1 == -1)
            break;
        end
        splitstring = regexp(tline_1,'\s+','split');
        if(size(splitstring,2) >= num_sym)
            for i = 1:28
        %         if(length(str2num(splitstring{i})) == 0)
        %             continue;
        %         end
                temp(n,i) = str2num(splitstring{i});

            end
        end
        n = n+1;
    end
    demod_sym_stack = temp;
    save(filename1,'demod_sym_stack');
end
end
