%%
pkt = 1;
ser = ~pkt;
% a1 = [0 0.0034 0.0109 0.0552 0.0385];
% a2 = [0    0.0116    0.0379    0.1118    0.0744];
% a = [a1; a2];

% a6 = [0.05 0.18 0.318 0.39];
% a5 = [0.39 0.2112 0.1619 0.0926 0.0735];
% a4 = [0.318 0.2012 0.1192 0.0994 0.0895];
% a3 = [0.18 0.1047 0.0519 0.0256 0.0256];
% a2 = [0.05 0.0233 0.0442 0.0093 0.0093];
% a = [a2; a3; a4; a5];

% a1 = [0.3591 0.3009 0.2656 0.2391 0.2153 0.2014 0.1893 0.1809 0.1865 0.1888]';
% a2 = [0.2051 0.1567 0.1242 0.1181 0.1065 0.1242 0.1340 0.1456]';
% a3 = [0.1963 0.1637 0.1665 0.1837 0.1847 0.1902 0.1953 0.1977]';
% 
% a4 = [0.1612 0.0977 0.0659 0.0473 0.0426 0.0419 0.0426 0.0481]';
% a5 = [0.1558 0.1221 0.1291 0.1233 0.1297 0.1279]';
% b = bar([5:5:60],a,'grouped');
b = bar(15.*[0.5:0.5:6],a,'grouped');
% 
set(gca,'linewidth',1.5,'fontsize',25,'fontname','Times New Roman');
% xlabel('x - value','FontSize',30);
% xlabel('# of TX collision','FontSize',30);
xlabel('Aggregate Rate (Pkts/sec)','FontSize',30);
% xlabel('Number of Nodes','FontSize',30);


% view(0,90);
set(gca,'YDir','normal');
set(gcf,'Color','w');
legend('CT','FTRACK')
% legend('SNR = 20dB','SNR = 15dB','SNR = 10dB','SNR = 5dB','SNR = 0dB')
% legend('SNR = 20dB')
% legend('FTrack','CT','CT + offset-filter','CT + PWR-filter','CT + offset-filter + PWR-filter')
% legend('SNR = 20dB','SNR = 15dB','SNR = 10dB','SNR = 5dB')
%
if(pkt == 1)
    ylim([0 100])
    title('bargraph of number of Pkts located','FontSize',30);
    ylabel('Packets located (%age)','FontSize',30);
elseif(ser == 1)
    ylim([0 1])
     title('bargraph of SER','FontSize',30);
     ylabel('SER','FontSize',30);
end
% pbaspect([1 1 1])
grid minor