kp= 0;
for i = 1:16
    subplot(8,2,kp+1)
    plot(video_raw_data(:,i))
    title(['Noisy signal ',num2str(i)])
    kp = kp + 1;
end

level = 4;
wname = 'db4';
npc = 'heur';
[x_sim, qual, npc] = wmspca(video_raw_data,level,wname,npc);
qual

figure(2)
kp = 0;
for i = 1:16
    subplot(8,2,kp+1)
    plot(x_sim(:,i))
    title(['First PCA ',num2str(i)])
    kp = kp+1;
end