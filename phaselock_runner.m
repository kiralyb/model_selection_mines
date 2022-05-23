N = 8;
d = [1:0.2:2];
for i = 4:length(d)
[dAIC(i),dBIC(i)] = phaselock_modes_simulation_param(N(1),d(i));
end
figure
scatter ([2:0.2:3],flip(dAIC(1:6)),'g')
hold on
scatter ([2:0.2:3],flip(dBIC(1:6)),'b')
line([2,3],[0,0])



%%

N = [1,11,21,31,41];
d = [1:0.1:2];
for i = 1:length(N)
[dAIC(i),dBIC(i)] = phaselock_modes_simulation_param(N(i),1.15);
end
figure
scatter (250*[1:10:41],dAIC(1:5),'g')
hold on
scatter (250*[1:10:41],dBIC(1:5),'b')
line(250*[1,41],[0,0])
