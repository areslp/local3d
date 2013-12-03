clear all;
close all;

pnum=77678;

ours=[20 139 0.0521995;
    30 148 0.0500138;
    40 138 0.0434336;
    50 135 0.0426198;
    60 130 0.0420869;
    70 137 0.041386;
    80 145 0.0409924;
    90 146 0.0402313;
    100 151 0.0410302];

lrr=[20 316 0.0381556;
    30 414 0.0373343;
    40 677 0.0384946;
    50 828 0.0370635;
    60 1036 0.0374528;
    70 1393 0.0379977;
    80 1691 0.0374946;
    90 2245 0.0377065;
    100 2800 0.0378641];

linear_lrr=[20 204 0.0523885;
    30 232 0.0497409;
    40 306 0.0437856;
    50 414 0.0423275;
    60 507 0.0421078;
    70 666 0.0416793;
    80 856 0.0409824;
    90 1156 0.04054
    100 1615 0.0405087];

x=[20 30 40 50 60 70 80 90 100];
% plot time
% figure;
% title('Fandisk model, 77678 points, 0.5 noise');
% hold on;

% y=ours(:,2);
% y=y/pnum;
% plot(x,y,'--rs');

% y=lrr(:,2);
% y=y/pnum;
% plot(x,y,'-g+');

% y=linear_lrr(:,2);
% y=y/pnum;
% plot(x,y,':bx');
% hold off;
% axis tight;
% set(gcf, 'Color', 'w');
% xlabel('k');
% ylabel('time per point(second)');

% plot acc
% figure;
% title('Fandisk model, 77678 points, 0.5 noise');
% hold on;

% y=ours(:,3);
% plot(x,y,'--rs');

% y=lrr(:,3);
% plot(x,y,'-g+');

% y=linear_lrr(:,3);
% plot(x,y,':bx');
% hold off;
% axis tight;
% set(gcf, 'Color', 'w');
% xlabel('k');
% ylabel('accuracy');
% set(gca, 'ylim', [0, 0.08]);

acc_per_lambda=[1.0 0.0459191;
                0.9 0.0459023;
                0.8 0.0457552;
                0.7 0.0456338;
                0.6 0.0453701;
                0.5 0.0450245;
                0.4 0.0447262;
                0.3 0.0445822;
                0.2 0.0447109;
                0.1 0.0475806];

x=acc_per_lambda(:,1);
figure;
title('Fandisk model, 77678 points, 0.5 noise');
hold on;

y=acc_per_lambda(:,2);
plot(x,y,'--rs');
axis tight;
set(gcf, 'Color', 'w');
xlabel('lambda');
ylabel('accuracy');
set(gca, 'ylim', [0, 0.1]);
