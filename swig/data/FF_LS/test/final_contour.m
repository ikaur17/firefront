fy = fopen('/home/ikaur/firefront/swig/data/FF_LS/final_points_100.txt');
B = textscan(fy, '%f %f %f');
X = reshape(B{1}, [], 1);
Y = reshape(B{2}, [], 1);
Z = reshape(B{3}, [], 1);
p = [X, Y];
t = [X, Y, Z];
% tricontour(p, t, Z,  4);
%  gscatter(X, Y, Z);
% % line(X, Y)
% xlabel('X (m)')
% ylabel('Y (m)')
% grid on;
% legend('0', '20', '40', '60', '80', '100', '120', '140')
% legend boxoff;
% % axis([1500, 6500, 1500, 6500]);