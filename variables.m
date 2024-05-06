%% Look into:
%https://www.zora.uzh.ch/id/eprint/111164/1/ICRA15_Faessler.pdf
%https://www.youtube.com/watch?v=pGU1s6Y55JI

%% Generate Model
% Dumbass, don't forget gravity!
g = 9.81;

motor_thrust = 5;
% Individual arm mass and total mass
mass_arm = 0.18; % kg (180 grams)
total_mass = 1; % kg

% Distance from center of mass to motor (assuming equal for all arms)
arm_length = 0.2; % meters (20 cm)

% Tilt Force Length
tilt_force_length = arm_length*sqrt(2);

% Inertia of each arm (assuming uniform rod along arm length)
I_arm = (mass_arm * arm_length^2) / 3; % kg*m^2

% Center mass inertia (assuming point mass)
I_center = (total_mass - 4 * mass_arm) * 0; % kg*m^2 (assuming negligible inertia for center mass)

% Offsets due to mass distribution (assuming all masses are at the end of arms)
offsetX = arm_length * mass_arm; % meters
offsetY = arm_length * mass_arm; % meters

% Product of inertia terms due to mass distribution (considering x and y symmetry)
Ixy = -2 * mass_arm * offsetX * offsetY; % kg*m^2

% Inertia matrix
I_c = [
    I_center + 4 * I_arm, Ixy, 0;
    Ixy, I_center + 4 * I_arm, 0;
    0, 0, I_center + 2 * I_arm;
];

M_rb = [total_mass*eye(3),zeros(3,3);
        zeros(3,3),I_c];



%% Generate Waypoints:

W = [ 
    0 , 25 , 40 , 50 ; 
    0 , 10 , 40 , 45; 
    5 , 80 , 100 , 110 ];


% point 0
px_0 = W(1,0+1);
py_0 = W(2,0+1);
pz_0 = W(3,0+1);

% point 1
px_1 = W(1,1+1);
py_1 = W(2,1+1);
pz_1 = W(3,1+1);

% point 2
px_2 = W(1,2+1);
py_2 = W(2,2+1);
pz_2 = W(3,2+1);

% point 3
px_3 = W(1,3+1);
py_3 = W(2,3+1);
pz_3 = W(3,3+1);

% velocities
v_0 = 100; % longitudinal path "speed"
v_0_ = [pz_1-pz_0,py_1-py_0,px_1-px_0]/norm([pz_1-pz_0,py_1-py_0,px_1-px_0]);
vx_0 = v_0_(3) * v_0;
vy_0 = v_0_(2) * v_0;
vz_0 = v_0_(1) * v_0;

v_1 = 100; % longitudinal path "speed"
v_1_ = [pz_2-pz_1,py_2-py_1,px_2-px_1]/norm([pz_2-pz_1,py_2-py_1,px_2-px_1]);
vx_1 = v_1_(3) * v_1;
vy_1 = v_1_(2) * v_1;
vz_1 = v_1_(1) * v_1;

v_2 = 100; % longitudinal path "speed"
v_2_ = [pz_3-pz_2,py_3-py_2,px_3-px_2]/norm([pz_3-pz_2,py_3-py_2,px_3-px_2]);
vx_2 = v_2_(3) * v_2;
vy_2 = v_2_(2) * v_2;
vz_2 = v_2_(1) * v_2;

v_3 = 0; % longitudinal path "speed"
% v_3_ = [pz_4-pz_3,py_4-py_3,px_4-px_3]/norm([pz_4-pz_3,py_4-py_3,px_4-px_3]);
% vx_3 = v_3_(3) * v_3;
% vy_3 = v_3_(2) * v_3;
% vz_3 = v_3_(1) * v_3;
vx_3 = v_3;
vy_3 = v_3;
vz_3 = v_3;

% computation alpha coefficient vectors
O = zeros(1,4);
Bmatrix = [ b(0)'       ,  O            ,  O            ;
            b_prime(0)' ,  O            ,  O            ;
            b(1)'       ,  O            ,  O            ;
            b_prime(1)' , -b_prime(1)'  ,  O            ;
            b_pprime(1)', -b_pprime(1)' ,  O            ;
            O           ,  b(1)'        ,  O            ;
            O           ,  b(2)'        ,  O            ;
            O           ,  b_prime(2)'  , -b_prime(2)'  ;
            O           ,  b_pprime(2)' , -b_pprime(2)' ;
            O           ,  O            ,  b(2)'        ;
            O           ,  O            ,  b(3)'        ;
            O           ,  O            ,  b_prime(3)'  ];

BCx = [ px_0 ; vx_0 ; px_1 ; 0 ; 0 ; px_1 ; px_2 ; 0 ; 0 ; px_2 ; px_3 ; vx_3 ];
alpha_x = Bmatrix\BCx;

BCy = [ py_0 ; vy_0 ; py_1 ; 0 ; 0 ; py_1 ; py_2 ; 0 ; 0 ; py_2 ; py_3 ; vy_3 ];
alpha_y = Bmatrix\BCy;

BCz = [ pz_0 ; vz_0 ; pz_1 ; 0 ; 0 ; pz_1 ; pz_2 ; 0 ; 0 ; pz_2 ; pz_3 ; vz_3 ];
alpha_z = Bmatrix\BCz;

% computing paths
sigma1 = 0:0.05:1;
sigma2 = 1.05:0.05:2;
sigma3 = 2.05:0.05:3;
sigma = [sigma1 , sigma2 , sigma3];
path_x = zeros(length(sigma),1);
path_y = zeros(length(sigma),1);
path_z = zeros(length(sigma),1);

for i = 1:length(sigma1) % for the first segment
    path_x(i) = b(sigma(i))'*alpha_x(1:4);
    path_y(i) = b(sigma(i))'*alpha_y(1:4);
    path_z(i) = b(sigma(i))'*alpha_z(1:4);
end
for i = length(sigma1)+1:length(sigma1)+length(sigma2) % for the second segment
    path_x(i) = b(sigma(i))'*alpha_x(5:8);
    path_y(i) = b(sigma(i))'*alpha_y(5:8);
    path_z(i) = b(sigma(i))'*alpha_z(5:8);
end
for i = length(sigma1)+length(sigma2)+1:length(sigma1)+length(sigma2)+length(sigma3) % for the third segment
    path_x(i) = b(sigma(i))'*alpha_x(9:12);
    path_y(i) = b(sigma(i))'*alpha_y(9:12);
    path_z(i) = b(sigma(i))'*alpha_z(9:12);
end

% combine both paths into a single path array
path = [ path_x' ; path_y' ; path_z'];

% plot path alone (positions depending on sigma)
figure, plot3(path_x,path_y,path_z,'LineWidth',2)
title('3D position')
axis equal
grid on
hold on

plot3(W(1,:),W(2,:),W(3,:),'o',...
    'LineWidth',2,...
    'MarkerSize',5,...
    'MarkerEdgeColor','r')
hold off


%-------------------------
function b_out = b(sigma)

b_out = [ 1 ; sigma ; sigma^2 ; sigma^3 ];

end

%-------------------------
function bp_out = b_prime(sigma)

bp_out = [ 0 ; 1 ; 2*sigma ; 3*sigma^2 ];

end

%-------------------------
function bp_out = b_pprime(sigma)

bp_out = [ 0 ; 0 ; 2 ; 2*3*sigma ];

end