clear all
% Vertical swg, patterning the width of the wg
% swg strips are parallel to propagation
w_0=0.400;
L_0=15;
h=0.220;
swg_pitch=0.1;
dc=0.5;
disp('Starting Horizental swg patterning');
v_n_period=ceil(w_0/swg_pitch);
w=v_n_period*swg_pitch;
sprintf('v_n_period=%d  swg_pitch= %2.4f',v_n_period,w)
v_z_teeth_min={};
v_z_teeth_max={};
for i=1:v_n_period
    v_z_teeth_min{i}=(i-1)*swg_pitch;
    v_z_teeth_max{i}=(i-1+dc)*swg_pitch;
    v_z_teeth{i}=[v_z_teeth_min{i} v_z_teeth_max{i}];
    %
    v_z_groove_min{i}=(i-1+dc)*swg_pitch;
    v_z_groove_max{i}=(i)*swg_pitch;
    v_z_groove{i}=[v_z_groove_min{i} v_z_groove_max{i}];
end
disp('Vertical Teeth locations');
for i=1:v_n_period
    sprintf('v_z_teeth{%d} = %2.4f   %2.4f',i,v_z_teeth{i}(1),v_z_teeth{i}(2))
end
disp('Vertical Groove locations');
for i=1:v_n_period
    sprintf('v_z_groove{%d} = %2.4f   %2.4f',i,v_z_groove{i}(1),v_z_groove{i}(2))
end
%%%%%%%%%
% Horizental swg, patterning the length of the wg
% swg strips are perpandecular to propagation
disp('Starting Horizental swg patterning');
h_n_period=ceil(L_0/swg_pitch);
L=h_n_period*swg_pitch;
sprintf('h_n_period=%d  swg_pitch= %2.4f',h_n_period,L)
h_z_teeth_min={};
h_z_teeth_max={};
%
for i=1:h_n_period
    h_x_teeth_min{i}=(i-1)*swg_pitch;
    h_x_teeth_max{i}=(i-1+dc)*swg_pitch;
    h_x_teeth{i}=[h_x_teeth_min{i} h_x_teeth_max{i}];
    %
    h_x_groove_min{i}=(i-1+dc)*swg_pitch;
    h_x_groove_max{i}=(i)*swg_pitch;
    h_x_groove{i}=[h_x_groove_min{i} h_x_groove_max{i}];
end
%
disp('Horizental Teeth locations');
for i=1:h_n_period
    sprintf('h_x_teeth{%d} = %2.4f   %2.4f',i,h_x_teeth{i}(1),h_x_teeth{i}(2))
end
disp('Horizental Groove locations');
for i=1:h_n_period
    sprintf('h_x_groove{%d} = %2.4f   %2.4f',i,h_x_groove{i}(1),h_x_groove{i}(2))
end
vfileId=fopen('v_swg_patterns.txt','w');
hfileId=fopen('h_swg_patterns.txt','w');
for i=1:v_n_period
    fprintf(vfileId,'%2.4f, %2.4f, %2.4f, %2.4f\n', v_z_teeth{i}(1),v_z_teeth{i}(2),v_z_groove{i}(1),v_z_groove{i}(2));
end
for i=1:h_n_period
    fprintf(hfileId,'%2.4f, %2.4f, %2.4f, %2.4f\n', h_x_teeth{i}(1),h_x_teeth{i}(2),h_x_groove{i}(1),h_x_groove{i}(2));
end
