
clear all;
close all;

%read in constants and system-specific information
my_constants;

filename = 'Conditions.xlsx';
data = string(table2cell(readtable(filename)));

num_rows = size(data,1);

start_index = 6;
solid_index = 1;
for i = start_index:num_rows
    
    %ignore this entry if it corresponds to a fluid state or is an empty line
    if (data(i,8) == "Fluid" || data(i,1) == "")
  
    else
        dos.T_array(solid_index) = str2num(data(i,3));
        dos.a_array(solid_index) = str2num(data(i,1));
        dos.rho_array(solid_index) = str2num(data(i,6));
        dos.V_array(solid_index) = mw/dos.rho_array(solid_index);
        solid_index = solid_index + 1;
    end
end
num_solids = solid_index - 1;

%location of current directory
this_directory = fullfile(pwd);

for j = 1:num_solids
    
    %location of the files
    filename = strcat(this_directory,'/DOS/',num2str(dos.T_array(j)),'K/S_',...
        num2str(dos.a_array(j)*10),'.dat');
    
    %read in DOS
    dos_data = string(table2cell(readtable(filename)));
    dos_matrix = str2double(dos_data(2:end,:));
    
    %convert frequencies from Hz to cm^-1
    dos_matrix(:,1) = dos_matrix(:,1)/c;
    
    %convert intensities from s to cm for each atom type
    for i = 2:size(dos_matrix,2)
        dos_matrix(:,i) = dos_matrix(:,i)*c;
    end
    
    %fit a simple smooth b-spline for each atom type
    smoothing = 0;
    for i = 1:num_atom_types
        sp(i) = spaps(dos_matrix(:,1),dos_matrix(:,i+1),smoothing);
    end

    %check that the splines adequately represent the data for each atom type 
    %(change smoothing if necessary)
%     figure(1)
%     for i = 1:num_atom_types
%         fnplt(sp(i),color(i));
%         hold on;
%         plot(dos_matrix(:,1),dos_matrix(:,i+1),'o','color',color(i));
%     end
%     hold off;
    
    %check that the integrals are normalized for each atom type
%     figure(2)
%     for i = 1:num_atom_types
%         spI(i) = fnder(sp(i),-1); 
%         fnplt(spI(i),color(i));
%         hold on;
%     end
%     hold off;
%     title(['DOS integral for each atom type']);
    
    %compute the internal energy due to ionic motion (Evib)
    Evib = 0;
    for i = 1:num_atom_types
        Evib = Evib + 3*num_atom(i)*calculate_Evib(sp(i),...
            dos_matrix(1,1),dos_matrix(end,1),dos.T_array(j));
    end
    dos.Evib_array(j) = Evib;
end

Vo = 12.7218;   %cm^3/mol
dos.VovVo = dos.V_array/Vo;





% figure(3);
% fnplt(dos.Evib_sp); 
% hold on;
% plot3(x(1,:),x(2,:),dos.Evib_array,'wo','markerfacecolor','k');

