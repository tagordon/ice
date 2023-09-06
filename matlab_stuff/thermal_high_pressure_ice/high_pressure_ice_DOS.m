
clear all;
close all;

%collect vibrational density-of-states data and compute the corresponding
%internal energy Evib

%read in constants and system-specific information
my_constants;

%read in Excel sheet specifying PVT conditions where vibrational
%density-of-states (DOS) data is available
filename = 'Conditions.xlsx';
% probably new in new matlab version -- change something... 
data = string(table2cell(readtable(filename)));

%number of different conditions (each on one row in the Excel sheet)
num_rows = size(data,1);

%starting from the row specified by start_index, read through all entries
start_index = 6;
solid_index = 1;
for i = start_index:num_rows
    
    %ignore this entry if it corresponds to a fluid state or is an empty line
    if (data(i,8) == "Fluid" || data(i,1) == "")
  
    %if not empty or fluid, then it must be a solid
    else
        dos.T_array(solid_index) = str2num(data(i,3));
        dos.a_array(solid_index) = str2num(data(i,1));
        dos.rho_array(solid_index) = str2num(data(i,6));
        dos.V_array(solid_index) = mw/dos.rho_array(solid_index);
        solid_index = solid_index + 1;
    end
end
%number of solid PVT points where the DOS is available
num_solids = solid_index - 1;   

%create an array representing the volumes normalized by Vo
dos.VovVo_array = dos.V_array/Vo;

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

    %compute the internal energy due to ionic motion (Evib) in units of J/kg
    Evib = 0;
    for i = 1:num_atom_types
        Evib = Evib + 3*num_atom(i)*calculate_Evib(sp(i),...
            dos_matrix(1,1),dos_matrix(end,1),dos.T_array(j));
    end
    dos.Evib_array(j) = Evib;
end

%plot to make sure that the spline fit gives a good representation
% figure(3);
% fnplt(dos.Evib_sp); 
% hold on;
% plot3(x(1,:),x(2,:),dos.Evib_array,'wo','markerfacecolor','k');

