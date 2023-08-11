%Function to read data from input files Taylor.in
function [eulerangles,iopt,nslip,ns1,ns2,bs,CRSS,nsteps,tdot,L,ihard,lin_hard_const,tau0_1,tau1_1,theta0_1,theta1_1,tau0_2,tau1_2,theta0_2,theta1_2] = read_data()

fid = fopen('Taylor.in');
tline = fgetl(fid);
tline2 = fgetl(fid);
eulerangles = load(tline2);     %Read the euler angles file

delimiterIn = ' ';
headerlinesIn = 3;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
iopt = input.data;              %Choice of crystal structure

headerlinesIn = headerlinesIn + 2;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
nslip = 2*input.data;           %Total number of slip systems

headerlinesIn = headerlinesIn + 2;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
n = input.data;                 %Planes
n_size = size(n); 
n_rows = n_size(1,1);
if iopt == 1                    %If FCC
    ns1 = n(1:n_rows,1:3);
    ns2 = n(1:n_rows,1:3);
elseif iopt == 2                %If BCC
    ns1 = n(1:n_rows,4:6);
    ns2 = n(1:n_rows,7:9);
else
    ns1 = n(1:n_rows,4:6);
    ns2 = n(1:n_rows,4:6);
end

headerlinesIn = headerlinesIn + n_rows + 1;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
b = input.data;                 %Directions
n_size = size(b); 
n_rows = n_size(1,1);
if iopt == 1                    %If FCC
    bs = b(1:n_rows,1:3);
else                            %If BCC
    bs = b(1:n_rows,4:6);
end

headerlinesIn = headerlinesIn + n_rows + 1;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
CRSS = input.data;              %CRSS for the slip system

headerlinesIn = headerlinesIn + 2;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
nsteps = input.data;            %Number of increments

headerlinesIn = headerlinesIn + 2;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
tdot = input.data;              %Time tdot

headerlinesIn = headerlinesIn + 2;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
L = input.data;                 %Velocity gradient tensor

headerlinesIn = headerlinesIn + 4;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
ihard = input.data;             %Choice of hardening

headerlinesIn = headerlinesIn + 2;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
lin_hard_const = input.data;    %Linear hardening constant

headerlinesIn = headerlinesIn + 2;
input = importdata('Taylor.in',delimiterIn,headerlinesIn);
voce_hard = input.data;         %Voce hardening parameters
tau0_1 = voce_hard(1,1);    
tau1_1 = voce_hard(1,2);
theta0_1 = voce_hard(1,3);
theta1_1 = voce_hard(1,4);
tau0_2 = voce_hard(2,1);    
tau1_2 = voce_hard(2,2);
theta0_2 = voce_hard(2,3);
theta1_2 = voce_hard(2,4);
end