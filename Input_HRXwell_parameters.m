clear all
% ------import the coordinates of nodes in the background mesh
lx=input('x length of the pit=');%width of the pit
ly=input('y length of the pit='); %length of the pit
lz=input('z length of the pit='); %depth of the pit

fprintf('Choose position of the well\n')
xwc=input('well center x='); %x coordinate of the well center
ywc=input('well center y='); %y coordinate of the well center
depth=input('well center depth='); %depth of the well center
zwc=lz/2-depth; %z-coordinate of the well center

%fprintf('length of the horizontal part of the well')

diam=input('diameter of the well='); %diameter of the well
r=diam/2; %radius of the well 

i=0;
%left hand side boundary face
node(i+1,:)=[-lx/2,-ly/2,-lz/2];
node(i+2,:)=[lx/2,-ly/2,-lz/2];
node(i+3,:)=[lx/2,-ly/2,lz/2];
node(i+4,:)=[-lx/2,-ly/2,lz/2];

%right hand side bounary face
node(i+5,:)=[-lx/2,ly/2,-lz/2];
node(i+6,:)=[lx/2,ly/2,-lz/2];
node(i+7,:)=[lx/2,ly/2,lz/2];
node(i+8,:)=[-lx/2,ly/2,lz/2];

lhw=input('length of the horizontal part of the well casing(including media)='); 
node(i+9,:)=[-lx/2,-lhw/2,-lz/2];
node(i+10,:)=[lx/2,-lhw/2,-lz/2];
node(i+11,:)=[lx/2,-lhw/2,lz/2];
node(i+12,:)=[-lx/2,-lhw/2,lz/2];

node(i+13,:)=[-lx/2,lhw/2,-lz/2];
node(i+14,:)=[lx/2,lhw/2,-lz/2];
node(i+15,:)=[lx/2,lhw/2,lz/2];
node(i+16,:)=[-lx/2,lhw/2,lz/2];

i=16;

fbd1_mark=1; %dirichelet boudary condition face marker
fbd2_mark=2; %Neumann boundary condition face marker
screenbd_mark=3;
media_mark=4;

face{1,1}=[1 2 3 4];
face{2,1}=[5 6 7 8];
face{3,1}=[1 2 10 9];
face{4,1}=[4 3 11 12];
face{5,1}=[3 2 10 11];
face{6,1}=[4 1 9 12];
face{7,1}=[9 10 14 13];
face{8,1}=[12 11 15 16];
face{9,1}=[11 10 14 15];
face{10,1}=[12 9 13 16];
face{11,1}=[13 14 6 5];
face{12,1}=[16 15 7 8];
face{13,1}=[15 14 6 7];
face{14,1}=[16 13 5 8];

fbd_marker(1:2,1)=fbd1_mark;
fbd_marker(3:14,1)=fbd2_mark;

j=14;

lhs=input('length of the left horizontal screen='); 
rhs=input('length of the right horizontal screen='); 

ylhw=-lhw/2-lhs; %y coordinate of the left end of the horizontal part
yrhw=lhw/2+rhs; %y coordinate of the right end of the horizontal part

nr=20;%number of subdivision of circle
dradian=(2*pi)/nr;
radian=0:dradian:(2*pi);
rad=radian(1:nr)';
cx=r*cos(rad);
cz=r*sin(rad)+zwc;

nt=5; %number subdivision of torus with some drilling degree
rt=diam; %radius of the torus

%coordinate of center of the left torus
xlt_center=xwc; 
ylt_center=ywc+ylhw;
zlt_center=zwc+rt;

degree=input('Drilling degree='); %equal to torus degree 7*pi/180
lu=(3*pi/2-degree):(degree/nt):(3*pi/2);%degree division of left torus
%%
%llw=input('length of the left hand steel casing=');%9;%length of the left hand part of the well
%z coordinate of the left torus top
zltop=max((rt+r*cos(rad))*sin(lu(1))+zlt_center);
maxlscreen=(lz/2-zltop)/sin(degree) %length of whole left cylinder screen riser
lscreen=input('length of the left angled screen='); %should be less than maxlscreen
node=[node;
r*sin(rad)+xlt_center (rt+r*cos(rad))*cos(lu(1))+ylt_center-lscreen*cos(degree) (rt+r*cos(rad))*sin(lu(1))+zlt_center+lscreen*sin(degree)];
face{j+1,1}=(i+1):(i+nr);%face 15, node17-36 circle cross section of left end face of left screen riser
fbd_marker(j+1)=fbd2_mark; 

i=i+nr;%node index i=36
j=j+1;% face index j=15
nfls=j; %j=15

for jj=1:nr-1 %face 16-34
      face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj];
      j=j+1;
end
face{j+1,1}=[i i+1-nr i+1 i+nr]; %j=34 face 35
j=j+1;  %j=35


%%
%left angled screen riser
i=i+nr; %i=56;

%% left torus screen 
for ii=1:nt+1
      node=[node;r*sin(rad)+xlt_center (rt+r*cos(rad))*cos(lu(ii))+ylt_center (rt+r*cos(rad))*sin(lu(ii))+zlt_center];   
end
for ii=1:nt
    for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj]; %face 36-54, 56-74, 76-94, 96-114, 116-134
        j=j+1;
    end
    face{j+1,1}=[i i+1-nr i+1 i+nr];%face 55, 75, 95, 115, 135
        j=j+1; %j=55, 75, 95, 115, 135
        i=i+nr; %i=76, 96, 116, 136, 156    
end

%% left horizontal screen
node=[node;node((i-nr+1):i,1) node((i-nr+1):i,2)+lhs node((i-nr+1):i,3)];%node:157-176
for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj];%face 136-154
        j=j+1;
end
face{j+1,1}=[i i+1-nr i+1 i+nr];%face 155
j=j+1; %
i=i+nr;%i=176
fbd_marker(nfls+1:j)=screenbd_mark;

%% left media part
face{j+1,1}=(i-nr+1):i; %face 156, left boundary face of the left media
j=j+1; %j=156
fbd_marker(j)=media_mark;
lml=input('length of left media=');
rml=input('length of right media=');
node=[node;node((i-nr+1):i,1) node((i-nr+1):i,2)+lml node((i-nr+1):i,3)];
for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj];%face 157-175
        j=j+1;
end
face{j+1,1}=[i i+1-nr i+1 i+nr];%face 176
j=j+1;%j=176
i=i+nr;%i=196
fbd_marker(j-nr+1:j)=fbd2_mark;%face 157-176
face{j+1,1}=(i-nr+1):i; %face 177, right boundary face of the left media
j=j+1;%j=177
fbd_marker(j)=media_mark;
%% hollow part
lcw=lhw-lml-rml;%length of hollow well
node=[node;node(i-nr+1:i,1) node(i-nr+1:i,2)+lcw node(i-nr+1:i,3)];%node 197-216
for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj];%face 178-196
        j=j+1;
end
face{j+1,1}=[i i+1-nr i+1 i+nr];%face 197
j=j+1; %j=197
i=i+nr; %i=216    
fbd_marker(j-nr+1:j)=fbd2_mark;
%% right media part
face{j+1,1}=(i-nr+1):i; %face 198
j=j+1; %j=198
fbd_marker(j)=media_mark; %left boundary face of the right media
node=[node;node(i-nr+1:i,1) node(i-nr+1:i,2)+rml node(i-nr+1:i,3)];%node 217-236
for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj];%face 199-217
        j=j+1;
end
face{j+1,1}=[i i+1-nr i+1 i+nr];%face 218
j=j+1;%j=218
i=i+nr; %i=236
fbd_marker(j-nr+1:j)=fbd2_mark;%face 199-218
face{j+1,1}=(i-nr+1):i; %face 219
j=j+1;
fbd_marker(j)=media_mark; %face 219, right boundary face of the right media
nfrm=j; 

%% right horizontal screen
node=[node;node(i-nr+1:i,1) node(i-nr+1:i,2)+rhs node(i-nr+1:i,3)];%node 237-256
for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj];%face 220-238
        j=j+1;
end
face{j+1,1}=[i i+1-nr i+1 i+nr];%face 239
j=j+1;%j=239
i=i+nr;%i=256

%% right torus screen
xrt_center=xwc;
yrt_center=ywc+yrhw;
zrt_center=zwc+rt;
ru=(3*pi/2):(degree/nt):(3*pi/2+degree); 
for ii=2:nt+1 %node 257-356
      node=[node;r*sin(rad)+xrt_center (rt+r*cos(rad))*cos(ru(ii))+yrt_center (rt+r*cos(rad))*sin(ru(ii))+zrt_center];   
end
for ii=1:nt
    for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj];%face 240-258,...,320-338
        j=j+1;
    end
    face{j+1,1}=[i i+1-nr i+1 i+nr];%face 259,...,339
        j=j+1; %j=259,...,339
        i=i+nr; %i=276,...,356    
end

%%  right angled screen
rscreen=input('length of the right angled screen='); 
%node 357-376
node=[node;
r*sin(rad)+xrt_center (rt+r*cos(rad))*cos(ru(end))+yrt_center+rscreen*cos(degree) (rt+r*cos(rad))*sin(ru(end))+zrt_center+rscreen*sin(degree)];
for jj=1:nr-1
        face{j+1,1}=[i+jj-nr i+jj+1-nr i+jj+1 i+jj]; %face 340-358
        j=j+1;
end
face{j+1,1}=[i i+1-nr i+1 i+nr];%face 359
j=j+1;%j=359
i=i+nr;%i=376
nfrs=j;
fbd_marker(nfrm+1:nfrs)=screenbd_mark;%face 220-359

face{j+1,1}=(i-nr+1):i; %face 360
j=j+1;%j=360;
fbd_marker(j,1)=fbd2_mark;
nnode1=length(node);  
nf1=length(face);

%% All the nodes and faces of the casing
i=16+nr*((nt+1)+1);%the number of all the points before the left end of the horizontal well casing
thick=0.1;
rout=r+thick;
ylc=-lhw/2;
node=[node;rout*sin(rad)+xlt_center ylc*ones(nr,1) (rt+rout*cos(rad))*sin(lu(end))+zlt_center];   
nn=length(node);
node=[node;node(nn-nr+1:nn,1) node(nn-nr+1:nn,2)+lhw node(nn-nr+1:nn,3)];
nnode=length(node);

lefthole=(nnode-2*nr+1):(nnode-nr);
righthole=(nnode-nr+1):(nnode);

    for jj=1:nr-1
        face{j+1,1}=[nnode1+jj nnode1+jj+1 i+jj+1 i+jj];%face 361-379
        j=j+1;
    end
    face{j+1,1}=[nnode1+nr nnode1+1 i+1 i+nr];%face 380
        j=j+1;     
        
    for jj=1:nr-1
        face{j+1,1}=[nnode1+jj nnode1+jj+1 nnode1+jj+1+nr nnode1+jj+nr];
        j=j+1;
    end
    face{j+1,1}=[nnode1+nr nnode1+1 nnode1+1+nr nnode1+nr+nr];
        j=j+1;  
i=16+nr*(nt+6);%the number of all the points before the right end of the horizontal well 
 for jj=1:nr-1
        face{j+1,1}=[nnode1+nr+jj nnode1+nr+jj+1 i+jj+1 i+jj];
        j=j+1;
 end
 face{j+1,1}=[nnode1+nr+nr nnode1+nr+1 i+1 i+nr];
        j=j+1;       
fbd_marker(nf1+1:j)=fbd2_mark;
nf=length(face);
%% 
face{nf+1,1}=[11 10 9 12];
face{nf+2,1}=[15 14 13 16];
fbd_marker(nf+1:nf+2,1)=0;
nf=length(face);
%% generate poly file
dlmwrite('input_project_noriser.poly',[nnode 3 0 1],'-append','delimiter',' ','newline','pc'); 
wellstart=7*nr+16+1; wellend=11*nr+16;
for row = 1:(wellstart-1)
    dlmwrite('input_project_noriser.poly',[row node(row,:) 0],'-append','delimiter',' ','newline','pc'); 
    %dlmwrite('input2.txt',face{row,:},'-append','delimiter','\t','newline'
    %,'pc'); 
end
for row = wellstart:wellend
    dlmwrite('input_project_noriser.poly',[row node(row,:) 1],'-append','delimiter',' ','newline','pc'); 
    %dlmwrite('input2.txt',face{row,:},'-append','delimiter','\t','newline'
    %,'pc'); 
end
for row = (wellend+1):nnode
    dlmwrite('input_project_noriser.poly',[row node(row,:) 0],'-append','delimiter',' ','newline','pc'); 
    %dlmwrite('input2.txt',face{row,:},'-append','delimiter','\t','newline'
    %,'pc'); 
end

dlmwrite('input_project_noriser.poly',[nf 1],'-append','delimiter',' ','newline','pc'); 
%left face tangent to well casing left end
dlmwrite('input_project_noriser.poly',[2 1 fbd_marker(nf-1)],'-append','delimiter',' ','newline','pc');
  dlmwrite('input_project_noriser.poly',[length(face{nf-1,:}) face{nf-1,:}],'-append','delimiter',' ','newline','pc'); 
  dlmwrite('input_project_noriser.poly',[nr lefthole],'-append','delimiter',' ','newline','pc'); 
  dlmwrite('input_project_noriser.poly',[1 xwc ywc-lhw/2 zwc],'-append','delimiter',' ','newline','pc'); 
 
  %right face tangent to well casing right end
  dlmwrite('input_project_noriser.poly',[2 1 fbd_marker(nf)],'-append','delimiter',' ','newline','pc');
  dlmwrite('input_project_noriser.poly',[length(face{nf,:}) face{nf,:}],'-append','delimiter',' ','newline','pc'); 
  dlmwrite('input_project_noriser.poly',[nr righthole],'-append','delimiter',' ','newline','pc'); 
  dlmwrite('input_project_noriser.poly',[1 xwc ywc+lhw/2 zwc],'-append','delimiter',' ','newline','pc'); 

for row = 1:nf-2
    dlmwrite('input_project_noriser.poly',[1 0 fbd_marker(row)],'-append','delimiter',' ','newline','pc');
    dlmwrite('input_project_noriser.poly',[length(face{row,:}) face{row,:}],'-append','delimiter',' ','newline','pc'); 
end
  
dlmwrite('input_project_noriser.poly',1,'-append','delimiter',' ','newline','pc');% 3 hole
dlmwrite('input_project_noriser.poly',[1 xwc ywc zwc+r+thick/4],'-append','delimiter',' ','newline','pc');%

nregion=8; 
region_sand=1;
region_well=2;
region_screen=3;
region_media=4;

dlmwrite('input_project_noriser.poly',nregion,'-append','delimiter',' ','newline','pc');
node_aqu=[xwc+r+thick+1 ywc zwc];
dlmwrite('input_project_noriser.poly',[1 node_aqu region_sand 0.1],'-append','delimiter',' ','newline','pc');
node_aqu=[xwc+r+thick+1  ywc-lhw/2-1 zwc];
dlmwrite('input_project_noriser.poly',[2 node_aqu region_sand 0.1],'-append','delimiter',' ','newline','pc');
node_aqu=[xwc+r+thick+1 ywc+lhw/2+1 zwc];
dlmwrite('input_project_noriser.poly',[3 node_aqu region_sand 0.1],'-append','delimiter',' ','newline','pc');
node_well=[xwc ywc zwc];
dlmwrite('input_project_noriser.poly',[4 node_well region_well 0.1],'-append','delimiter',' ','newline','pc');
node_screen=[xwc ywc-lhw/2-lhs/2 zwc];
dlmwrite('input_project_noriser.poly',[5 node_screen region_screen 0.1],'-append','delimiter',' ','newline','pc');
node_screen=[xwc ywc+lhw/2+rhs/2 zwc];
dlmwrite('input_project_noriser.poly',[6 node_screen region_screen 0.1],'-append','delimiter',' ','newline','pc');
node_media=[xwc ywc-lcw/2-lml/2 zwc];
dlmwrite('input_project_noriser.poly',[7 node_media region_media 0.1],'-append','delimiter',' ','newline','pc');
node_media=[xwc ywc+lcw/2+rml/2 zwc];
dlmwrite('input_project_noriser.poly',[8 node_media region_media 0.1],'-append','delimiter',' ','newline','pc');
%%