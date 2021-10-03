clear
clc
% Generating Cross Sectional Images of TPMS Geometries 

%Inputs:

method = 1; 
%1:Simple 
%2:CubicRadialGradient 
%3:CylindricalRadialGradient 
%4:CylindricalLognitudalGradient
%5:CAD Gradient [Save your CAD cross-sectional images at D:/A]

xrange = 1; %number of unit cells in x direction
yrange = 1; %number of unit cells in y direction
zrange = 1; %number of unit cells in z direction

VF = 0.5; %VF is volume fraction

MeshNumPerCell = 10000;
step = (VF/MeshNumPerCell)^(1/3);%domain discretization step (resolution)
step=0.005; %domain discretization step (resolution)

VF1 = 0.5; %Start volume fraction in gradient designs
VF2 = 0.5; %End volume fraction in gradient designs

% Change f function under each "mothod" below based on the model geometry of your interest.
%P-Surfac: cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)
%C = (VF-0.5)/0.2862;
%ca = (VF1-0.5)/0.2862;
%cb = (VF2-0.5)/0.2862;

%D-Surface: (cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)-sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z))
%C = (VF-0.5)/0.5845660;
%ca = (VF1-0.5)/0.5845660;
%cb = (VF2-0.5)/0.5845660;

%G-Surface: sin(2*pi*x)*cos(2*pi*y) + sin(2*pi*z)*cos(2*pi*x) + sin(2*pi*y)*cos(2*pi*z)
C = (VF-0.5)/0.32851063;
ca = (VF1-0.5)/0.32851063;
cb = (VF2-0.5)/0.32851063;

%Diamand-Surface: sin(2*pi*x) *sin(2*pi*y) *sin(2*pi*z) +sin(2*pi*x) * cos(2*pi*y) * cos(2*pi*z) +cos(2*pi*x) * sin(2*pi*y) * cos(2*pi*z) + cos(2*pi*x) * cos(2*pi*y) * sin(2*pi*z)
%C = (VF-0.5)/0.41265;
%ca = (VF1-0.5)/0.41265;
%cb = (VF2-0.5)/0.41265;

%I-WP Surface: 2*(cos(2*pi*x)*cos(2*pi*y) + cos(2*pi*y)*cos(2*pi*z) + cos(2*pi*z)*cos(2*pi*x) ) - ( cos(4*pi*x) + cos(4*pi*y) + cos(4*pi*z))
%C = (VF-0.4692)/0.1359;
%ca = (VF1-0.4692)/0.1359;
%cb = (VF2-0.4692)/0.1359;

%F-RD Surface: 4*cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z) - ( cos(4*pi*x)*cos(4*pi*y) + cos(4*pi*x)*cos(4*pi*z)+cos(4*pi*y)*cos(4*pi*z) )
%C = (VF-0.4317)/0.215475;
%ca = (VF1-0.4317)/0.215475;
%cb = (VF2-0.4317)/0.215475;



%% CODE FUNCTIONS

nx = floor(xrange/step);
ny = floor(yrange/step);
nz = floor(zrange/step);

img(nx,ny,[1 2 3])=[0 0 0];
meshnum = 0;

if method == 1 

for k=1:nz
    

    
    for j=1:ny
    for i=1:nx
        x = step*(i-0.5);
        y = step*(j-0.5);
        z = step*(k-0.5);
        
		%copy the function  f from the top list
        f = sin(2*pi*x)*cos(2*pi*y) + sin(2*pi*z)*cos(2*pi*x) + sin(2*pi*y)*cos(2*pi*z);
        f = f - C;

        if f >= 0
           img(i,j,[1 2 3])= [255 255 255];
           meshnum = meshnum +1;
        end
        

        
    end 
    end
    
  
imwrite(img,['D:\',num2str(k),'.bmp'])

clear img
img(nx,ny,[1 2 3])=[0 0 0];

end
end



if method == 2

for k=1:nz
    
    for j=1:ny
    for i=1:nx
        x = step*(i-0.5);
        y = step*(j-0.5);
        z = step*(k-0.5);
        r = sqrt(x*x+y*y+z*z);
        tet = atan(z/x);
        alf = atan(y/x);
        lx = xrange*sqrt(1+(x/z)^2+(y/z)^2);
        ly = yrange*sqrt(1+(z/y)^2+(x/y)^2);
        lz = zrange*sqrt(1+(z/x)^2+(y/x)^2);
        L = min([lx,ly,lz]);
        
        c = (cb-ca)*(r)/(L)  + ca;
		
		%copy the function  f from the top list
        f = (cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)-sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z));
        f = f-c;

        if f <= 0
           img(i,j,[1 2 3])= [255 255 255];          
           meshnum = meshnum +1;
        end

      
    end 
    end
    
imwrite(img,['D:\',num2str(k),'.bmp'])
clear img
img(nx,ny,[1 2 3])=[0 0 0];

end
    
end



if method == 3
    r0 = min(xrange,yrange)/2;
    x0 = xrange/2;
    y0 = yrange/2;
    z0 = zrange/2;
for k=1:nz
    
    for j=1:ny
    for i=1:nx
        x = step*(i-0.5)-x0;
        y = step*(j-0.5)-y0;
        z = step*(k-0.5);
        r = sqrt(x*x+y*y);
        
        if r<=r0
        c = (cb-ca)*(r)/(r0)  + ca;
        
		%copy the function  f from the top list
        f = sin(2*pi*x)*cos(2*pi*y) + sin(2*pi*z)*cos(2*pi*x) + sin(2*pi*y)*cos(2*pi*z);
        f = f-c;

        if f <= 0
           img(i,j,[1 2 3])= [255 255 255];    
           meshnum = meshnum +1;
        end
        end
      
    end 
    end
    
imwrite(img,['D:\',num2str(k),'.bmp'])
clear img
img(nx,ny,[1 2 3])=[0 0 0];

end
end



if method == 4
    r0 = min(xrange,yrange)/2;
    %r0 = max(xrange,yrange)*10;
    x0 = xrange/2;
    %x0 = 0;
    y0 = yrange/2;
    %y0 = 0;
for k=1:nz
    
    for j=1:ny
    for i=1:nx
        x = step*(i-0.5)-x0;
        y = step*(j-0.5)-y0;
        z = step*(k-0.5);
        r = sqrt(x*x+y*y);
        %r = -1;
        
        if r<=r0
        c = (cb-ca)*(z)/(zrange)  + ca;
        
		%copy the function  f from the top list
		f = (cos(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)-sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z));
        f = f-c; 
		
        if f <= 0
           img(i,j,[1 2 3])= [255 255 255];       
           meshnum = meshnum +1;
        end
        end
      
    end 
    end
    
imwrite(img,['D:\',num2str(k),'.bmp'])
clear img
img(nx,ny,[1 2 3])=[0 0 0];

end
end


if method == 5
    image0 = imread('D:\A\1.bmp');
    nx = size(image0,1);
    ny = size(image0,2);
    nz = input('Please save CAD cross-sectional Images in D:\\A and enter the number of images: ');
    img(nx,ny,[1 2 3])=[0 0 0];
    
    for k=1:nz
        
        image0 = imread(['D:\A\',num2str(k),'.bmp']);
       
        boundaries = bwboundaries(image0,'noholes');
        [B,L,N,A] = bwboundaries(image0);
        for zz = 1:N
            if (nnz(A(:,zz)) > 0)
                for l = find(A(:,zz))'
                    boundary = B{l};
                    boundaries{zz}=[boundaries{zz};boundary];
                end
            end
        end

        max_dist = zeros(1,size(boundaries,1));
        
             
        for j=1:ny
            for i=1:nx

                if image0(i,j)==1;
                    min_dist_vec = zeros(1,size(boundaries,1));
                    max_dist_vec = zeros(1,size(boundaries,1));
                    for m=1:size(boundaries,1)
                        aaa = 1:size(boundaries{m},1);
                        bbb = [i,j];
                        XX = meshgrid(bbb,aaa);
                        vec = XX-boundaries{m};
                        dist_vec = sqrt(vec(:,1).^2 + vec(:,2).^2);
                        min_dist_vec(m) = min(dist_vec);
                        max_dist_vec(m) = max(dist_vec);                     
                    end         
                    
                    if max_dist( min_dist_vec==min(min_dist_vec) ) < min_dist_vec(min_dist_vec==min(min_dist_vec))
                        max_dist( min_dist_vec==min(min_dist_vec) ) = min_dist_vec(min_dist_vec==min(min_dist_vec));
                    end
                end

            end
        end
        
       
        
        for j=1:ny
            for i=1:nx
                if image0(i,j)== 1;
				
                    min_dist_vec = zeros(1,size(boundaries,1));
                    max_dist_vec = zeros(1,size(boundaries,1));
                    for m=1:size(boundaries,1)
                        aaa = 1:size(boundaries{m},1);
                        bbb = [i,j];
                        XX = meshgrid(bbb,aaa);
                        vec = XX-boundaries{m};
                        dist_vec = sqrt(vec(:,1).^2 + vec(:,2).^2);
                        min_dist_vec(m) = min(dist_vec);
                        max_dist_vec(m) = max(dist_vec);
                    end                
                    distance0 = min(min_dist_vec);
                
                
                x = step*(i-0.5);
                y = step*(j-0.5);
                z = step*(k-0.5);
                
				%copy the function  f from the top list
                f = cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z);
                
	
                C = ca + (cb-ca)/(max_dist(min_dist_vec==min(min_dist_vec)))*distance0;
                f = f-C;
                
                
                if f <= 0
                    img(i,j,[1 2 3])= [255 255 255];
                    meshnum = meshnum +1;
                end
                
                end
            end
        end
        imwrite(img,['D:\',num2str(k),'.bmp']);
        clear img
        img(nx,ny,[1 2 3])=[0 0 0];
        clc
        disp(['Modeling Image Number: ', num2str(k), '/', num2str(nz), ' Meshes Generater: ', num2str(meshnum)]);

    end
end


clc
disp([num2str(meshnum), ' Hex Meshes will be generated with step size: ', num2str(step)]);


 
