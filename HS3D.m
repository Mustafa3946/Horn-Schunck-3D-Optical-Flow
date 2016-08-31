function [ ux, uy, uz ] = HS3D( image1, image2, alpha, iterations, ...
    uxInitial, uyInitial, uzInitial)
%This function estimates deformations between two subsequent 3-D images
%using Horn-Schunck optical flow method. 
%
%   Description :  
%
%   -image1, image2 :   two subsequent images or frames
%   -alpha :    smoothness parameter, default value is 1.
%   -iterations :   number of iterations, default value is 10.
%   -uxInitial, uyInitial, uzInitial : initial flow vectors, default value
%                                     is 0;
%
%   Reference :
%   B. K. Horn and B. G. Schunck, Determining optical ow, Cambridge, MA,
%   USA, Tech. Rep., 1980.
%
%   Author : Mohammad Mustafa
%   By courtesy of The University of Nottingham and Mirada Medical Limited,
%   Oxford, UK
%
%   Published under a Creative Commons Attribution-Non-Commercial-Share Alike
%   3.0 Unported Licence http://creativecommons.org/licenses/by-nc-sa/3.0/
%   
%   June 2012


ux=zeros(size(image1)); uy=ux; uz=ux;

if nargin==2
    alpha=1; iterations=10; 
    ux=zeros(size(image1)); uy=ux; uz=ux;
elseif nargin==3
    iterations=10;
    ux=zeros(size(image1)); uy=ux; uz=ux;
elseif nargin==6
    ux=uxInitial; uy=uyInitial; uz=uzInitial;    
end


[Ix,Iy,Iz,It]=imageDerivatives3D(image1,image2);

mask = zeros(3,3,3);
mask(:,:,1) = [0 0 0;0 1 0;0 0 0]/6;
mask(:,:,2) = [0 1 0;1 0 1;0 1 0]/6;
mask(:,:,3) = mask(:,:,1);

for i=1:iterations
    uxAvg=convn(ux,mask,'same');
    uyAvg=convn(uy,mask,'same');
    uzAvg=convn(uz,mask,'same');
    ux=uxAvg - ( Ix.*( (Ix.*uxAvg) + (Iy.*uyAvg) + (Iz.*uzAvg) + It))...
        ./ ( alpha.^2 + Ix.^2 + Iy.^ 2 + Iz.^ 2);
    uy=uyAvg - ( Iy.*( (Ix.*uxAvg) + (Iy.*uyAvg) + (Iz.*uzAvg) + It))...
        ./ ( alpha.^2 + Ix.^2 + Iy.^ 2 + Iz.^ 2);
    uz=uzAvg - ( Iz.*( (Ix.*uxAvg) + (Iy.*uyAvg) + (Iz.*uzAvg) + It))...
        ./ ( alpha.^2 + Ix.^2 + Iy.^ 2 + Iz.^ 2);
end

ux(isnan(ux))=0;
uy(isnan(uy))=0;
uz(isnan(uz))=0;
end

