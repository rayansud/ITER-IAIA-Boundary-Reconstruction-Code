
%% function [Br,Bz,Flux]=VectorizedGreenFunction(measurementR,measurementZ,sourceR,sourceZ)
% 
%	Green functions
%
%	INPUT:		measurementR,measurementZ	-> list of co-ordinates of the measurement coils
%				sourceR,sourceZ             -> list of co-ordinates of the current sources
%
%	OUTPUT:     Br		-> matrix of r component of the field
%			    Bz		-> matrix of z component of the field
%               Flux    -> matrix of poloidal flux function
%
%  Author:	Zabeo Luca
%  Version: 1.0
%  Data:    04/16/2012
%  Modified: Rayan Sud, 18/07/2019
% 
% Output matrix schema:
% +-----------------------+-----------------------+-----+-----------------------+
% | source1, measurement1 | source2, measurement1 | ... | sourceN, measurement1 |
% +-----------------------+-----------------------+-----+-----------------------+
% | source1, measurement2 | source2, measurement2 |  .. |           ..          |
% +-----------------------+-----------------------+-----+-----------------------+
% |           ..          |           ..          |  .. |           ..          |
% +-----------------------+-----------------------+-----+-----------------------+
% | source1, measurementN | source2, measurementN | ... | sourceN, measurementN |
% +-----------------------+-----------------------+-----+-----------------------+						

function [Br,Bz,Flux]=VectorizedGreenFunction(measurementR,measurementZ,sourceR,sourceZ)

% Copies the vectors of position, to turn them into a matrix of the appropriate size for solving for every combination
Measurement_mat_R= repmat(measurementR,[1,size(sourceR)]);
Measurement_mat_Z= repmat(measurementZ,[1,size(sourceR)]);
Source_mat_R= repmat(transpose(sourceR),[size(measurementR),1]);
Source_mat_Z= repmat(transpose(sourceZ),[size(measurementZ),1]);

% Calculating the Green function in closed form, using elliptic integrals

h=(Measurement_mat_Z-Source_mat_Z);
u = sqrt((Measurement_mat_R+Source_mat_R).^2+h.^2);
k2 = 4.0.*Measurement_mat_R.*Source_mat_R./(u.^2);
v2 = Measurement_mat_R.^2 + Source_mat_R.^2 + h.^2;
w2 = Source_mat_R.^2 - Measurement_mat_R.^2 - h.^2;
d2 = (Source_mat_R -Measurement_mat_R).^2 + h.^2;
[K,E]=ellipke(k2,eps*1000);
Br = (2.0e-07.* h./(Measurement_mat_R.*u)).*(((v2./d2).*E)-K);
Bz = (2.0e-07./u).*(((w2./d2).*E)+K);
Flux=(4.0e-7.*sqrt(Measurement_mat_R.*Source_mat_R)./sqrt(k2)).*((1-k2.*0.5).*K-E);