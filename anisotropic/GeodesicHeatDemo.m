% Copyright (c) 2015 Jean-Marie Mirebeau

% Important note : Matlab uses an image-as-matrix convention,
% that is not much compatible with anisotropic PDEs.
% Hence this code contains a lot of transpositions.

n=99; %n must be odd, or zero divide.
[x,y]=meshgrid(-1:2/n:1,-1:2/n:1);
x=x'; y=y'; % !! transpose !!
s = size(x);

% -------- Generate a field of diffusion tensors --------
eVal1 = 0.01*ones(s); eVal2 = 1*ones(s); %eigenvalues

eVec1 = zeros([s,2]); eVec2=eVec1; %eigenvectors
r=sqrt(x.*x+y.*y);
eVec1(:,:,1)= x./r; eVec1(:,:,2)=y./r;
eVec2(:,:,1)=-y./r; eVec2(:,:,2)=x./r;

%diffusion tensors with the above eigenvectors and eigenvalues.
%Here : diffusion in a circular fashion, around the image center.
tensors = ones([s,3]);
tensors(:,:,1)=eVal1.*eVec1(:,:,1).*eVec1(:,:,1)+eVal2.*eVec2(:,:,1).*eVec2(:,:,1);
tensors(:,:,2)=eVal1.*eVec1(:,:,1).*eVec1(:,:,2)+eVal2.*eVec2(:,:,1).*eVec2(:,:,2);
tensors(:,:,3)=eVal1.*eVec1(:,:,2).*eVec1(:,:,2)+eVal2.*eVec2(:,:,2).*eVec2(:,:,2);

%tensors(:,:,1)=1; tensors(:,:,2)=0; tensors(:,:,3)=1; %Uniform diffusion
%tensors(:,:,1)=1; tensors(:,:,2)=0; tensors(:,:,3)=0.01; %x-axis diffusion


% ---------- Generate the operator matrix ----------
A=DiffusionSparseMatrix(tensors);
maxTimeStep = 1./full(max(A(:)));
dt=0.5*maxTimeStep;

% ----------- Demo : anisotropic diffusion, circularly. ------
% image = random('Uniform',0,1,s(2),s(1)); % diffused image, matlab convention
I = rand(s(1),s(2));

%Matrix vector product requires reshape
I=reshape(I',[prod(s),1]); % !! transpose !!
clf;
for i=1:100
    I=I-dt*A*I;
    imagesc( reshape(I,s)' ); drawnow;
end
I=reshape(I,s)'; % !! transpose !!

imagesc(I)
% pause()

% ---------- Demo : distance map from heat. --------
seed=[10,40]; % seed for distance computations. (xCoord, yCoord), in pixels.
I=0.*I;
I(seed(2),seed(1))=1;

%Matrix vector product requires reshape
I=reshape(I',[prod(s),1]); % !! transpose !!
I = (0.01*speye(prod(s)) + dt*A)\I;
I=reshape(I,s)'; % !! transpose !!

I=-log(I);
I=I-min(I(:));
I=I/max(I(:));

I(seed(2),seed(1))=1;
% seed is a white dot.
% black : close to seed, for the riemannian distance.
% white : far from seed, for the riemannian distance.
imagesc(I);
