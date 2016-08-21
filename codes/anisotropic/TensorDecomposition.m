function [b,weights]=TensorDecomposition(tensors)

% Copyright (c) 2015 Jean-Marie Mirebeau

% Check data, reshape tensors
s=size(tensors);
assert(length(s)==2)
assert(s(2)==3)

% Initialize superbases
s(2)=2;
b={ones(s),ones(s),-ones(s)};
b{1}(:,2)=0;
b{2}(:,1)=0;

% Selling's algorithm
posScal=0; %number of consecutive scalar products seen
counter=0;
maxIter=50; %increase for highly anisotropic tensor fields
while counter<maxIter && posScal<3
    scal = ScalarProduct(b{1},b{2},tensors);
    scal = (scal>0);
    if max(scal(:))
      b{3}(scal,:) = b{2}(scal,:)-b{1}(scal,:);
      b{2}(scal,:) = -b{2}(scal,:);
      posScal=0;
    else
      posScal = posScal+1;
    end
    counter=counter+1;
    b={b{2},b{3},b{1}};
end

if posScal<3
    'Warning : Selling-s algorithm unterminated'
end

weights={
   -ScalarProduct(b{2},b{3},tensors)/2,
   -ScalarProduct(b{3},b{1},tensors)/2,
   -ScalarProduct(b{1},b{2},tensors)/2
};

for i=1:3 %Compute orthogonal vectors
    c=b{i};
    b{i}(:,1)=-c(:,2);
    b{i}(:,2)= c(:,1);
end

end

function scal = ScalarProduct(u,v,D)
      scal = D(:,1).*u(:,1).*v(:,1) +  D(:,2).*(u(:,1).*v(:,2) + u(:,2).*v(:,1) ) + D(:,3).*u(:,2).*v(:,2);
end
