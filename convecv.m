function H = convecv(u1,u2,d1,d2)
% Calculates the convective terms in step 1 of the mom eqn. Boundary
% values are left blank(=0). Just to control, whether the inversion
% statement for v is correct.
H = zeros(size(u1));
for i = 2:(length(H(:,1))-1) % iterating over x
    for j = 2:(length(H(1,:))-1) % iterating over y
        H(i,j) = -1/(4*d1)*(u1(i,j+1)^2 + 2*u1(i,j+1)*u1(i,j) - ...
            2*u1(i,j)*u1(i,j-1) - u1(i,j-1)^2);
        H(i,j) = H(i,j) - 1/(4*d2)*((u1(i+1,j)+u1(i,j))*(u2(i+1,j)+u2(i+1,j-1)) - ...
            (u1(i,j)+u1(i-1,j))*(u2(i,j)+u2(i,j-1)));
    end
end
end