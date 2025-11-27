function dydt = mol_rhs_non_conserve(t,y,xi,dxi,Lfun,dLdt,Du,Dv,Dw,De,Ds,a,b,c,e,s,scale_factor)
    % Unpack
    N = numel(xi);
    Y = reshape(y,[N,5]);  % columns: U,V,W,Y,Z
    U1 = Y(:,1); V1 = Y(:,2); W1 = Y(:,3); Y1 = Y(:,4); Z1 = Y(:,5);

    % Growth scaling
    L   = Lfun(t);
    g   = dLdt(t)/L;           % stretching rate
    invL2 = 1/(L*L);

    % Derivatives (Neumann BC via mirrored ghost points)
    lap = @(u) lap_neumann(u, dxi);
    d1  = @(u) dxi_neumann(u, dxi);  % first derivative u_xi

    % Diffusion terms (scaled by 1/L(t)^2)
    Uxx = lap(U1); Vxx = lap(V1); Wxx = lap(W1); Yxx = lap(Y1); Zxx = lap(Z1);

    % Advection (stretching) term  - (dL/L) * xi * u_xi
    xiad = xi-0.5; % xi positions
    Ux  = d1(U1); Vx = d1(V1); Wx = d1(W1); Yx = d1(Y1); Zx = d1(Z1);

%     advU = - g .* xiad .* Ux .* 0;
%     advV = - g .* xiad .* Vx .* 0;
%     advW = - g .* xiad .* Wx .* 0;
%     advY = - g .* xiad .* Yx .* 0;
%     advZ = - g .* xiad .* Zx .* 0;

    advU = - g .* xiad .* Ux;
    advV = - g .* xiad .* Vx;
    advW = - g .* xiad .* Wx;
    advY = - g .* xiad .* Yx;
    advZ = - g .* xiad .* Zx;

    % dilution
    diluU = - g .* U1 .* 0;
    diluV = - g .* V1 .* 0;
    diluW = - g .* W1 .* 0;
    diluY = - g .* Y1 .* 0;
    diluZ = - g .* Z1 .* 0;

    % Reactions
    R1 = (a - U1 + (U1.^2)./V1) .* scale_factor;
    R2 = (b - V1 + U1.^2) .* scale_factor;
    R3 = (c - W1 + Y1 ./ (0.1 + 0.9*1.5*Z1)) .* scale_factor;
    R4 = (e - Y1 + U1) .* scale_factor;
    R5 = (s - Z1 + U1)*0.001 .* scale_factor;

    % Time derivatives
    dU = Du*invL2*Uxx .* scale_factor + advU + R1 + diluU;
    dV = Dv*invL2*Vxx .* scale_factor + advV + R2 + diluV;
    dW = Dw*invL2*Wxx .* scale_factor + advW + R3 + diluW;
    dY = De*invL2*Yxx .* scale_factor + advY + R4 + diluY;
    dZ = Ds*invL2*Zxx .* scale_factor + advZ + R5 + diluZ;

    dydt = [dU; dV; dW; dY; dZ];
end

function lap = lap_neumann(u,dx)
    % Second derivative with zero-flux (Neumann) BC via mirrored ghosts
    N = numel(u);
    lap = zeros(N,1);
    % interior
    lap(2:N-1) = (u(3:N) - 2*u(2:N-1) + u(1:N-2)) / (dx*dx);
    % boundaries (mirror neighbors)
    lap(1)   = 2*(u(2)   - u(1))   / (dx*dx);
    lap(N)   = 2*(u(N-1) - u(N))   / (dx*dx);
end

function dudx = dxi_neumann(u,dx)
    % First derivative with zero-flux BC via mirrored ghosts
    N = numel(u);
    dudx = zeros(N,1);
    % interior (central)
    dudx(2:N-1) = (u(3:N) - u(1:N-2)) / (2*dx);
    % boundaries: symmetric (zero normal derivative) -> use mirrored values
    % which makes first derivative zero at boundaries
    dudx(1) = 0;
    dudx(N) = 0;
end