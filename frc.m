% sepvecs = frc(theta, N, nruns)
% 
% This function generates nruns configurations of a length N polymer, 
% using the freely rotating chain model with fixed bond angle theta.
% It is returned as a single 3 x N x nruns tensor of separation vectors.
% The first dimension is the x,y,z coordinates; the second specifies which 
% monomer on the chain; and the third specifies which configuration.

function sepvecs = frc(theta, N, nruns)
    % Initialize a 3 x (N+1) x nruns matrix of zeros. We will throw away
    % the first vector later.
    sepvecs = zeros([3,N+1,nruns]);
    
    % We create two random, normalized vectors, and do not worry about the
    % bond angle between them. We throw away the first one later.
    sepvecs(:,1,:) = rand([3,1,nruns])-ones([3,1,nruns])/2;
    lng=repmat(sqrt(sum(sepvecs(:,1,:).^2)), [3,1,1]);
    sepvecs(:,1,:)=sepvecs(:,1,:)./lng;
    sepvecs(:,2,:) = rand([3,1,nruns])-ones([3,1,nruns])/2;
    lng=repmat(sqrt(sum(sepvecs(:,2,:).^2)), [3,1,1]);
    sepvecs(:,2,:)=sepvecs(:,2,:)./lng;
    
    %sepvecs(:,1,:) = repmat([0,0,1], [1,1,nruns]);
    %sepvecs(:,2,:) = repmat([0,sin(theta),cos(theta)], [1,1,nruns]);
    
    % We iterate through the polymers, generating the i-th separation
    % vector for all configurations at the same time
    for i=3:N+1
        % Quick check for the case theta==0
        if theta==0
            sepvecs(:,i,:) = sepvecs(:,i-1,:);
            continue
        elseif theta==pi
            sepvecs(:,i,:) = -sepvecs(:,i-1,:);
            continue
        end
           
        
        % We need three basis vectors: 
        % r_(i-1), r_(i-2) x r_(i-1), and r_(i-1) x (r_(i-2) x r_(i-1))
        % We call these a,b, and c
        a = sepvecs(:,i-1,:);
        x2 = sepvecs(:,i-2,:);
        b = cross(x2, a); % this is not normalized!
        bnorms = sqrt(sum(b.*b,1));
        b = b./ bnorms([1 1 1],:,:); % Normalize the b vector
        c = cross(a, b);
        % Generate a dihedral angle for each configuration, for the i-th
        % monomer
        phi = rand([1,1,nruns]) * 2 * pi;
        phi = phi([1 1 1],:,:); %  <-- repeat it for all three cordinates
        % Now the i-th separation vector is 
        % a * cos(theta) + b * sin(phi) * sin(theta) + c * cos(phi) * sin(theta)
        sepvecs(:,i,:) = a.*cos(theta) ...
            + b.*sin(phi).*sin(theta) ...
            + c.*cos(phi).*sin(theta);
        % Now we increment i, and make the separation vector for the next
        % i...
    end
    
    % Throw away the first vector, which has the wrong bond angle
    sepvecs=sepvecs(:,2:end,:);
end