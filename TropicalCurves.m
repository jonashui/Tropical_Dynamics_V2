function TL = TropicalCurves(F,G,dF,dG)

%%% Error if system is too small
if size(F,1) + size(G,1) < 3
    error('There are less than three terms.')
end

%%% Initialize variables
T = [F;G];
dT = [dF,dG];
cT = T(:,1);
GT = T(:,2:3);
NT = [-GT, ones(length(dT),1)];
nF = size(F,1);
nG = size(G,1);
nT = size(T,1);
nP = 0;
nV = 0;
nE = 0;
tol = 1e-12;

%%% Line classification constants
BENT_CROSSING = 1;
STRAIGHT_CROSSING = 2;
SLIDING = 3;
NULLCLINE = 4;
CROSS_N_SLIDE = 5;

%%% Output structure
TL.pointIndices = zeros(3,nT);
TL.vectorIndices = zeros(2,nT);
TL.tropicalPointIndices = zeros(4,nT);
TL.vectors = zeros(2,nT);
TL.points = zeros(2,nT);
TL.numPoints = zeros(nT,nT);
TL.pointMap = zeros(nT,nT,2);
TL.vectorMap = zeros(nT,nT);
TL.lineClassification = zeros(nT,nT);
TL.nullClines = zeros(nT,nT,nT);
TL.tropicalPoints = zeros(2,nT);
TL.tropicalPointMap = zeros(nT,nT,nT);
TL.numTropicalPoints = zeros(nT,nT);
TL.distanceToTropicalPoint = zeros(nT,nT,nT);
TL.isEquilibrium = zeros(1,nT);

%%% Warn if system is generic
for i = 1:nT-1
    for j = i+1:nT
        if isequal(T(i,:), T(j,:))
            if j <= nF
                iT = i;
                jT = j;
                Ti = 'F';
                Tj = 'F';
            elseif i > nF
                iT = i - nF;
                jT = j - nF;
                Ti = 'G';
                Tj = 'G';
            elseif i <= nF
                iT = i;
                jT = j - nF;
                Ti = 'F';
                Tj = 'G';
            else
                iT = i - nF;
                jT = j;
                Ti = 'G';
                Tj = 'F';
            end
            warning('Generic system: %c[%d] and %c[%d] are equal.\n',Ti,iT,Tj,jT)
        end         
    end
end

%%% Find intersection points
for i = 1:nT-2
    for j = i+1:nT-1
        for k = j+1:nT
            I = [i,j,k];
            detNT = det(NT(I,:));
            if abs(detNT) < tol
                continue
            end
            p = NT(I,:) \ cT(I);
            Tmax = max(GT * p(1:2) + cT);
            if nT == 3 || p(3)+tol >= Tmax
                nP = nP + 1;
                TL.pointIndices(:,nP) = I;
                TL.points(:,nP) = p(1:2);
                for m = 1:3
                    iP = I(mod(m,3)+1);
                    jP = I(mod(m+1,3)+1);                    
                    if TL.numPoints(iP,jP) == 1 && all(abs(p(1:2) - TL.points(:,TL.pointMap(iP,jP,1))) < tol)
                        % The point has already been added for this line.
                        continue
                    elseif TL.numPoints(iP,jP) == 2
                        % A line segment cannot have more than 2 endpoints
                        continue
                    end
                    TL.numPoints(iP,jP) = TL.numPoints(iP,jP) + 1;
                    TL.numPoints(jP,iP) = TL.numPoints(iP,jP);
                    TL.pointMap(iP,jP,TL.numPoints(iP,jP)) = nP;
                    TL.pointMap(jP,iP,TL.numPoints(iP,jP)) = nP;
                end
            end
        end
    end
end
TL.points = TL.points(:,1:nP);

%%% Find vectors
for i = 1:nT-1
    for j = i+1:nT
        iP = TL.pointMap(i,j,1);
        if ~iP
            continue
        end
        I = TL.pointIndices(:,iP);
        detNT = det(NT(I,:));
        for m = 1:3
            nV = nV + 1;
            iV = I(mod(m,3)+1);
            jV = I(mod(m+1,3)+1);
            r = detNT * cross(NT(iV,:), NT(jV,:));
            TL.vectors(:,nV) = 1/norm(r(1:2)) * r(1:2);
            TL.vectorMap(iV,jV) = nV;
            TL.vectorMap(jV,iV) = nV;
        end
        
        %%% Line Classification
        v = TL.vectors(:,TL.vectorMap(i,j));
        if (i <= nF && j <= nF) || (i > nF && j > nF)
            % Dominating directions are parallel
            if sign(dT(i)) == sign(dT(j))
                TL.lineClassification(i,j) = STRAIGHT_CROSSING;
            else
                TL.lineClassification(i,j) = NULLCLINE;
                if i <= nF
                    ks = nF+1;
                    ke = nT;
                else
                    ks = 1;
                    ke = nF;
                end
                % Find tropical points and equilibria
                for k = ks:ke-1
                    for l = ks+1:ke
                        A = [diff(GT([i,j],:)); diff(GT([k,l],:))];
                        b = [diff(cT([j,i])); diff(cT([l,k]))];
                        if abs(det(A)) < tol
                            continue
                        end
                        p = A \ b;
                        if TL.numPoints(i,j) == 2
                            Q = TL.points(:,TL.pointMap(i,j,:));
                            qmin = min(Q,[],2);
                            qmax = max(Q,[],2);
                            pqmin = p - qmin;
                            pqmax = p - qmax;
                            if (pqmin(1) < tol && pqmax(1) < tol) ||...
                                    (pqmin(2) < tol && pqmax(2) < tol) ||...
                                    (pqmin(1) > -tol && pqmax(1) > -tol) ||...
                                    (pqmin(2) > -tol && pqmax(2) > -tol)
                                % Point is outside the lineseqment
                                continue
                            end
                        elseif TL.numPoints(i,j) == 1
                            q = TL.points(:,TL.pointMap(i,j,1));
                            rp = p-q;
                            rq = TL.vectors(:,TL.vectorMap(i,j));
                            if ~all(rp.*rq > 0)
                                % Point is outside the half line
                                continue
                            end
                        else
                            continue
                        end
                        TmaxK = max(GT(ks:ke,:) * p + cT(ks:ke));
                        TmaxI = max(GT * p + cT);
                        pk3 = GT(k,:) * p + cT(k);
                        pi3 = GT(i,:) * p + cT(i);
                        if length(ks:ke) == 2 || ...
                                (pk3 + tol >= TmaxK && pi3 + tol >= TmaxI)
                            nE = nE + 1;
                            TL.tropicalPoints(:,nE) = p;
                            TL.tropicalPointIndices(:,nE) = [i j k l];
                            TL.numTropicalPoints(i,j) = TL.numTropicalPoints(i,j) + 1;
                            TL.tropicalPointMap(i,j,TL.numTropicalPoints(i,j)) = nE;
                            TL.distanceToTropicalPoint(i,j,TL.numTropicalPoints(i,j)) = norm(p - TL.points(:,TL.pointMap(i,j,1)));
                            TL.isEquilibrium(nE) = dT(k) * dT(l) < 0;
                        end
                    end
                end
                distances = TL.distanceToTropicalPoint(i,j,1:TL.numTropicalPoints(i,j));
                [~,id] = sort(distances,'ComparisonMethod','abs');
                TL.tropicalPointMap(i,j,1:TL.numTropicalPoints(i,j)) = TL.tropicalPointMap(i,j,id);
                TL.distanceToTropicalPoint(i,j,1:TL.numTropicalPoints(i,j)) = distances(id);
            end
        else
            if i <= nF % && j > nF
                sigma = -dT(i)*v(2) * dT(j)*v(1);
            else
                sigma = -dT(j)*v(2) * dT(i)*v(1);
            end
            if sigma > 0
                TL.lineClassification(i,j) = BENT_CROSSING;
            elseif sigma < 0
                TL.lineClassification(i,j) = SLIDING;
            else
                TL.lineClassification(i,j) = CROSS_N_SLIDE;
            end
        end
    end
end
TL.tropicalPoints = TL.tropicalPoints(:,1:nE);
TL.tropicalPointIndices = TL.tropicalPointIndices(:,1:nE);
TL.vectors = TL.vectors(:,1:nV);
end

