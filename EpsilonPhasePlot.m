function [h,w1] = EpsilonPhasePlot(F,G,dF,dG,w0,epsilon,tspan,maxseconds,color)

% Initialize variables
lineWidth = 2;
ax = gca;
hold(ax,'on')

if ~exist('epsilon','var')
    epsilon = 0.2;
end
if ~exist('tspan','var')
    tspan = [0 1e4];
end
if ~exist('maxseconds','var')
    maxseconds = 1;
end
if ~exist('color','var')
    color = [0,0.6,0];
end

% % ODE system
% udot = @(w) dF * exp(1/epsilon * F * [1;w]);
% vdot = @(w) dG * exp(1/epsilon * G * [1;w]);
% wdot = @(t,w) 1/epsilon* [udot(w); vdot(w)];
% % Jacobian
% J11 = @(w) F(:,2)'.*dF * exp(1/epsilon * F * [1;w]);
% J12 = @(w) F(:,3)'.*dF * exp(1/epsilon * F * [1;w]);
% J21 = @(w) G(:,2)'.*dG * exp(1/epsilon * G * [1;w]);
% J22 = @(w) G(:,3)'.*dG * exp(1/epsilon * G * [1;w]);
% J = @(t,w) 1/epsilon^2*[J11(w), J12(w); J21(w), J22(w)];

% ODE system
udot = @(w) epsilon * dF * exp(1/epsilon * F * [1;w]);
vdot = @(w) epsilon * dG * exp(1/epsilon * G * [1;w]);
wdot = @(t,w) [udot(w); vdot(w)];
% Jacobian
J11 = @(w) F(:,2)'.*dF * exp(1/epsilon * F * [1;w]);
J12 = @(w) F(:,3)'.*dF * exp(1/epsilon * F * [1;w]);
J21 = @(w) G(:,2)'.*dG * exp(1/epsilon * G * [1;w]);
J22 = @(w) G(:,3)'.*dG * exp(1/epsilon * G * [1;w]);
J = @(t,w) [J11(w), J12(w); J21(w), J22(w)];

% Solve and plot
initialTime = cputime;
options = odeset('Jacobian',J,...
    'Events',@(t,w)event(w,ax.XLim,ax.YLim,initialTime,maxseconds));
% options = odeset('Events',@(t,w)event(w,ax.XLim,ax.YLim,initialTime,maxseconds));
[~,wID1] = lastwarn;
warning('off','all')
[~,w] = ode23s(wdot,tspan,w0,options);
warning('on','all')
[~,wID2] = lastwarn;
if wID1 ~= wID2
    warning(lastwarn)
end
handle = plot(w(:,1),w(:,2),'Color',color,'LineWidth',lineWidth);
plot(w0(1),w0(2),'ok','MarkerSize',6,'MarkerFaceColor',color,'LineWidth',1)
if nargout >= 1
    h = handle;
end
if nargout >= 2
    w1 = w(end,:); 
end
hold(ax,'off')
end

function [value, isterminal, direction] = event(w,ulim,vlim,initialTime,maximumTime)
dt = cputime - initialTime;
timeleft = maximumTime - dt;
du = (ulim - w(1)) .* [-1 1];
dv = (vlim - w(2)) .* [-1 1];
wmin = zeros(1,2);
wmin(1) = min(du);
wmin(2) = min(dv);
value = min(timeleft,min(wmin));
isterminal = 1;
direction = 0;
end

