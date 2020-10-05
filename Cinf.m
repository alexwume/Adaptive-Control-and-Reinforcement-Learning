function Cinf
% computation of an invariant set
%
% M. Herceg, Automatic Control Laboratory, ETH Zurich 2014
%
% requires MPT3

close all
clear all
%% printing parameters
label_font_size = 14;
tick_font_size = 10;
line_width = 0.8;
axeswidth=0.2;
figure_name = ['.',filesep,'figures',filesep,'invariantSet_ex_12_4'];

%% core code
% matrices of LTI dynamics 
% x(k+1) = A*x(k) + B*u(k)
T = 0.1;
A = [1, T; 0,1];
B = [(T^2) / 2; T];


% create model in MPT3 interface
model = LTISystem('A',A,'B',B);

% constraints on inputs and states
x_bar = 5;
u_bar = 1;
model.u.min = -u_bar;
model.u.max = u_bar;
model.x.min = [-x_bar;-x_bar];
model.x.max = [ x_bar; x_bar];

% constraint sets represented as polyhedra
X = Polyhedron('lb',model.x.min,'ub',model.x.max);
U = Polyhedron('lb',model.u.min,'ub',model.u.max);

% compute the invariant set including the intermediate steps
maxIterations = 50;
Piter = []; % save intermediate set during computation
X0 = X; % initial set constraint
for i = 1:maxIterations
    % compute backward reachable set
    P = model.reachableSet('X', X0, 'U', U, ...
        'direction', 'backward');
    % intersect with the state constraints
    P = P.intersect(X0).minHRep();
    Piter = [Piter, P];
    if P==X
        break
    else
        X0 = P;
    end
end

% direct command for computing the invariant set (but here we want to show
% intermediate sets)
% Cinf = model.invariantSet(); 

% the invariant set
Cinf = P;

figure
% plot the constraint set
plot(X,'color',[0.7 0.7 0.7],'linewidth',line_width);
hold on
grid on
% plot the intermediate sets
hp = plot(Piter,'color',[0.7 0.7 0.7]);
% plot the invariant set
plot(Cinf,'color',[0.5 0.5 0.5],'linewidth',line_width)
axis([model.x.min(1),model.x.max(1),model.x.min(2),model.x.max(2)])

set(gca,'LineWidth',axeswidth)
set(gca,'FontSize', tick_font_size);
xt = transpose(-10:5:10);
yt = transpose(-10:plot(1:N, squeeze(x_star(1,1,1:N) - x_noisy(1,1,1:N)));5:10);
set(gca,'XTick',xt);
set(gca,'YTick',yt);
set(gca,'YTickLabel',num2str(xt));
set(gca,'XTickLabel',num2str(yt));

hx1 = xlabel('$x_1$', 'Interpreter','latex');
set(hx1, 'FontSize', label_font_size);
hy1 = ylabel('$x_2$','Interpreter','latex');
set(hy1, 'FontSize', label_font_size);

ht1=text(7,-8,'$\mathcal{X}$','Interpreter','latex');
set(ht1, 'FontSize', label_font_size);
ht2=text(-1,-0.5,'$\mathcal{C}_{\infty}$','Interpreter','latex');
set(ht2, 'FontSize', label_font_size);

border = 1;
axis([model.x.min(1)-border,model.x.max(1)+border,model.x.min(2)-border,model.x.max(2)+border])

%% save book figures
% generate FIG file
% disp('Hit any key to print the figure and continue.');
% pause
% saveas(gcf, [figure_name,'_matlab'], 'fig');

% print 
% laprint(gcf, [figure_name,'_tex'],'scalefonts','off');

end