%% DETERMINISTIC OPTIMIZATION PROBLEM USING fmincon %%
%% Created by Swaminath VENKATESWARAN (14.07.2018)

function X = optimization(P0)

%% Function Input(s): P0- Random guess points provided by the user (Ex: 10*rand(8,2))
%% Function Ouput(s): X- Optimal positions of guess points to cover short distance 
                      ...from A to B and avoid obstacles (Results of "fmincon")


%% INITIALIZATION

    global A;
    global B;
    global C;
    global r;
    global D;
    
    A= [1 9];
    B= [9 1];
    C= [2 3;2 8;8 2;8 9];
    r=[0.7;0.7;0.6;0.7];
   

    options=optimset('Display','iter','Tolx',1e-10,'MaxIter',20000,...
                'MaxFunEvals',40000);
            
    tic;
    [lb,ub,D]= bounds(P0);
    X= fmincon(@(P0)obj(P0),P0,[],[],[],[],lb,ub,@(P0)const(P0),options);
    toc;

    %% OBJECTIVE SUB-FUNCTION %%
    
    function f= obj(P0)
        
    %% Function Input(s): P0- Random guess points provided by the user (Passed from optimization.m)
    %% Function Ouput(s): f- Distance travelled from A to B through the intermediate points
        
    %%
        f1= norm(A-P0(1,:));
        f2= norm(B-P0(end,:));
        f3= zeros();

        for i= 1:length(P0)-1
            f3(i)= norm(P0(i,:)-P0(i+1,:));
        end

        f= f1+f2+sum(f3);

    end

%% CONSTRAINT SUB-FUNCTION %%

    function [c,ceq]= const(P0)
        
    %% Function Input(s): P0- Random guess points provided by the user (Passed from optimization.m)
    %% Function Ouput(s): c- Values of inequality constraints
                      ... ceq- Values of equality constraints
                          
    %%
        seg= [A;P0;B];
        n=0:0.01:1;
        c= zeros();
        data= struct();
        for i= length(seg)-1:-1:1
            for j= length(C):-1:1
                for k= length(n):-1:1
                    data(i).disc(k,:)= seg(i,:)+ (n(k)*(seg(i+1,:)-seg(i,:)));
                    data(i).dist(j,k)= norm(data(i).disc(k,:)-C(j,:));
                end
                data(i).minid(j,:)= min(data(i).dist(j,:));
                data(i).g(j,:)=r(j)-data(i).minid(j,:);
                c(j,i)= data(i).g(j,:);
            end  
        end
        ceq=[];
    end

%% BOUNDS SUB-FUNCTION%%

    function [lb,ub,D]= bounds(P0)
        
    %% Function Input(s): P0- Random guess points provided by the user (Passed from optimization.m)
    %% Function Ouput(s): lb- Lower bound for a particular guess point
                      ... ub- Upper bound for a particular guess point
                      ... D- An array to fix lb and ub zone for a particular guess point
                          
    %%
        lb= zeros();
        ub= zeros();
        D= [0 0; 4 10; 4 4; 4 6; 6 4; 6 6; 6 0; 10 10];

        for i = 1:length(P0)

            for j= 1:2

                if i<= round(length(P0)/2)-2

                    lb(i,j)= D(1,j);
                    ub(i,j)= D(2,j);

                elseif i == round(length(P0)/2)-1

                    lb(i,j)= D(3,j);
                    ub(i,j)= D(4,j);

                elseif i == round(length(P0)/2)

                    lb(i,j)= D(5,j);
                    ub(i,j)= D(6,j);

                elseif i > round(length(P0)/2)

                    lb(i,j)= D(7,j);
                    ub(i,j)= D(8,j);

                end

            end

        end
    end

%% PLOTTING %%

    M= [A;X;B];
    M1= [A;P0;B];
    
    Results = struct();
    disc= 0:0.01:1;
    
    for o = 1 : length(M)-1
        
        for m= 1 : length(disc)
            
            Results(o).pos_nonopt(m,:) = M1(o,:) + (disc(m) * (M1(o+1,:) - M1(o,:)));
            Results(o).pos_opt(m,:) = M(o,:) + (disc(m) * (M(o+1,:) - M(o,:)));
            
        end
        
    end
    
    dist_nonopt=0; dist_opt=0;
    x0=50;
    y0=15;
    largeur=1250;
    hauteur=700;
    set(gcf,'units','points','position',[x0,y0,largeur,hauteur])
    
    for o = 1 : length(M)-1
        
        for m= 1 : length(disc)-1
            
            subplot(1,2,1)
            xlabel('x (mm)')
            ylabel('y (mm)')
            xlim([-5 15]);
            ylim([-5 15]);
            grid on;
            [H1,dist_nonopt]= roomdesc(Results(o).pos_nonopt(m,:),Results(o).pos_nonopt(m+1,:),dist_nonopt,1);
            subplot(1,2,2)
            xlabel('x (mm)')
            ylabel('y (mm)')
            xlim([-5 15]);
            ylim([-5 15]);
            grid on;
            [H2,dist_opt]= roomdesc(Results(o).pos_opt(m,:),Results(o).pos_opt(m+1,:),dist_opt,2);
            drawnow;
            F(o,m)= getframe(gcf);
            pause(0.001);
            delete(H1);
            delete(H2);
            
        end
        
        if (o< length(M)-1)
                
                subplot(1,2,1)
                plot(Results(o).pos_nonopt(end,1),Results(o).pos_nonopt(end,2),'r*','LineWidth',2);
                subplot(1,2,2)
                plot(Results(o).pos_opt(end,1),Results(o).pos_opt(end,2),'b*','LineWidth',2);
                
        end
              
    end  
    
    video = VideoWriter('Optimization','MPEG-4');
    open(video);
    [s,t] = size(F);
    for u = 1:s    
        for v = 1:t
            writeVideo(video,F(u,v));
        end
    end
    close(video)
  
%% SUB-FUNCTION TO PLOT TEST SPACE WITH OBSTACLES%%

    function [H,dist]= roomdesc(X,Y,dist,g)
        
    %% Function Input(s): X- Position of a particular point on a particular segment
                      ... Y- Position of successive point of X on a particular segment
                      ... dist- Scalar to calculate distance travelled from X to Y
                      ... g- A scalar condition to calculate user-defined and optimal distance seperately
    %% Function Ouput(s): H- Vector to save plots
                      ... dist- Calculated distance travelled from X to Y for incrementation in loop
                          
    %%
        
        viscircles(C,r,'Color','r');
        set(gca,'FontSize',20,'FontName','Times New Roman','FontWeight','Bold')
        hold on;
        rectangle('Position',[4 0 2 4],'Facecolor',[0 0.45 0.74],'EdgeColor','w');
        rectangle('Position',[4 6 2 4],'Facecolor',[0 0.45 0.74],'EdgeColor','w');
        H(1)= text(A(1,1),A(1,2), 'A' ,'VerticalAlignment','bottom', ...
                              'HorizontalAlignment','right',....
                              'FontSize', 18,'FontName','Times New Roman','FontWeight','Bold');
        H(2)= text(B(1,1),B(1,2), 'B' ,'VerticalAlignment','top', ...
                              'HorizontalAlignment','left',....
                              'FontSize', 18,'FontName','Times New Roman','FontWeight','Bold'); 
        H(3)= plot(A(1,1),A(1,2),'kx','LineWidth',2);
        H(4)= plot(B(1,1),B(1,2), 'kx' ,'LineWidth',2); 
        H(5)= plot([0,0],[0,10],'Color','k','LineWidth',2);
        H(6)= plot([0,10],[10,10],'Color','k','LineWidth',2);
        H(7)= plot([10,10],[10,0],'Color','k','LineWidth',2);
        H(8)= plot([10,0],[0,0],'Color','k','LineWidth',2);
        
        dist= dist + norm(X - Y);
                
        if (g==1)
            plot([X(1,1),Y(1,1)],[X(1,2),Y(1,2)],'r','LineWidth',2);
            H(9)= text(-2.5,12.5, ['- User defined path (random)= ' num2str(dist) ' mm'],...
                                'FontSize', 16,'FontName','Times New Roman','FontWeight','Bold'); 

        elseif (g==2)
            plot([X(1,1),Y(1,1)],[X(1,2),Y(1,2)],'b','LineWidth',2);
            H(9)= text(-2.5,12.5, ['- Optimal path= ' num2str(dist) ' mm'],...
                                'FontSize', 16,'FontName','Times New Roman','FontWeight','Bold'); 
        end
            H(10)=text(-2.5,14, 'Distance travelled from A to B :','FontSize', 16,...
                                'FontName','Times New Roman','FontWeight','Bold');

    end   
    
end
