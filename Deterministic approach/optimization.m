%% DETERMINISTIC OPTIMIZATION PROBLEM USING fmincon %%
%% Created by Swaminath VENKATESWARAN (14.07.2018)

function X = optimization(P0)

%% Function Input(s): Random guess points provided by the used (Ex: 10*rand(8,2))
%% Function Ouput(s): Optimal positions of guess points to cover short distances and avoid obstacles


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

        f1= norm(A-P0(1,:));
        f2= norm(B-P0(end,:));
        f3= zeros();

        for i= 1:length(P0)-1
            f3(i)= norm(P0(i,:)-P0(i+1,:));
        end

        f= f1+f2+sum(f3);

    end

%% CONSTRAINT SUB-FUNCTION %%
    function [c,ceq,data]= const(P0)

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
    
    x0=50;
    y0=15;
    largeur=850;
    hauteur=700;
    set(gcf,'units','points','position',[x0,y0,largeur,hauteur])
    title('Optimization of a trajectory inside a 10 x 10 room'); 
    xlabel('x (mm)')
    ylabel('y (mm)')
    xlim([-5 15]);
    ylim([-5 15]);
    grid on;
    
    h1=viscircles(C,r,'Color','r');
    set(gca,'FontSize',20,'FontName','Times New Roman','FontWeight','Bold')
    hold on;
    h2=rectangle('Position',[4 0 2 4],'Facecolor',[0 0.45 0.74],'EdgeColor','w');
    h3=rectangle('Position',[4 6 2 4],'Facecolor',[0 0.45 0.74],'EdgeColor','w');
    
    h4=text(A(1,1),A(1,2), 'A' ,'VerticalAlignment','bottom', ...
                              'HorizontalAlignment','right',....
                              'FontSize', 18,'FontName','Times New Roman','FontWeight','Bold'); 
                          
    h5=text(B(1,1),B(1,2), 'B' ,'VerticalAlignment','top', ...
                              'HorizontalAlignment','left',....
                              'FontSize', 18,'FontName','Times New Roman','FontWeight','Bold'); 
                          
    h6= plot(A(1,1),A(1,2),'kx','LineWidth',2);
                          
    h7=plot(B(1,1),B(1,2), 'kx' ,'LineWidth',2); 
                          
    h8=plot([0,0],[0,10],'Color','k','LineWidth',2);
    h9=plot([0,10],[10,10],'Color','k','LineWidth',2);
    h10=plot([10,10],[10,0],'Color','k','LineWidth',2);
    h11=plot([10,0],[0,0],'Color','k','LineWidth',2);
    
    disc= 0:0.01:1;
    
    pos_nonopt= zeros(length(disc),2); pos_opt=zeros(length(disc),2);
    
    Dist_nonopt=0; Dist_opt=0;
    
    for o= 1 : length(M)-1
        
        for m= 1 : length(disc)
            
            pos_nonopt(m,:) = M1(o,:) + (disc(m) * (M1(o+1,:) - M1(o,:)));
            pos_opt(m,:)= M(o,:) + disc(m) * (M(o+1,:) - M(o,:));
                      
        end
        
        for m= 1 : length(disc)-1
            Dist_nonopt= Dist_nonopt + norm(pos_nonopt(m,:) - pos_nonopt(m+1,:));
            Dist_opt= Dist_opt + norm(pos_opt(m,:) - pos_opt(m+1,:));
            h12=plot([pos_nonopt(m,1),pos_nonopt(m+1,1)],[pos_nonopt(m,2),pos_nonopt(m+1,2)],'r--','LineWidth',0.1);
            h13=plot([pos_opt(m,1),pos_opt(m+1,1)],[pos_opt(m,2),pos_opt(m+1,2)],'Linewidth',2,'Color','b');
            h14=text(-2.5,14, 'Distance travelled from A to B :','FontSize', 18,'FontName','Times New Roman','FontWeight','Bold');
            h15=text(-2.5,12.5, ['- User defined path= ' num2str(Dist_nonopt) ' mm'],...
						'FontSize', 18,'FontName','Times New Roman','FontWeight','Bold'); 
            h16=text(-2.5,11, ['- Optimal path= ' num2str(Dist_opt) ' mm'],...
					    'FontSize', 18,'FontName','Times New Roman','FontWeight','Bold'); 
            drawnow;
            pause(0.01); 
            delete(h14);
            delete(h15);
            delete(h16);            
        end
        
        if (o< length(M)-1)
            
            h17= plot(pos_nonopt(end,1),pos_nonopt(end,2),'r*','LineWidth',2);
            h18= plot(pos_opt(end,1),pos_opt(end,2),'b*','LineWidth',2);
            
        end
                                                         
    end    
    
    legend([h17 h18],'User-defined points','Relocated points after optimization','Location','Southwest')
    
       
end
