function X = optimization(P0)

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
viscircles(C,r,'Color','r');
[x,y,z]=cylinder(0.7);
x0=C(1,1);y0=C(1,2);
x=x+x0;
y=y+y0;
hold on;
surf(x,y,z);
rotate3d;
rectangle('Position',[4 0 2 4],'Facecolor',[0 0.45 0.74],'EdgeColor','w')
rectangle('Position',[4 6 2 4],'Facecolor',[0 0.45 0.74],'EdgeColor','w')
plot(M(:,1),M(:,2),'Linewidth',2);
plot(M1(:,1),M1(:,2),'--');
xlim([0 10]);
ylim([0 10]);
grid on;

end