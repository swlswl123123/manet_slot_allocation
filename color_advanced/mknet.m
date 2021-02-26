function [X,Y,D,f] = mknet(N,hop,width,height)
%myFun - Description
%
% Syntax: output = mknet(N,hop,width,height)
%
% Long description

while 1

    X = rand(N,1)*width - width/2;
    Y = rand(N,1)*height - height/2;
    D = sqrt((((repmat(X,[1,N]) - repmat(X,[1,N])').^2) + ((repmat(Y,[1,N]) - repmat(Y,[1,N])').^2)));
    D(D==0)=1;
    D(D>1)=0;
    D(D>0)=1;
    [p,q,v]=find(D);
    b=sparse(p,q,v,N,N);
    hop_count=zeros(N,N);
    for i=1:N-1
        for j=i+1:N
            [hop_count(i,j),~,~]=graphshortestpath(b,i,j,'Directed',false);
            hop_count(j,i) = hop_count(i,j);
        end
    end

    hop_count(hop_count==Inf)=-1;

    hop_max=max(max(hop_count));
    hop_min=min(min(hop_count));

    % disp(['HOP = ',num2str(hop_max)]);
    if(hop_max == hop && hop_min~=-1)

        if max(sum(D))>20
            continue;
        end
        break;
    end

    % figure(1);
    % scatter();

    % figure(2);
end

D = uint8(D);

f = figure;
f.Units = 'centimeters';
f.PaperSize=[8,6];
scatter(X,Y,'filled','r');
axis equal;

names = cell(N,1);
for i =1:N
names{i} = num2str(i);
end

text(X,Y,names,'FontSize',12);

hold on;
for i=1:N-1
    for j = i+1:N
        if D(i,j)==1
        line([X(i),X(j)]',[Y(i),Y(j)]','LineStyle',':');
        end
    end
end
hold off;

xlabel('x/D_{hop}');
ylabel('y/D_{hop}');

% G=graph(triu(D==1,1),'upper');
% figure(2);
% h2 = plot(G);

% for i=1:N
%     labelnode(h2,i,num2str(i));
% end

    
end