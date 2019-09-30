function plot_signal( X_orig, X, X_rec )

[n,T] = size(X_orig);
% this script works only for height = width = n = 1

% if image_rec is [], then we plot only original and with the noise
if isempty(X_rec)
    toplot = 2;
else
    toplot = 3;
end

figure
for i=1:n
    subplot(n,1,i);
    hold on

    plot(1:T,X(i,:),'b','Color',[0.8,0.8,1.0],'LineWidth',1.0)
    plot(1:T,X_orig(i,:),'g','Color',[0,0.6,0],'LineWidth',2.0)
    
    xlabel('$t$','Interpreter','latex')
    ylabel(['$x_{' num2str(i) '}(t)$'],'Interpreter','latex')

    if toplot < 3
        if i==1
            legend('data','true')
        end
    else
        plot(1:T,X_rec(i,:),'r','LineWidth',2.0)
        if i==1
            legend('data','true','reconstructed')
        end
    end
    
    hold off
end

end

