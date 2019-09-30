function plot_signal_clustering( Gamma, Theta )

[K,T] = size(Gamma);

figure
for k=1:K
    subplot(K,1,k);
    hold on
    title(['$\theta_' num2str(k) ' = [' vec2str(Theta(:,k)) ']$'], 'Interpreter', 'latex');

    plot(1:T,Gamma(k,:),'b')

    xlabel('$t$','Interpreter','latex')
    ylabel(['$\gamma_{' num2str(k) '}(t)$'],'Interpreter','latex')
    
    axis([1 T -0.2 1.2])
    hold off
end

end

function out = vec2str(v)
    out = '';
    for i = 1:length(v)
        if i > 1
            out = [out ','];
        end
        out = [out num2str(v(i))]; 
    end
end