function plot_signal_clustering_compare( Gamma1, Gamma2 )

[K,T] = size(Gamma1);

for k=1:K
    subplot(K,1,k);
    hold on

    plot(1:T,Gamma1(k,:),'b')
    plot(1:T,Gamma2(k,:),'r')

    xlabel('$t$','Interpreter','latex')
    ylabel(['$\Gamma_{' num2str(k) ',:}$'],'Interpreter','latex')
    
    axis([1 T -0.2 1.2])
    hold off
end

end

